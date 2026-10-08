"""Lightweight local RPC and singleton startup; importing this never loads a model."""

from __future__ import annotations

import asyncio
import fcntl
import hashlib
import json
import os
import socket
import subprocess
import sys
import time
import uuid
from functools import lru_cache
from importlib.metadata import version
from pathlib import Path

CHECK_INTERVAL = 300
HEARTBEAT = 30
LEASE = 120
IDLE_TIMEOUT = 60
MAX_MESSAGE = 8 * 1024 * 1024
REQUEST_TIMEOUT = 900


@lru_cache(maxsize=1)
def implementation_id() -> str:
    """Share identical code/locked runtimes across clones, not incompatible branches."""
    root = Path(__file__).resolve().parents[2]
    digest = hashlib.sha256()
    for path in sorted((root / "tools/agent_search").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    digest.update((root / ".opencode/search/uv.lock").read_bytes())
    digest.update(str(sys.version_info[:3]).encode())
    for package in (
        "lancedb",
        "mcp",
        "sentence-transformers",
        "torch",
        "transformers",
        "tokenizers",
        "pyarrow",
        "tree-sitter",
        "tree-sitter-rust",
    ):
        digest.update(f"{package}={version(package)}".encode())
    return digest.hexdigest()


def runtime_directory(db_path: Path) -> Path:
    # A short user-owned path also works for deeply nested checkouts/test databases
    # on macOS, whose Unix socket path limit is only 104 bytes.
    key = hashlib.sha256(str(db_path.expanduser().resolve()).encode()).hexdigest()[:20]
    return Path.home() / ".cache/qmodeling-search" / key


def socket_path(db_path: Path) -> Path:
    path = runtime_directory(db_path) / "service.sock"
    if len(os.fsencode(path)) >= 104:
        raise ValueError("Search service socket path exceeds the POSIX path limit")
    return path


def probe(db_path: Path) -> dict:
    with socket.socket(socket.AF_UNIX) as connection:
        connection.settimeout(5)
        connection.connect(str(socket_path(db_path)))
        connection.sendall(b'{"operation":"health"}\n')
        with connection.makefile("rb") as stream:
            return json.loads(stream.readline(MAX_MESSAGE))["result"]


def ensure_running(db_path: Path) -> None:
    directory = runtime_directory(db_path)
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    directory.chmod(0o700)
    with (directory / "startup.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            health = probe(db_path)
        except (FileNotFoundError, ConnectionRefusedError):
            root = Path(__file__).resolve().parents[2]
            with (directory / "service.log").open("w") as log:
                process = subprocess.Popen(
                    [
                        sys.executable,
                        "-m",
                        "tools.agent_search.service",
                        "serve",
                        "--db",
                        str(db_path),
                    ],
                    cwd=root,
                    stdin=subprocess.DEVNULL,
                    stdout=log,
                    stderr=log,
                    start_new_session=True,
                )
            deadline = time.monotonic() + 20
            while True:
                try:
                    health = probe(db_path)
                    break
                except (FileNotFoundError, ConnectionRefusedError):
                    if process.poll() is not None or time.monotonic() >= deadline:
                        raise RuntimeError(
                            f"Search service failed to start; read {directory / 'service.log'}"
                        ) from None
                    time.sleep(0.05)
        if health["implementation"] != implementation_id():
            raise RuntimeError(
                "Shared search service runs different code/dependencies; disconnect older mounts, "
                "then run python -m tools.agent_search.service stop --db " + str(db_path)
            )


async def rpc(db_path: Path, operation: str, **arguments) -> dict:
    async with asyncio.timeout(REQUEST_TIMEOUT):
        reader, writer = await asyncio.open_unix_connection(
            str(socket_path(db_path)), limit=MAX_MESSAGE
        )
        try:
            writer.write((json.dumps({"operation": operation, **arguments}) + "\n").encode())
            await writer.drain()
            response = json.loads(await reader.readline())
            if "error" in response:
                raise RuntimeError(response["error"])
            return response["result"]
        finally:
            writer.close()
            await writer.wait_closed()


class SearchClient:
    """One leased MCP mount; reads recover automatically after a service exit."""

    def __init__(self, source, db_path: Path, table: str, interval: float):
        self.db_path = db_path
        self.mount = {
            "root": str(source.root),
            "remote": source.remote,
            "profile": source.profile_path,
            "corpus": source.corpus,
            "owner": source.owner,
            "table": table,
            "interval": interval,
            "client_id": uuid.uuid4().hex,
            "proxy_pid": os.getpid(),
        }

    async def call(self, operation: str, **arguments) -> dict:
        try:
            return await rpc(
                self.db_path,
                operation,
                implementation=implementation_id(),
                mount=self.mount,
                **arguments,
            )
        except (
            FileNotFoundError,
            ConnectionRefusedError,
            ConnectionResetError,
            BrokenPipeError,
            json.JSONDecodeError,
        ):
            await asyncio.to_thread(ensure_running, self.db_path)
            return await rpc(
                self.db_path,
                operation,
                implementation=implementation_id(),
                mount=self.mount,
                **arguments,
            )

    async def close(self) -> None:
        try:
            await rpc(
                self.db_path,
                "unmount",
                implementation=implementation_id(),
                client_id=self.mount["client_id"],
            )
        except (OSError, RuntimeError, json.JSONDecodeError):
            pass
