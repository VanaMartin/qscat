"""Host-shared search runtime for independent main corpora and lightweight MCP mounts."""

from __future__ import annotations

import argparse
import asyncio
import fcntl
import json
import logging
import math
import os
import re
import signal
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from tools.agent_search import __main__ as indexer
from tools.agent_search.git_source import MainSource
from tools.agent_search.service_client import (
    IDLE_TIMEOUT,
    LEASE,
    MAX_MESSAGE,
    REQUEST_TIMEOUT,
    implementation_id,
    rpc,
    runtime_directory,
    socket_path,
)


@dataclass
class Corpus:
    table: str
    owner: dict
    mounts: dict = field(default_factory=dict)
    first_attempt: asyncio.Event = field(default_factory=asyncio.Event)
    wake: asyncio.Event = field(default_factory=asyncio.Event)
    task: asyncio.Task | None = None
    error: str | None = None

    @property
    def interval(self) -> float:
        return min(mount["interval"] for mount in self.mounts.values())


class SearchService:
    def __init__(self, db_path: Path, *, idle_timeout: float = IDLE_TIMEOUT):
        self.db_path = db_path
        self.idle_timeout = idle_timeout
        self.implementation = implementation_id()
        self.corpora: dict[str, Corpus] = {}
        self.active_requests = 0
        self.last_used = time.monotonic()
        self.stop = asyncio.Event()

    def health(self) -> dict:
        return {
            "pid": os.getpid(),
            "implementation": self.implementation,
            "socket": str(socket_path(self.db_path)),
            "db_uri": str(self.db_path),
            "loaded_models": len(indexer._EMBEDDINGS),
            "proxy_pids": sorted(
                {
                    mount["proxy_pid"]
                    for corpus in self.corpora.values()
                    for mount in corpus.mounts.values()
                }
            ),
            "corpora": {
                table: {"mounts": len(corpus.mounts), "maintainer_error": corpus.error}
                for table, corpus in self.corpora.items()
            },
        }

    async def resolve(self, mount: dict) -> MainSource:
        source = await asyncio.to_thread(
            MainSource.discover,
            Path(mount["root"]),
            mount["remote"],
            mount["profile"],
            corpus=mount["corpus"],
        )
        if source.owner != mount["owner"]:
            raise ValueError("Configured upstream changed since mount")
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", mount["table"]):
            raise ValueError("Table name must be an identifier")
        interval = mount["interval"]
        if not math.isfinite(interval) or interval <= 0:
            raise ValueError("Main refresh interval must be positive")
        return source

    async def register(self, mount: dict) -> tuple[MainSource, Corpus]:
        source = await self.resolve(mount)
        table = mount["table"]
        corpus = self.corpora.setdefault(table, Corpus(table, source.owner))
        if corpus.owner != source.owner:
            raise ValueError("Corpus belongs to a different upstream or profile")
        new = mount["client_id"] not in corpus.mounts
        corpus.mounts[mount["client_id"]] = {
            **mount,
            "source": source,
            "expires": time.monotonic() + LEASE,
        }
        if corpus.task is None or corpus.task.done():
            corpus.first_attempt.clear()
            corpus.task = asyncio.create_task(self.maintain(corpus))
        elif new:
            corpus.wake.set()
        return source, corpus

    async def maintain(self, corpus: Corpus) -> None:
        while corpus.mounts:
            corpus.wake.clear()
            started = time.monotonic()
            interval = corpus.interval
            # A removed worktree or changed remote must not strand other mounts.
            sources = {mount["source"] for mount in corpus.mounts.values()}
            for attempt, source in enumerate(sorted(sources, key=lambda source: str(source.root))):
                try:
                    await asyncio.to_thread(
                        indexer.ensure_main,
                        source,
                        self.db_path,
                        corpus.table,
                        interval=interval,
                        force=attempt > 0,
                    )
                    corpus.error = None
                    break
                except Exception as error:
                    corpus.error = str(error)
                    logging.getLogger(__name__).exception("Main-index refresh failed")
            corpus.first_attempt.set()
            if not corpus.mounts:
                break
            try:
                await asyncio.wait_for(
                    corpus.wake.wait(), timeout=max(0.05, interval - (time.monotonic() - started))
                )
            except TimeoutError:
                pass

    async def dispatch(self, request: dict) -> dict:
        operation = request["operation"]
        if operation == "health":
            return self.health()
        if operation == "stop":
            self.stop.set()
            return {"stopping": True, "pid": os.getpid()}
        if request.get("implementation") != self.implementation:
            raise ValueError(
                "Shared search service runs different code/dependencies; stop it first"
            )
        if operation == "unmount":
            for corpus in self.corpora.values():
                if corpus.mounts.pop(request["client_id"], None):
                    corpus.wake.set()
            return {"unmounted": True}
        if operation == "command":
            # One-shot CLI reads stay read-only; plans must not start a publisher.
            source = await self.resolve(request["mount"])
            table = request["mount"]["table"]
            command = request["command"]
            if command == "status":
                return await asyncio.to_thread(indexer.status, source, self.db_path, table)
            if command == "query":
                return await asyncio.to_thread(
                    indexer.query,
                    source,
                    self.db_path,
                    table,
                    request["text"],
                    request["mode"],
                    request["limit"],
                    paper_id=request.get("paper_id"),
                )
            if command in {"sync", "bootstrap", "plan"}:
                return await asyncio.to_thread(
                    indexer.ensure_main,
                    source,
                    self.db_path,
                    table,
                    force=True,
                    plan=command == "plan",
                )
            raise ValueError("Unsupported search command")
        source, corpus = await self.register(request["mount"])
        if operation == "mount":
            return self.health()
        freshness = await asyncio.to_thread(
            indexer.status,
            source,
            self.db_path,
            corpus.table,
            interval=corpus.interval,
        )
        if operation == "details":
            return {
                "name": corpus.table,
                "db_uri": str(self.db_path),
                "num_rows": freshness["chunks"],
                "freshness": freshness,
                "maintainer_error": corpus.error,
                "service": self.health(),
            }
        if operation == "query":
            if (
                request["query_type"] not in {"vector", "fts", "hybrid"}
                or not 1 <= request["top_k"] <= 50
            ):
                raise ValueError("Use vector/fts/hybrid and top_k between 1 and 50")
            if freshness["integrity"] != "complete":
                await corpus.first_attempt.wait()
                if corpus.error:
                    raise RuntimeError(f"Initial main indexing failed: {corpus.error}")
            return await asyncio.to_thread(
                indexer.query,
                source,
                self.db_path,
                corpus.table,
                request["query"],
                request["query_type"],
                request["top_k"],
                interval=corpus.interval,
                paper_id=request.get("paper_id"),
            )
        raise ValueError("Unsupported search service operation")

    async def handle(self, reader, writer) -> None:
        self.active_requests += 1
        try:
            async with asyncio.timeout(REQUEST_TIMEOUT):
                request = json.loads(await reader.readline())
                result = await self.dispatch(request)
            response = {"result": result}
        except Exception as error:
            response = {"error": str(error)}
        try:
            writer.write((json.dumps(response) + "\n").encode())
            await writer.drain()
        except (ConnectionError, OSError):
            pass
        finally:
            writer.close()
            try:
                await writer.wait_closed()
            finally:
                self.active_requests -= 1
                self.last_used = time.monotonic()

    async def reap(self) -> None:
        while not self.stop.is_set():
            await asyncio.sleep(0.5)
            now = time.monotonic()
            for table, corpus in list(self.corpora.items()):
                for client_id, mount in list(corpus.mounts.items()):
                    if mount["expires"] <= now:
                        del corpus.mounts[client_id]
                        corpus.wake.set()
                if not corpus.mounts and corpus.task.done():
                    del self.corpora[table]
            if self.corpora or self.active_requests:
                self.last_used = now
            elif now - self.last_used >= self.idle_timeout:
                self.stop.set()

    async def run(self) -> None:
        path = socket_path(self.db_path)
        path.unlink(missing_ok=True)
        server = await asyncio.start_unix_server(self.handle, path=path, limit=MAX_MESSAGE)
        path.chmod(0o600)
        loop = asyncio.get_running_loop()
        for signum in (signal.SIGTERM, signal.SIGINT):
            loop.add_signal_handler(signum, self.stop.set)
        reaper = asyncio.create_task(self.reap())
        try:
            async with server:
                await self.stop.wait()
        finally:
            reaper.cancel()
            tasks = [corpus.task for corpus in self.corpora.values()]
            for task in tasks:
                task.cancel()
            await asyncio.gather(reaper, *tasks, return_exceptions=True)
            path.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["serve", "status", "stop"])
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")),
    )
    parser.add_argument("--idle-timeout", type=float, default=IDLE_TIMEOUT)
    args = parser.parse_args()
    db_path = args.db.expanduser().resolve()
    if args.command != "serve":
        try:
            result = asyncio.run(rpc(db_path, "health" if args.command == "status" else "stop"))
        except (FileNotFoundError, ConnectionRefusedError):
            result = {"state": "stopped"}
        print(json.dumps(result, indent=2))
        return
    if not math.isfinite(args.idle_timeout) or args.idle_timeout <= 0:
        parser.error("Idle timeout must be positive")
    directory = runtime_directory(db_path)
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    directory.chmod(0o700)
    with (directory / "daemon.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            sys.exit(0)
        asyncio.run(SearchService(db_path, idle_timeout=args.idle_timeout).run())


if __name__ == "__main__":
    main()
