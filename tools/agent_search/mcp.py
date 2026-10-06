"""Serve committed main snapshots, refreshed on mount and at a bounded interval."""

from __future__ import annotations

import argparse
import asyncio
import logging
import math
import os
import re
import time
from contextlib import asynccontextmanager
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from tools.agent_search import __main__ as indexer
from tools.agent_search.git_source import PROFILE, MainSource


def create_server(
    source: MainSource, db_path: Path, table: str, *, interval: float = indexer.CHECK_INTERVAL
) -> FastMCP:
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError("Main refresh interval must be positive")
    task: asyncio.Task | None = None
    first_attempt = asyncio.Event()
    maintainer_error: str | None = None

    async def refresh() -> None:
        nonlocal maintainer_error
        try:
            await asyncio.to_thread(indexer.ensure_main, source, db_path, table, interval=interval)
            maintainer_error = None
        except Exception as error:
            maintainer_error = str(error)
            logging.getLogger(__name__).exception("Main-index refresh failed")
        finally:
            first_attempt.set()

    async def maintain() -> None:
        while True:
            started = time.monotonic()
            await refresh()
            await asyncio.sleep(max(0.05, interval - (time.monotonic() - started)))

    @asynccontextmanager
    async def lifespan(_server):
        nonlocal task
        task = asyncio.create_task(maintain())
        try:
            yield
        finally:
            task.cancel()
            await asyncio.gather(task, return_exceptions=True)

    server = FastMCP("qmodeling-main-code-search", lifespan=lifespan)

    @server.tool()
    async def query_table(query: str, top_k: int = 5, query_type: str = "vector") -> dict:
        """Find leads in a committed upstream-main snapshot using vector/fts/hybrid.

        Results name the indexed commit. Resolve anchors in the task checkout and
        search branch/dirty/new changes locally. A completed older main snapshot
        remains available during refresh or an upstream outage.
        """
        if not 1 <= top_k <= 50 or query_type not in {"vector", "fts", "hybrid"}:
            raise ValueError("Use vector/fts/hybrid and top_k between 1 and 50")
        # There is no prior main snapshot on first mount. Wait for that attempt;
        # later refreshes leave the last completed generation immediately usable.
        details = await asyncio.to_thread(indexer.status, source, db_path, table, interval=interval)
        if details["integrity"] != "complete" and task is not None and not task.done():
            await first_attempt.wait()
            if maintainer_error:
                raise RuntimeError(f"Initial main indexing failed: {maintainer_error}")
        return await asyncio.to_thread(
            indexer.query, source, db_path, table, query, query_type, top_k, interval=interval
        )

    @server.tool()
    async def table_details(table_name: str | None = None, db_uri: str | None = None) -> dict:
        """Inspect the shared main corpus, indexed commit, integrity, and refresh state."""
        if table_name is not None and table_name != table:
            raise ValueError("table_name differs from the mounted main corpus")
        if db_uri is not None and Path(db_uri).expanduser().resolve() != db_path:
            raise ValueError("db_uri differs from the mounted database")
        freshness = await asyncio.to_thread(
            indexer.status, source, db_path, table, interval=interval
        )
        return {
            "name": table,
            "db_uri": str(db_path),
            "num_rows": freshness["chunks"],
            "freshness": freshness,
            "maintainer_error": maintainer_error,
        }

    return server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--remote", default="origin")
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")),
    )
    parser.add_argument("--table", default=os.environ.get("TABLE_NAME"))
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--check-interval", type=float, default=indexer.CHECK_INTERVAL)
    args = parser.parse_args()
    source = MainSource.discover(args.root, args.remote, args.profile)
    table = args.table or source.table
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table):
        parser.error("Table name must be an identifier")
    create_server(source, args.db.expanduser().resolve(), table, interval=args.check_interval).run()


if __name__ == "__main__":
    main()
