"""Expose corpus-bound MCP tools through the host-shared local search runtime."""

from __future__ import annotations

import argparse
import asyncio
import logging
import math
import os
import re
from contextlib import asynccontextmanager
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from tools.agent_search.git_source import PROFILE, MainSource
from tools.agent_search.service_client import CHECK_INTERVAL, HEARTBEAT, SearchClient


def create_server(
    source: MainSource, db_path: Path, table: str, *, interval: float = CHECK_INTERVAL
) -> FastMCP:
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError("Main refresh interval must be positive")
    client = SearchClient(source, db_path, table, interval)
    connection_error: str | None = None

    async def renew() -> None:
        nonlocal connection_error
        while True:
            try:
                await client.call("mount")
                connection_error = None
            except Exception as error:
                connection_error = str(error)
                logging.getLogger(__name__).exception("Shared search connection failed")
            await asyncio.sleep(HEARTBEAT)

    @asynccontextmanager
    async def lifespan(_server):
        task = asyncio.create_task(renew())
        try:
            yield
        finally:
            task.cancel()
            await asyncio.gather(task, return_exceptions=True)
            await client.close()

    name = (
        "qmodeling-main-article-search"
        if source.corpus == "articles"
        else "qmodeling-main-code-search"
    )
    server = FastMCP(name, lifespan=lifespan)

    async def lookup(query, top_k, query_type, paper_id=None):
        if not 1 <= top_k <= 50 or query_type not in {"vector", "fts", "hybrid"}:
            raise ValueError("Use vector/fts/hybrid and top_k between 1 and 50")
        return await client.call(
            "query",
            query=query,
            query_type=query_type,
            top_k=top_k,
            paper_id=paper_id,
        )

    if source.corpus == "articles":

        @server.tool()
        async def query_table(
            query: str,
            top_k: int = 8,
            query_type: str = "hybrid",
            paper_id: str | None = None,
        ) -> dict:
            """Search committed literature notes using vector/fts/hybrid, optionally by paper_id.

            paper_id is the note filename stem. Results carry the note hash, paper
            edition declarations, literal locator clauses and printed-page labels.
            These are tracked-note excerpts, not independently verified full text;
            page clauses may mention other papers. Read the note before citing.
            extraction_page is null because PDFs are not processed by this index.
            """
            return await lookup(query, top_k, query_type, paper_id)
    else:

        @server.tool()
        async def query_table(query: str, top_k: int = 5, query_type: str = "vector") -> dict:
            """Find leads in a committed upstream-main snapshot using vector/fts/hybrid.

            Results name the indexed commit. Resolve anchors in the task checkout and
            search branch/dirty/new changes locally. A completed older main snapshot
            remains available during refresh or an upstream outage.
            """
            return await lookup(query, top_k, query_type)

    @server.tool()
    async def table_details(table_name: str | None = None, db_uri: str | None = None) -> dict:
        """Inspect the mounted corpus, source commit, integrity, FTS and refresh state.

        Explicit table/database selections must match this connection; mismatches
        raise errors rather than inspecting another corpus silently.
        """
        if table_name is not None and table_name != table:
            raise ValueError("table_name differs from the mounted main corpus")
        if db_uri is not None and Path(db_uri).expanduser().resolve() != db_path:
            raise ValueError("db_uri differs from the mounted database")
        details = await client.call("details")
        if connection_error:
            details["connection_error"] = connection_error
        return details

    return server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--remote", default="origin")
    parser.add_argument("--corpus", choices=["repository", "articles"], default="repository")
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")),
    )
    parser.add_argument("--table", default=os.environ.get("TABLE_NAME"))
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--check-interval", type=float, default=CHECK_INTERVAL)
    args = parser.parse_args()
    source = MainSource.discover(args.root, args.remote, args.profile, corpus=args.corpus)
    table = args.table or source.table
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table):
        parser.error("Table name must be an identifier")
    create_server(source, args.db.expanduser().resolve(), table, interval=args.check_interval).run()


if __name__ == "__main__":
    main()
