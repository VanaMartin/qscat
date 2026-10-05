"""Code-search MCP reader with one-time population in its mount lifecycle."""

from __future__ import annotations

import argparse
import asyncio
import json
import logging
import os
import re
from contextlib import asynccontextmanager
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from tools.agent_search import __main__ as indexer


def create_server(
    root: Path, db_path: Path, table: str, profile: dict, *, profile_path: Path | None = None
) -> FastMCP:
    """Mount immediately, populating an absent/empty owned corpus in the background."""
    task: asyncio.Task | None = None
    bootstrap = {"state": "pending"}

    def current_profile() -> dict:
        return json.loads(profile_path.read_text()) if profile_path is not None else profile

    async def populate() -> None:
        bootstrap["state"] = "running"
        try:
            result = await asyncio.to_thread(
                indexer.reconcile, root, db_path, table, current_profile(), [], only_if_empty=True
            )
            bootstrap.update(state="complete", result=result)
        except Exception as error:
            bootstrap.update(state="failed", error=str(error))
            logging.getLogger(__name__).exception("Code-index mount bootstrap failed")

    @asynccontextmanager
    async def lifespan(_server):
        nonlocal task
        task = asyncio.create_task(populate())
        try:
            yield
        finally:
            if not task.done():
                task.cancel()
            await asyncio.gather(task, return_exceptions=True)

    server = FastMCP("qmodeling-code-search", lifespan=lifespan)

    @server.tool()
    async def query_table(query: str, top_k: int = 5, query_type: str = "vector") -> dict:
        """Find code leads with provenance; wait for mount-time initial population.

        query_type accepts vector, fts, or hybrid. A stale result is a pointer:
        resolve its anchor and read current source before citing or editing it.
        """
        if not 1 <= top_k <= 50:
            raise ValueError("top_k must be between 1 and 50")
        if query_type not in {"vector", "fts", "hybrid"}:
            raise ValueError("query_type must be vector, fts, or hybrid")
        if task is not None:
            await asyncio.shield(task)
        if bootstrap["state"] == "failed":
            raise RuntimeError(f"Initial code indexing failed: {bootstrap['error']}")
        freshness = await asyncio.to_thread(indexer.status, root, db_path, table, current_profile())
        if freshness["state"] not in {"current", "stale"}:
            raise RuntimeError(f"Code index is {freshness['state']}; request checkpoint sync")
        results = await asyncio.to_thread(indexer.query, db_path, table, query, query_type, top_k)
        return {"freshness": freshness, "results": results}

    @server.tool()
    async def table_details(table_name: str | None = None, db_uri: str | None = None) -> dict:
        """Inspect the selected code table, bootstrap progress, and worktree freshness.

        Optional selectors must match this connection. Use a separate MCP mount
        for a different worktree or database.
        """
        if table_name is not None and table_name != table:
            raise ValueError("table_name differs from the mounted worktree table")
        if db_uri is not None and Path(db_uri).expanduser().resolve() != db_path:
            raise ValueError("db_uri differs from the mounted database")

        def inspect() -> dict:
            import lancedb

            freshness = indexer.status(root, db_path, table, current_profile())
            db = lancedb.connect(str(db_path))
            stored = db.open_table(table) if table in db.table_names() else None
            return {
                "name": table,
                "db_uri": str(db_path),
                "root": str(root),
                "num_rows": stored.count_rows() if stored is not None else 0,
                "schema": str(stored.schema) if stored is not None else None,
                "freshness": freshness,
                "bootstrap": dict(bootstrap),
            }

        return await asyncio.to_thread(inspect)

    return server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")),
    )
    parser.add_argument("--table", default=os.environ.get("TABLE_NAME"))
    parser.add_argument("--profile", type=Path, default=Path(".opencode/search/profiles.json"))
    args = parser.parse_args()
    root = Path(indexer.git(args.root, "rev-parse", "--show-toplevel").decode().strip()).resolve()
    db_path = args.db.expanduser().resolve()
    profile_path = (root / args.profile).resolve()
    profile = json.loads(profile_path.read_text())
    table = args.table or indexer.select_table(root, db_path, profile["repository"]["table"])
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table):
        parser.error("Table name must be an identifier")
    create_server(root, db_path, table, profile, profile_path=profile_path).run()


if __name__ == "__main__":
    main()
