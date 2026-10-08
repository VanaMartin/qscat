"""Shared-process, model-reuse and lifecycle checks against real Git/LanceDB/MCP."""

import asyncio
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from contextlib import AsyncExitStack
from pathlib import Path

import pytest

from tests import test_agent_search as fixtures
from tests.test_agent_search import (
    PROFILE,
    commit,
    decoded,
    git,
    parameters,
    publish,
)
from tools.agent_search.git_source import MainSource
from tools.agent_search.service_client import (
    SearchClient,
    probe,
    rpc,
    runtime_directory,
    socket_path,
)

runtime = fixtures.runtime
repos = fixtures.repos
article_repos = fixtures.article_repos
stop_test_search_service = fixtures.stop_test_search_service


def test_proxy_imports_no_model_or_database_runtime():
    pytest.importorskip("mcp")
    result = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "import sys; import tools.agent_search.mcp; "
            "print([name for name in ('torch', 'lancedb', 'sentence_transformers', 'pyarrow') "
            "if name in sys.modules])",
        ]
    )
    assert result.strip() == b"[]"


def test_unchanged_refresh_never_materializes_rows_or_loads_a_model(
    runtime, repos, tmp_path, monkeypatch
):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    first = publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    before = runtime.read_manifest(db_path, source.table)

    def forbidden(*args, **kwargs):
        pytest.fail("Unchanged refresh allocated rows/vectors or loaded a model")

    stored = runtime.open_published(db_path, before)
    monkeypatch.setattr(type(stored), "to_arrow", forbidden)
    monkeypatch.setattr(runtime, "embedding", forbidden)
    unchanged = publish(runtime, worker, db_path)
    assert unchanged["reused_rows"] == first["chunks"]
    assert runtime.read_manifest(db_path, source.table) == before
    (author / "docs/superpowers/plan.md").write_text("# Excluded change\n")
    commit(author)
    git(author, "push", "-q", "origin", "main")
    advanced = publish(runtime, worker, db_path)
    assert advanced["new_embeddings"] == advanced["updated_rows"] == 0
    assert advanced["source_commit"] != before["source_commit"]
    assert (
        runtime.read_manifest(db_path, source.table)["published_table"] == before["published_table"]
    )


def test_explicit_query_embeddings_match_sdk_ranking_and_reuse_one_model(
    runtime, article_repos, tmp_path, monkeypatch
):
    from lancedb.embeddings.sentence_transformers import SentenceTransformerEmbeddings

    _, worker, _ = article_repos
    db_path = tmp_path / "db"
    for corpus in ("repository", "articles"):
        source = MainSource.discover(worker, corpus=corpus)
        runtime.ensure_main(source, db_path, source.table, force=True)
        manifest = runtime.read_manifest(db_path, source.table)
        stored = runtime.open_published(db_path, manifest)
        text = "bilinear pairing" if corpus == "articles" else "symbolic ordering numeric refactor"
        for mode in ("vector", "fts", "hybrid"):
            builder = stored.search(text, query_type=mode)
            if mode != "fts":
                builder = builder.metric("cosine")
            expected = builder.select(["id"]).limit(3).to_list()
            # A published schema's implicit SDK function must never load a second model.
            with monkeypatch.context() as patch:
                patch.setattr(
                    SentenceTransformerEmbeddings,
                    "get_embedding_model",
                    lambda *args: pytest.fail("Implicit SDK model loader used"),
                )
                actual = runtime.query(source, db_path, source.table, text, mode, 3)["results"]
            assert [row["id"] for row in actual] == [row["id"] for row in expected]
            score = (
                "_distance"
                if mode == "vector"
                else "_score"
                if mode == "fts"
                else "_relevance_score"
            )
            assert [row[score] for row in actual] == pytest.approx([row[score] for row in expected])
    assert runtime.embedding(PROFILE)[0] is runtime.embedding(PROFILE)[0]
    assert len(runtime._EMBEDDINGS) == 1


def test_two_corpora_and_two_worktrees_share_one_daemon_and_model(runtime, article_repos, tmp_path):
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client

    author, worker, _ = article_repos
    worktree = tmp_path / "worktree"
    git(worker, "worktree", "add", "--detach", str(worktree), "HEAD")
    db_path = tmp_path / "db"
    implementation = Path(__file__).parents[1]
    for root in (worker, worktree):
        (root / "tools").mkdir()
        shutil.copy(implementation / "tools/__init__.py", root / "tools/__init__.py")
        shutil.copytree(implementation / "tools/agent_search", root / "tools/agent_search")
        shutil.copy(implementation / ".opencode/search/uv.lock", root / ".opencode/search/uv.lock")

    async def exercise():
        async with AsyncExitStack() as stack:
            clients = []
            for root in (worker, worktree):
                for corpus in ("repository", "articles"):
                    read, write = await stack.enter_async_context(
                        stdio_client(parameters(root, db_path, corpus=corpus, cwd=root))
                    )
                    client = await stack.enter_async_context(ClientSession(read, write))
                    await client.initialize()
                    clients.append((client, corpus))
            results = await asyncio.gather(
                *[
                    client.call_tool(
                        "query_table", {"query": "bilinear pairing", "query_type": "hybrid"}
                    )
                    for client, _ in clients
                ]
            )
            assert all(decoded(result)["results"] for result in results)
            details = [
                decoded(await client.call_tool("table_details", {})) for client, _ in clients
            ]
            assert len({detail["service"]["pid"] for detail in details}) == 1
            assert all(detail["service"]["loaded_models"] == 1 for detail in details)
            health = details[-1]["service"]
            assert len(health["proxy_pids"]) == 4
            assert len(health["corpora"]) == 2
            assert all(corpus["mounts"] == 2 for corpus in health["corpora"].values())
            assert socket_path(db_path).stat().st_mode & 0o777 == 0o600
            assert runtime_directory(db_path).stat().st_mode & 0o777 == 0o700
            for detail, (_, corpus) in zip(details, clients, strict=True):
                assert detail["freshness"]["corpus"] == corpus
                assert detail["freshness"]["source_commit"] == git(author, "rev-parse", "HEAD")
            # Sharing the runtime does not permit selecting another corpus/database.
            assert (
                await clients[0][0].call_tool("table_details", {"table_name": details[1]["name"]})
            ).isError

    asyncio.run(exercise())


def test_service_restarts_after_crash_and_rejects_incompatible_clients(runtime, repos, tmp_path):
    _, worker, _ = repos
    db_path = tmp_path / "db"
    source = MainSource.discover(worker)
    client = SearchClient(source, db_path, source.table, 300)

    async def exercise():
        first = await client.call("query", query="numeric solve", query_type="vector", top_k=3)
        old_pid = probe(db_path)["pid"]
        with pytest.raises(RuntimeError, match="different code/dependencies"):
            await rpc(db_path, "mount", implementation="incompatible", mount=client.mount)
        with pytest.raises(RuntimeError, match="different upstream"):
            wrong = SearchClient(
                MainSource.discover(worker, corpus="articles"), db_path, source.table, 300
            )
            await wrong.call("mount")
        os.kill(old_pid, signal.SIGKILL)
        second = await client.call("query", query="numeric solve", query_type="vector", top_k=3)
        assert probe(db_path)["pid"] != old_pid
        assert first["results"] == second["results"]
        assert first["freshness"]["snapshot_id"] == second["freshness"]["snapshot_id"]
        await client.close()

    asyncio.run(exercise())


def test_cli_status_and_plan_do_not_start_background_publication(runtime, repos, tmp_path):
    _, worker, _ = repos
    db_path = tmp_path / "db"
    source = MainSource.discover(worker)

    def command(name):
        return json.loads(
            subprocess.check_output(
                [
                    sys.executable,
                    "-m",
                    "tools.agent_search",
                    name,
                    "--root",
                    str(worker),
                    "--db",
                    str(db_path),
                ]
            )
        )

    assert command("status")["integrity"] == "absent"
    assert probe(db_path)["corpora"] == {}
    assert probe(db_path)["loaded_models"] == 0
    assert command("plan")["mode"] == "plan"
    assert runtime.read_manifest(db_path, source.table) is None
    assert not (db_path / f"{source.table}.refresh.json").exists()
    assert probe(db_path)["corpora"] == {}


def test_unmount_and_expired_lease_allow_idle_shutdown(runtime, repos, tmp_path):
    from tools.agent_search.service import SearchService

    _, worker, _ = repos
    db_path = tmp_path / "db"
    service = SearchService(db_path, idle_timeout=0.1)
    source = MainSource.discover(worker)
    client = SearchClient(source, db_path, source.table, 300)

    async def exercise():
        reaper = asyncio.create_task(service.reap())
        try:
            # Resolve/lease the mount normally, then let its real refresh finish.
            _, corpus = await service.register(client.mount)
            await corpus.first_attempt.wait()
            assert not service.stop.is_set()
            corpus.mounts[client.mount["client_id"]]["expires"] = time.monotonic() - 1
            await asyncio.wait_for(service.stop.wait(), timeout=3)
            assert service.corpora == {}
            assert corpus.task.done()
        finally:
            reaper.cancel()
            await asyncio.gather(reaper, return_exceptions=True)

    asyncio.run(exercise())


def test_refresh_fails_over_from_a_removed_checkout(runtime, repos, tmp_path):
    from tools.agent_search.service import SearchService

    author, worker, _ = repos
    db_path = tmp_path / "db"
    service = SearchService(db_path)
    first = MainSource.discover(author)
    second = MainSource.discover(worker)
    client_a = SearchClient(first, db_path, first.table, 0.2)
    client_b = SearchClient(second, db_path, second.table, 0.2)

    async def exercise():
        _, corpus = await service.register(client_a.mount)
        await service.register(client_b.mount)
        try:
            await corpus.first_attempt.wait()
            (author / "docs/guide.md").write_text("# New main flux\nNew outgoing channel.\n")
            newest = commit(author)
            git(author, "push", "-q", "origin", "main")
            author.rename(tmp_path / "removed-checkout")
            corpus.wake.set()
            async with asyncio.timeout(10):
                while True:
                    report = await asyncio.to_thread(runtime.status, second, db_path, second.table)
                    if (
                        report["source_commit"] == newest
                        and report["refresh"]["state"] == "complete"
                    ):
                        break
                    await asyncio.sleep(0.05)
            assert corpus.error is None
        finally:
            corpus.mounts.clear()
            corpus.wake.set()
            await corpus.task

    asyncio.run(exercise())
