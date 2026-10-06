"""Exercise main-only provenance and publication with real Git, embeddings, and MCP."""

import asyncio
import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from tools.agent_search.chunks import markdown_units, python_units
from tools.agent_search.git_source import MainSource, repository_id

PROFILE = json.loads((Path(__file__).parents[1] / ".opencode/search/profiles.json").read_text())


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args]).decode().strip()


def commit(root, message="Change"):
    git(root, "add", ".")
    return git(
        root,
        "-c",
        "user.name=Index Test",
        "-c",
        "user.email=index@example.invalid",
        "commit",
        "-qm",
        message,
    ) or git(root, "rev-parse", "HEAD")


@pytest.fixture
def runtime():
    pytest.importorskip("lancedb")
    pytest.importorskip("sentence_transformers")
    pytest.importorskip("tree_sitter_rust")
    from tools.agent_search import __main__ as indexer

    return indexer


@pytest.fixture
def repos(tmp_path):
    upstream = tmp_path / "upstream.git"
    subprocess.run(
        ["git", "init", "--bare", "-q", "--initial-branch=main", str(upstream)], check=True
    )
    author = tmp_path / "author"
    subprocess.run(["git", "clone", "-q", str(upstream), str(author)], check=True)
    sources = {
        "libs/solver.py": (
            '"""Cached sparse symbolic analysis."""\n\nclass Solver:\n'
            "    def solve(self, value):\n"
            '        """Refactor the numeric matrix while reusing symbolic ordering."""\n'
            "        return value * 2\n"
        ),
        "docs/guide.md": "# Sparse solves\nReuse symbolic ordering across shifted matrices.\n",
        "native/kernel.rs": (
            "use std::f64;\n\n/// Euclidean norm.\nfn norm(x: f64) -> f64 { x.abs() }\n"
        ),
        "docs/superpowers/plan.md": "# Temporary plan\nExcluded.\n",
        ".opencode/search/profiles.json": json.dumps(PROFILE),
    }
    for name, text in sources.items():
        path = author / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    commit(author, "Fixture")
    git(author, "push", "-q", "origin", "main")
    worker = tmp_path / "worker"
    subprocess.run(["git", "clone", "-q", str(upstream), str(worker)], check=True)
    git(worker, "switch", "-qc", "feature")
    return author, worker, upstream


def publish(runtime, root, db_path):
    source = MainSource.discover(root)
    return runtime.ensure_main(source, db_path, source.table, force=True)


def rows(runtime, root, db_path):
    source = MainSource.discover(root)
    manifest = runtime.read_manifest(db_path, source.table)
    return runtime.open_published(db_path, manifest).to_arrow().to_pylist()


def test_python_method_anchor_survives_unrelated_line_shifts():
    before = "class Solver:\n    def solve(self, value):\n        return value * 2\n"
    after = "def helper():\n    return 1\n\n" + before
    old = next(u for u in python_units(before) if u.anchor == "Solver.solve")
    new = next(u for u in python_units(after) if u.anchor == "Solver.solve")
    assert old.content == new.content
    assert new.start_line == old.start_line + 3


def test_markdown_fenced_headings_stay_inside_their_section():
    source = "# Guide\n## Example\n```python\n# this is code\n```\n## Next\nDone.\n"
    units = markdown_units(source)
    assert [u.anchor for u in units] == ["Guide", "Guide > Example", "Guide > Next"]
    assert "# this is code" in units[1].content


def test_repository_identity_normalizes_transports_and_removes_credentials(tmp_path):
    ssh = repository_id("git@github.com:VanaMartin/qscat.git", tmp_path)
    https = repository_id("https://token:secret@github.com/vanamartin/qscat.git", tmp_path)
    assert ssh == https == "github.com/vanamartin/qscat"
    assert repository_id("ssh://git@github.com/VanaMartin/qscat", tmp_path) == ssh
    assert repository_id("ssh://git@example.com:443/repo", tmp_path) != repository_id(
        "https://example.com/repo", tmp_path
    )


def test_branch_commits_staged_dirty_new_and_policy_edits_never_enter_main(
    runtime, repos, tmp_path
):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    (worker / "libs/branch.py").write_text("def branch_only():\n    return 99\n")
    commit(worker, "Branch only")
    (worker / "libs/solver.py").write_text("def broken(\n")
    git(worker, "add", "libs/solver.py")
    (worker / "docs/guide.md").write_text("# Uncommitted guidance\nTransient.\n")
    (worker / "libs/untracked.py").write_text("def untracked():\n    pass\n")
    (worker / ".opencode/search/profiles.json").write_text("not JSON")
    before = git(worker, "status", "--porcelain=v1")
    result = publish(runtime, worker, db_path)
    assert result["source_commit"] == git(author, "rev-parse", "HEAD")
    current = rows(runtime, worker, db_path)
    assert "libs/branch.py" not in {row["path"] for row in current}
    assert "libs/untracked.py" not in {row["path"] for row in current}
    assert any("value * 2" in row["doc"] for row in current)
    assert git(worker, "status", "--porcelain=v1") == before
    assert git(worker, "branch", "--show-current") == "feature"
    source = MainSource.discover(worker)
    manifest = runtime.read_manifest(db_path, source.table)
    assert manifest["profile_sha256"] == runtime.fingerprint(PROFILE)
    delta = git(worker, "diff", "--name-only", manifest["source_commit"])
    assert "libs/branch.py" in delta and "libs/solver.py" in delta
    assert "libs/untracked.py" in git(worker, "ls-files", "--others", "--exclude-standard")
    git(author, "fetch", "-q", str(worker), "feature")
    git(author, "merge", "--ff-only", "FETCH_HEAD")
    git(author, "push", "-q", "origin", "main")
    accepted = publish(runtime, worker, db_path)
    assert accepted["source_commit"] == git(worker, "rev-parse", "HEAD")
    assert "libs/branch.py" in {row["path"] for row in rows(runtime, worker, db_path)}
    assert git(worker, "status", "--porcelain=v1") == before


def test_main_reconciliation_reuses_vectors_and_retires_deleted_and_renamed_sources(
    runtime, repos, tmp_path
):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    first = publish(runtime, worker, db_path)
    assert first["files"] == 4 and first["new_embeddings"] > 0
    before = rows(runtime, worker, db_path)
    original = next(row for row in before if row["anchor"] == "Solver.solve")
    source = MainSource.discover(worker)
    manifest_before = runtime.read_manifest(db_path, source.table)
    unchanged = publish(runtime, worker, db_path)
    assert (unchanged["parsed_files"], unchanged["new_embeddings"], unchanged["updated_rows"]) == (
        0,
        0,
        0,
    )
    assert runtime.read_manifest(db_path, source.table) == manifest_before

    path = author / "libs/solver.py"
    path.write_text("def helper():\n    return 1\n\n" + path.read_text())
    commit(author)
    git(author, "push", "-q", "origin", "main")
    shifted_result = publish(runtime, worker, db_path)
    assert shifted_result["parsed_files"] == 1
    shifted = next(row for row in rows(runtime, worker, db_path) if row["anchor"] == "Solver.solve")
    for key in ("id", "text", "payload_sha256", "embedding_key", "vector"):
        assert shifted[key] == original[key]
    assert shifted["start_line"] == original["start_line"] + 3
    assert shifted["source_blob"] != original["source_blob"]
    path.write_text(path.read_text().replace("value * 2", "value * 3"))
    commit(author)
    git(author, "push", "-q", "origin", "main")
    edited = publish(runtime, worker, db_path)
    assert edited["new_embeddings"] == edited["deleted_rows"] == 1
    (author / "docs/guide.md").unlink()
    path.rename(author / "libs/renamed.py")
    commit(author)
    git(author, "push", "-q", "origin", "main")
    publish(runtime, worker, db_path)
    current = rows(runtime, worker, db_path)
    assert "libs/solver.py" not in {row["path"] for row in current}
    assert "docs/guide.md" not in {row["path"] for row in current}
    assert "libs/renamed.py" in {row["path"] for row in current}
    assert len({row["id"] for row in current}) == len(current)
    for mode in ("vector", "fts", "hybrid"):
        found = runtime.query(
            source, db_path, source.table, "symbolic ordering numeric refactor", mode, 3
        )
        assert any(row["path"] == "libs/renamed.py" for row in found["results"])
        assert all("vector" not in row for row in found["results"])
        assert found["freshness"]["source_commit"] == git(author, "rev-parse", "HEAD")


def test_commit_without_selected_changes_advances_provenance_without_rebuilding(
    runtime, repos, tmp_path
):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    before = runtime.read_manifest(db_path, source.table)
    (author / "docs/superpowers/plan.md").write_text("# Excluded update\n")
    commit(author)
    git(author, "push", "-q", "origin", "main")
    result = publish(runtime, worker, db_path)
    after = runtime.read_manifest(db_path, source.table)
    assert result["new_embeddings"] == result["parsed_files"] == result["updated_rows"] == 0
    assert after["published_table"] == before["published_table"]
    assert after["source_commit"] != before["source_commit"]


def test_main_parse_failure_and_offline_clone_preserve_last_published_snapshot(
    runtime, repos, tmp_path
):
    author, worker, upstream = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    before = runtime.read_manifest(db_path, source.table)
    old_rows = rows(runtime, worker, db_path)
    (author / "libs/solver.py").write_text("def incomplete(\n")
    commit(author)
    git(author, "push", "-q", "origin", "main")
    with pytest.raises(ValueError, match=r"Cannot index libs/solver\.py"):
        publish(runtime, worker, db_path)
    assert runtime.read_manifest(db_path, source.table) == before
    assert rows(runtime, worker, db_path) == old_rows
    report = runtime.status(source, db_path, source.table)
    assert report["integrity"] == "complete" and report["upstream"]["state"] == "behind"
    assert report["refresh"]["state"] == "failed"
    assert runtime.query(source, db_path, source.table, "numeric solve", "fts", 2)["results"]

    # The clone still has old local refs; an outage must not republish those refs.
    upstream.rename(tmp_path / "offline.git")
    with pytest.raises(RuntimeError, match="Cannot fetch"):
        publish(runtime, worker, db_path)
    assert runtime.read_manifest(db_path, source.table) == before
    report = runtime.status(source, db_path, source.table)
    assert report["upstream"]["state"] == "unknown"
    assert report["refresh"]["state"] == "failed"


def test_clones_share_corpus_and_changed_upstream_is_rejected(runtime, repos, tmp_path):
    author, worker, upstream = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    clone = tmp_path / "another-clone"
    subprocess.run(["git", "clone", "-q", str(upstream), str(clone)], check=True)
    other = MainSource.discover(clone)
    assert source.table == other.table
    before = runtime.read_manifest(db_path, source.table)
    assert publish(runtime, clone, db_path)["new_embeddings"] == 0
    assert runtime.read_manifest(db_path, source.table) == before
    git(clone, "remote", "set-url", "origin", str(author))
    with pytest.raises(ValueError, match="different upstream"):
        runtime.ensure_main(MainSource.discover(clone), db_path, source.table, force=True)


def test_mounted_upstream_identity_change_and_absolute_policy_are_rejected(repos):
    author, worker, _ = repos
    source = MainSource.discover(worker)
    git(worker, "remote", "set-url", "origin", str(author))
    with pytest.raises(ValueError, match="upstream changed"):
        source.fetch()
    with pytest.raises(ValueError, match="repository-relative"):
        MainSource.discover(worker, profile=str(author / ".opencode/search/profiles.json"))


def test_publication_keeps_old_readers_pinned_and_failed_candidate_is_invisible(
    runtime, repos, tmp_path, monkeypatch
):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    before = runtime.read_manifest(db_path, source.table)
    (author / "docs/guide.md").write_text(
        "# New quantum algorithm\nStable accepted implementation.\n"
    )
    commit(author)
    git(author, "push", "-q", "origin", "main")
    write_json = runtime.write_json
    observations = []

    def fail_publication(path, value):
        if path.name == f"{source.table}.manifest.json":
            found = runtime.query(source, db_path, source.table, "symbolic ordering", "fts", 3)
            observations.append(found["freshness"]["source_commit"])
            raise OSError("Interrupted publication")
        write_json(path, value)

    with monkeypatch.context() as context:
        context.setattr(runtime, "write_json", fail_publication)
        with pytest.raises(OSError, match="Interrupted publication"):
            publish(runtime, worker, db_path)
    assert observations == [before["source_commit"]]
    assert runtime.read_manifest(db_path, source.table) == before
    publish(runtime, worker, db_path)
    newest = runtime.read_manifest(db_path, source.table)
    assert newest["published_table"] != before["published_table"]
    pinned = runtime.query(
        source, db_path, source.table, "symbolic ordering", "fts", 3, manifest=before
    )
    assert pinned["freshness"]["source_commit"] == before["source_commit"]
    assert any("Reuse symbolic" in row["doc"] for row in pinned["results"])


def test_external_mutation_is_repaired_in_a_new_generation(runtime, repos, tmp_path):
    import lancedb

    _, worker, _ = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    before = runtime.read_manifest(db_path, source.table)
    table = lancedb.connect(str(db_path)).open_table(before["published_table"])
    table.delete("path = 'docs/guide.md'")
    assert runtime.status(source, db_path, source.table)["integrity"] == "incomplete"
    with pytest.raises(RuntimeError, match="incomplete"):
        runtime.query(source, db_path, source.table, "numeric solve", "fts", 2)
    publish(runtime, worker, db_path)
    assert runtime.status(source, db_path, source.table)["state"] == "current"
    assert (
        runtime.read_manifest(db_path, source.table)["published_table"] != before["published_table"]
    )


def test_committed_policy_changes_rechunk_and_branch_policy_is_ignored(runtime, repos, tmp_path):
    author, worker, _ = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    changed = copy.deepcopy(PROFILE)
    changed["embedding"]["chunk_target_tokens"] = 120
    (author / ".opencode/search/profiles.json").write_text(json.dumps(changed))
    commit(author)
    git(author, "push", "-q", "origin", "main")
    result = publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    assert result["parsed_files"] == 4
    assert runtime.read_manifest(db_path, source.table)["profile_sha256"] == runtime.fingerprint(
        changed
    )


def test_real_tokenizer_budget_duplicate_declarations_and_rust_context(runtime, repos):
    from tools.agent_search.chunks import rust_units, split_units

    rust = (
        "struct Solver;\nimpl Solver {\n    /// Numeric solve.\n"
        "    #[inline]\n    fn solve(&self) {}\n}\n"
    )
    method = next(unit for unit in rust_units(rust) if unit.anchor == "Solver::solve")
    assert "#[inline]" in method.content and "/// Numeric solve." in method.content
    model, spec = runtime.embedding(PROFILE)
    tokenizer = model.embedding_model.tokenizer
    source = ('def same():\n    return "' + "long scientific content " * 200 + '"\n') * 2
    _, chunks = split_units("libs/example.py", source, tokenizer, PROFILE["embedding"])
    assert {chunk.declaration for chunk in chunks} == {0, 1}
    assert all(chunk.tokens == len(tokenizer(chunk.text)["input_ids"]) <= 256 for chunk in chunks)
    spaced = "A scientific paragraph about nuclear scattering.\n\n" * 100
    _, paragraphs = split_units("docs/example.md", spaced, tokenizer, PROFILE["embedding"])
    lines = spaced.splitlines(keepends=True)
    assert any(chunk.content.endswith("\n\n") for chunk in paragraphs)
    assert all(
        chunk.content in "".join(lines[chunk.start_line - 1 : chunk.end_line])
        for chunk in paragraphs
    )
    snapshot = MainSource.discover(repos[1]).fetch()
    records, _, _ = runtime.make_records(snapshot, tokenizer, runtime.fingerprint(spec))
    assert all(
        "source_sha256=" not in row["text"]
        and "lines=" not in row["text"]
        and "source_commit=" not in row["text"]
        for row in records
    )


def test_legacy_migration_reuses_payloads_but_regenerates_main_provenance(runtime, repos, tmp_path):
    import lancedb

    _, worker, _ = repos
    db_path = tmp_path / "db"
    publish(runtime, worker, db_path)
    source = MainSource.discover(worker)
    manifest = runtime.read_manifest(db_path, source.table)
    current = rows(runtime, worker, db_path)
    db = lancedb.connect(str(db_path))
    db.create_table("legacy", data=current)
    legacy = {
        "format": "qmodeling-search-v1",
        "root": str(worker),
        "table": "legacy",
        "state": "ready",
        "table_version": db.open_table("legacy").version,
        "chunks": len(current),
        "embedding": manifest["embedding"],
    }
    runtime.write_json(db_path / "legacy.manifest.json", legacy)
    (db_path / f"{source.table}.manifest.json").unlink()
    (worker / "libs/solver.py").write_text("def incomplete(\n")
    result = publish(runtime, worker, db_path)
    assert result["new_embeddings"] == 0
    assert any("value * 2" in row["doc"] for row in rows(runtime, worker, db_path))
    assert runtime.read_manifest(db_path, source.table)["format"] == runtime.VERSION


def parameters(root, db_path, interval=300):
    from mcp import StdioServerParameters

    return StdioServerParameters(
        command=sys.executable,
        args=[
            "-m",
            "tools.agent_search.mcp",
            "--root",
            str(root),
            "--db",
            str(db_path),
            "--check-interval",
            str(interval),
        ],
    )


def decoded(result):
    assert not result.isError, result.content
    return json.loads(result.content[0].text)


def test_concurrent_mcp_mounts_share_one_initial_generation(runtime, repos, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client

    author, worker, _ = repos
    db_path = tmp_path / "db"

    async def mount(root):
        async with stdio_client(parameters(root, db_path)) as (read, write):
            async with ClientSession(read, write) as client:
                await asyncio.wait_for(client.initialize(), timeout=10)
                catalog = await client.list_tools()
                assert {tool.name for tool in catalog.tools} == {"query_table", "table_details"}
                found = decoded(
                    await client.call_tool(
                        "query_table",
                        {
                            "query": "symbolic ordering numeric refactor",
                            "top_k": 3,
                            "query_type": "hybrid",
                        },
                    )
                )
                assert found["freshness"]["source_commit"] == git(author, "rev-parse", "HEAD")
                return decoded(await client.call_tool("table_details", {}))

    async def exercise():
        return await asyncio.wait_for(asyncio.gather(mount(worker), mount(author)), timeout=90)

    first, second = asyncio.run(exercise())
    assert first["name"] == second["name"] == MainSource.discover(worker).table
    assert first["freshness"]["published_table"] == second["freshness"]["published_table"]
    import lancedb

    generations = [
        name
        for name in lancedb.connect(str(db_path)).table_names()
        if name.startswith(first["name"] + "_g_")
    ]
    assert len(generations) == 1


def test_mcp_refreshes_live_main_and_serves_last_good_snapshot_after_fetch_failure(
    runtime, repos, tmp_path
):
    pytest.importorskip("mcp")
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client

    author, worker, upstream = repos
    db_path = tmp_path / "db"

    async def wait_for(client, predicate):
        for _ in range(200):
            report = decoded(await client.call_tool("table_details", {}))
            if predicate(report["freshness"]):
                return report
            await asyncio.sleep(0.1)
        pytest.fail("Periodic main refresh did not reach expected state")

    async def exercise():
        async with stdio_client(parameters(worker, db_path, interval=0.5)) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                initial = decoded(await client.call_tool("query_table", {"query": "numeric solve"}))
                assert initial["freshness"]["integrity"] == "complete"
                (author / "docs/guide.md").write_text(
                    "# Published quantum flux\nOutgoing flux is stable.\n"
                )
                newest = commit(author)
                git(author, "push", "-q", "origin", "main")
                report = await wait_for(
                    client,
                    lambda info: (
                        info["source_commit"] == newest and info["refresh"]["state"] == "complete"
                    ),
                )
                assert report["freshness"]["integrity"] == "complete"
                found = decoded(
                    await client.call_tool(
                        "query_table",
                        {
                            "query": "Published quantum flux",
                            "query_type": "fts",
                        },
                    )
                )
                assert any("Outgoing flux" in hit["doc"] for hit in found["results"])
                upstream.rename(tmp_path / "offline.git")
                await wait_for(client, lambda info: info["refresh"]["state"] == "failed")
                offline = decoded(await client.call_tool("query_table", {"query": "numeric solve"}))
                assert offline["freshness"]["source_commit"] == newest
                assert offline["freshness"]["upstream"]["state"] == "unknown"

    asyncio.run(exercise())


def test_first_mount_fetch_failure_is_inspectable(runtime, repos, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client

    _, worker, upstream = repos
    upstream.rename(tmp_path / "offline.git")

    async def exercise():
        async with stdio_client(parameters(worker, tmp_path / "db")) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                result = await client.call_tool("query_table", {"query": "numeric solve"})
                assert result.isError
                report = decoded(await client.call_tool("table_details", {}))
                assert report["freshness"]["integrity"] == "absent"
                assert report["freshness"]["refresh"]["state"] == "failed"
                assert "Cannot fetch" in report["freshness"]["refresh"]["error"]

    asyncio.run(exercise())


def test_first_mount_writer_failure_does_not_leave_query_waiting_forever(runtime, repos, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client

    db_path = tmp_path / "not-a-directory"
    db_path.write_text("Existing unrelated file")

    async def exercise():
        async with stdio_client(parameters(repos[1], db_path)) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                result = await asyncio.wait_for(
                    client.call_tool("query_table", {"query": "numeric solve"}), timeout=10
                )
                assert result.isError
                assert "Initial main indexing failed" in result.content[0].text
                report = decoded(await client.call_tool("table_details", {}))
                assert report["maintainer_error"]
                assert report["freshness"]["integrity"] == "absent"

    asyncio.run(exercise())
