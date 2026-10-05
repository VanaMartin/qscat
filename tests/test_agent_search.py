"""Exercise derived-index freshness against real parsers, embeddings, and LanceDB."""

import asyncio
import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from tools.agent_search.chunks import markdown_units, python_units

PROFILE = json.loads((Path(__file__).parents[1] / ".opencode/search/profiles.json").read_text())


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


@pytest.fixture
def runtime():
    pytest.importorskip("lancedb")
    pytest.importorskip("sentence_transformers")
    pytest.importorskip("tree_sitter_rust")
    from tools.agent_search import __main__ as indexer

    return indexer


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "checkout"
    root.mkdir()
    subprocess.run(["git", "init", "-q", str(root)], check=True)
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
    }
    for name, text in sources.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "-c",
            "user.name=Index Test",
            "-c",
            "user.email=index@example.invalid",
            "commit",
            "-qm",
            "Fixture",
        ],
        check=True,
    )
    return root


def rows(db_path):
    import lancedb

    return lancedb.connect(str(db_path)).open_table("code").to_arrow().to_pylist()


def test_rust_item_anchors_include_attributes_and_impl_context(runtime):
    from tools.agent_search.chunks import rust_units

    source = (
        "struct Solver;\nimpl Solver {\n    /// Numeric solve.\n"
        "    #[inline]\n    fn solve(&self) {}\n}\n"
    )
    units = rust_units(source)
    method = next(u for u in units if u.anchor == "Solver::solve")
    assert "#[inline]" in method.content
    assert "/// Numeric solve." in method.content
    assert method.start_line == 3


def test_population_reuse_deletion_rename_and_search(runtime, repo, tmp_path):
    db_path = tmp_path / "db"
    first = runtime.reconcile(repo, db_path, "code", PROFILE, [])
    assert first["files"] == 3
    assert first["new_embeddings"] > 0
    assert first["max_tokens"] <= 256
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "current"
    original = next(r for r in rows(db_path) if r["anchor"] == "Solver.solve")
    unchanged = runtime.reconcile(repo, db_path, "code", PROFILE, [])
    assert (unchanged["new_embeddings"], unchanged["updated_rows"], unchanged["deleted_rows"]) == (
        0,
        0,
        0,
    )

    path = repo / "libs/solver.py"
    path.write_text("def helper():\n    return 1\n\n" + path.read_text())
    assert runtime.status(repo, db_path, "code", PROFILE)["changed_paths"] == ["libs/solver.py"]
    runtime.reconcile(repo, db_path, "code", PROFILE, [])
    shifted = next(r for r in rows(db_path) if r["anchor"] == "Solver.solve")
    for key in ("id", "text", "payload_sha256", "embedding_key", "vector"):
        assert shifted[key] == original[key]
    assert shifted["start_line"] == original["start_line"] + 3
    assert shifted["source_sha256"] != original["source_sha256"]

    path.write_text(path.read_text().replace("value * 2", "value * 3"))
    edited = runtime.reconcile(repo, db_path, "code", PROFILE, [])
    assert edited["new_embeddings"] == 1
    assert edited["deleted_rows"] == 1
    assert original["id"] not in {r["id"] for r in rows(db_path)}
    (repo / "docs/guide.md").unlink()
    path.rename(repo / "libs/renamed.py")
    runtime.reconcile(repo, db_path, "code", PROFILE, ["libs/renamed.py"])
    current = rows(db_path)
    assert {r["path"] for r in current} == {"libs/renamed.py", "native/kernel.rs"}
    assert len({r["id"] for r in current}) == len(current)
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "current"
    for mode in ("vector", "fts", "hybrid"):
        hits = runtime.query(db_path, "code", "symbolic ordering numeric refactor", mode, 3)
        assert any(h["path"] == "libs/renamed.py" for h in hits)
        assert all("vector" not in h for h in hits)


def test_parse_failure_retains_committed_rows_and_reports_stale(runtime, repo, tmp_path):
    db_path = tmp_path / "db"
    runtime.reconcile(repo, db_path, "code", PROFILE, [])
    before = rows(db_path)
    (repo / "libs/solver.py").write_text("def incomplete(\n")
    with pytest.raises(ValueError, match=r"Cannot index libs/solver\.py"):
        runtime.reconcile(repo, db_path, "code", PROFILE, [])
    assert rows(db_path) == before
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "stale"


def test_worktree_ownership_and_task_owned_untracked_files(runtime, repo, tmp_path):
    db_path = tmp_path / "db"
    (repo / "libs/new.py").write_text("def owned_new():\n    return 1\n")
    (repo / "libs/output.py").write_text("def untracked_output():\n    return 99\n")
    result = runtime.reconcile(repo, db_path, "code", PROFILE, ["libs/new.py"])
    assert result["files"] == 4
    assert "libs/output.py" not in {r["path"] for r in rows(db_path)}
    other = tmp_path / "other-checkout"
    other.mkdir()
    with pytest.raises(ValueError, match="different worktree"):
        runtime.reconcile(other, db_path, "code", PROFILE, [])
    with pytest.raises(ValueError, match="different worktree"):
        runtime.status(other, db_path, "code", PROFILE)
    policy = copy.deepcopy(PROFILE)
    policy["embedding"]["chunk_target_tokens"] = 120
    assert runtime.status(repo, db_path, "code", policy)["profile_changed"]


def test_token_budget_and_duplicate_declarations(runtime, repo):
    from tools.agent_search.chunks import split_units

    model, spec = runtime.embedding(PROFILE)
    tokenizer = model.embedding_model.tokenizer
    source = ('def same():\n    return "' + "long scientific content " * 200 + '"\n') * 2
    _, chunks = split_units("libs/example.py", source, tokenizer, PROFILE["embedding"])
    assert {c.declaration for c in chunks} == {0, 1}
    assert all(c.tokens == len(tokenizer(c.text)["input_ids"]) <= 256 for c in chunks)
    assert all(c.text.endswith(c.content) for c in chunks)
    # Splitting at a blank line must include that physical line in the locator.
    spaced = "A scientific paragraph about nuclear scattering.\n\n" * 100
    _, paragraphs = split_units("docs/example.md", spaced, tokenizer, PROFILE["embedding"])
    lines = spaced.splitlines(keepends=True)
    assert any(c.content.endswith("\n\n") for c in paragraphs)
    assert all(c.content in "".join(lines[c.start_line - 1 : c.end_line]) for c in paragraphs)
    records, _ = runtime.make_records(
        repo, PROFILE, ["libs/solver.py"], tokenizer, runtime.fingerprint(spec)
    )
    assert all("source_sha256=" not in r["text"] and "lines=" not in r["text"] for r in records)


def test_external_table_mutation_is_incomplete(runtime, repo, tmp_path):
    import lancedb

    db_path = tmp_path / "db"
    runtime.reconcile(repo, db_path, "code", PROFILE, [])
    table = lancedb.connect(str(db_path)).open_table("code")
    table.delete("path = 'docs/guide.md'")
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "incomplete"
    runtime.reconcile(repo, db_path, "code", PROFILE, [])
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "current"


def test_mount_bootstrap_adopts_empty_schema_and_skips_populated_index(runtime, repo, tmp_path):
    import lancedb
    import pyarrow as pa

    db_path = tmp_path / "db"
    db = lancedb.connect(str(db_path))
    db.create_table("code", schema=pa.schema([("doc", pa.string())]))
    runtime.reconcile(repo, db_path, "code", PROFILE, [], only_if_empty=True)
    assert "anchor" in db.open_table("code").schema.names
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "current"
    before = rows(db_path)
    manifest = runtime.read_manifest(db_path, "code")
    (repo / "libs/solver.py").write_text("def partial(\n")
    skipped = runtime.reconcile(repo, db_path, "code", PROFILE, [], only_if_empty=True)
    assert skipped["state"] == "populated"
    assert skipped["new_embeddings"] == 0
    assert rows(db_path) == before
    assert runtime.read_manifest(db_path, "code") == manifest
    assert runtime.status(repo, db_path, "code", PROFILE)["state"] == "stale"


def test_new_clone_selects_and_populates_an_independent_table(runtime, repo, tmp_path):
    db_path = tmp_path / "db"
    runtime.reconcile(repo, db_path, "qscat_knowledge", PROFILE, [])
    owner_before = runtime.read_manifest(db_path, "qscat_knowledge")
    assert runtime.select_table(repo, db_path) == "qscat_knowledge"
    clone = tmp_path / "new-clone"
    subprocess.run(["git", "clone", "-q", str(repo), str(clone)], check=True)
    table = runtime.select_table(clone, db_path)
    assert table.startswith("qscat_knowledge_")
    assert table != runtime.select_table(repo, db_path)
    runtime.reconcile(clone, db_path, table, PROFILE, [], only_if_empty=True)
    assert runtime.status(clone, db_path, table, PROFILE)["state"] == "current"
    assert runtime.read_manifest(db_path, "qscat_knowledge") == owner_before


def test_mcp_mount_initializes_and_concurrent_mounts_do_not_duplicate(runtime, repo, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    db_path = tmp_path / "db"
    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-m",
            "tools.agent_search.mcp",
            "--root",
            str(repo),
            "--db",
            str(db_path),
            "--profile",
            str(Path(__file__).parents[1] / ".opencode/search/profiles.json"),
        ],
    )

    def decoded(result):
        assert not result.isError, result.content
        return json.loads(result.content[0].text)

    async def mount():
        async with stdio_client(parameters) as (read, write):
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
                assert found["freshness"]["state"] == "current"
                assert any(hit["path"] == "libs/solver.py" for hit in found["results"])
                details = decoded(await client.call_tool("table_details", {}))
                assert details["bootstrap"]["state"] == "complete"
                assert details["num_rows"] == found["freshness"]["chunks"]
                return details

    async def exercise():
        return await asyncio.wait_for(asyncio.gather(mount(), mount()), timeout=90)

    first, second = asyncio.run(exercise())
    assert first["name"] == second["name"] == runtime.select_table(repo, db_path)
    states = {first["bootstrap"]["result"].get("state"), second["bootstrap"]["result"].get("state")}
    assert states == {None, "populated"}
    import lancedb

    records = lancedb.connect(str(db_path)).open_table(first["name"]).to_arrow().to_pylist()
    assert len({r["id"] for r in records}) == len(records) == first["num_rows"]


def test_mcp_bootstrap_failure_is_inspectable(runtime, repo, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    (repo / "libs/solver.py").write_text("def incomplete(\n")
    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-m",
            "tools.agent_search.mcp",
            "--root",
            str(repo),
            "--db",
            str(tmp_path / "db"),
            "--profile",
            str(Path(__file__).parents[1] / ".opencode/search/profiles.json"),
        ],
    )

    async def exercise():
        async with stdio_client(parameters) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                result = await client.call_tool("query_table", {"query": "numeric solve"})
                assert result.isError
                assert "Cannot index libs/solver.py" in result.content[0].text
                details = await client.call_tool("table_details", {})
                assert not details.isError
                report = json.loads(details.content[0].text)
                assert report["bootstrap"]["state"] == "failed"
                assert report["num_rows"] == 0

    asyncio.run(exercise())


def test_mcp_reads_live_policy_after_checkpoint_catchup(runtime, repo, tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    db_path = tmp_path / "db"
    profile_path = repo / "policy.json"
    profile_path.write_text(json.dumps(PROFILE))
    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-m",
            "tools.agent_search.mcp",
            "--root",
            str(repo),
            "--db",
            str(db_path),
            "--profile",
            str(profile_path),
        ],
    )

    async def exercise():
        async with stdio_client(parameters) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                result = await client.call_tool("query_table", {"query": "numeric solve"})
                assert not result.isError
                changed = copy.deepcopy(PROFILE)
                changed["embedding"]["chunk_target_tokens"] = 120
                profile_path.write_text(json.dumps(changed))
                details = await client.call_tool("table_details", {})
                assert json.loads(details.content[0].text)["freshness"]["profile_changed"]
                table = runtime.select_table(repo, db_path)
                await asyncio.to_thread(runtime.reconcile, repo, db_path, table, changed, [])
                caught_up = await client.call_tool("table_details", {})
                report = json.loads(caught_up.content[0].text)
                assert report["freshness"]["state"] == "current"
                assert not report["freshness"]["profile_changed"]

    asyncio.run(exercise())
