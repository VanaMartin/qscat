"""Keep repository search package-owned without a second Python search runtime."""

import fnmatch
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / ".opencode/search-package/runtime.json"
PROFILE = ROOT / ".opencode/search-package/profiles.json"


def selected(path, policy):
    """Check this repository's path policy, not the package's indexing engine."""
    return (
        any(fnmatch.fnmatchcase(path, pattern) for pattern in policy["include"])
        and not any(fnmatch.fnmatchcase(path, pattern) for pattern in policy["exclude"])
        and (
            Path(path).suffix in policy["extensions"]
            or any(
                fnmatch.fnmatchcase(path, pattern)
                for pattern in policy.get("extensionless_files", [])
            )
        )
    )


def test_search_registration_belongs_to_the_native_plugin():
    config = json.loads((ROOT / "opencode.json").read_text())
    servers = config.get("mcp", {}).get("servers", {})
    assert not {"lancedb", "lancedb-articles"}.intersection(servers)
    assert "tools.agent_search" not in json.dumps(config)


def test_native_settings_keep_upstream_source_and_explicit_local_policy():
    assert json.loads(SETTINGS.read_text()) == {
        "remote": "origin",
        "branch": "main",
        "profile": ".opencode/search-package/profiles.json",
        "policy_source": "local",
    }


def test_native_profile_keeps_independent_code_and_note_corpora():
    profile = json.loads(PROFILE.read_text())
    assert profile["format"] == "repo-search-profiles-v1"
    assert set(profile["corpora"]) == {"repository", "articles"}
    repository, articles = (profile["corpora"][name] for name in ("repository", "articles"))
    assert repository["enabled"] is articles["enabled"] is True
    assert repository["server"] == "lancedb"
    assert articles["server"] == "lancedb-articles"
    assert repository["parser"] == "source"
    assert articles["parser"] == "reference_note"
    assert articles["extensions"] == [".md"]
    embedding = profile["embedding"]
    assert embedding["model"] == "all-MiniLM-L6-v2"
    assert embedding["device"] == "cpu"
    assert embedding["dimensions"] == 384
    assert (
        0
        <= embedding["overlap_tokens"]
        < embedding["chunk_target_tokens"]
        <= embedding["chunk_max_tokens"]
        < embedding["total_input_max_tokens"]
        == 256
    )


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("libs/qscat/qscat/linalg/sparse_lu.py", {"repository"}),
        ("projects/ukrmol_co/model.py", {"repository"}),
        ("validation/n2/exact2d.py", {"repository"}),
        ("benchmarks/mumps_vs_superlu.py", {"repository"}),
        ("docs/physics/n2-resonance.md", {"repository"}),
        ("docker/base.Dockerfile", {"repository"}),
        (".claude/skills/knowledge-search/SKILL.md", {"repository"}),
        ("reference/literature/alt-houfek-2021-pra103-032829.md", {"articles"}),
        ("reference/literature/README.md", set()),
        ("reference/literature/article.pdf", set()),
        ("reference/eMoScat/solver.py", set()),
        ("docs/superpowers/plans/work.md", set()),
        ("validation/factory/results/output.py", set()),
        ("tests/fixtures/broken.py", set()),
    ],
)
def test_qmodeling_corpus_boundaries(path, expected):
    corpora = json.loads(PROFILE.read_text())["corpora"]
    assert {name for name, policy in corpora.items() if selected(path, policy)} == expected


def test_duplicate_search_engine_and_environment_are_retired():
    assert not list((ROOT / "tools/agent_search").glob("*.py"))
    for name in ("profiles.json", "pyproject.toml", "uv.lock"):
        assert not (ROOT / ".opencode/search" / name).exists()


@pytest.mark.parametrize(
    "module", ["tools.agent_search", "tools.agent_search.mcp", "tools.agent_search.service"]
)
def test_retired_entrypoints_cannot_start_a_legacy_service(module):
    result = subprocess.run(
        [sys.executable, "-B", "-m", module, "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=5,
    )
    assert result.returncode != 0
    assert "No module named" in result.stderr
    assert module in result.stderr


def test_active_guidance_uses_the_package_runtime_and_updater():
    for name in (
        "docs/agent-search.md",
        ".opencode/README.md",
        ".claude/skills/knowledge-search/SKILL.md",
    ):
        text = (ROOT / name).read_text()
        assert "tools.agent_search" not in text
        assert ".cache/qmodeling-search" not in text
        assert "repo-knowledge-search" in text
        assert "/search-status" in text
    guide = (ROOT / "docs/agent-search.md").read_text()
    assert "opencode plugin check 'github:VanaMartin/opencode-repo-search#latest'" in guide
    assert "opencode plugin update 'github:VanaMartin/opencode-repo-search#latest'" in guide
    assert "service stop --root" in guide


def test_native_specialists_can_load_the_package_skill():
    for path in (ROOT / ".opencode/agents").glob("*.md"):
        metadata = yaml.safe_load(path.read_text().split("---", 2)[1])
        assert {
            "action": "skill",
            "resource": "repo-knowledge-search",
            "effect": "allow",
        } in metadata["permissions"], path
    explore = json.loads((ROOT / "opencode.json").read_text())["agents"]["explore"]
    assert {"action": "skill", "resource": "repo-knowledge-search", "effect": "allow"} in explore[
        "permissions"
    ]
