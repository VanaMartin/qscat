"""The read side of the artifact store.

A published sweep is a function of (config, commit, image) that costs
minutes-to-hours to recompute, so it lives at `data.qscat.org` rather than in
git. What stays in git is a KB-sized `artifacts.json` pointer, and these tests
cover what that pointer has to guarantee: that a reader gets the bytes the
maintainer published, or a clear error -- never silently different bytes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from qscat_run.artifact_store import (
    ArtifactPointer,
    ArtifactStoreError,
    ChecksumMismatch,
    MissingPointer,
    fetch,
    load_pointer,
    pointer_path,
)


def _pointer(tmp_path: Path, payload: bytes = b"energy,sigma\n0.1,2.0\n") -> Path:
    d = tmp_path / "o2-ve"
    d.mkdir()
    (d / "artifacts.json").write_text(
        json.dumps(
            {
                "git_sha": "0" * 40,
                "url_prefix": "https://data.qscat.org/main/abc1234/o2-ve/",
                "artifacts": {
                    "cross_section.csv": {
                        "sha256": hashlib.sha256(payload).hexdigest(),
                        "bytes": len(payload),
                    }
                },
            }
        )
    )
    return d


def test_load_pointer_reads_the_published_locations(tmp_path: Path) -> None:
    p = load_pointer(_pointer(tmp_path))
    assert isinstance(p, ArtifactPointer)
    assert p.git_sha == "0" * 40
    assert p.url_for("cross_section.csv") == (
        "https://data.qscat.org/main/abc1234/o2-ve/cross_section.csv"
    )


def test_a_directory_without_a_pointer_says_so(tmp_path: Path) -> None:
    (tmp_path / "empty").mkdir()
    with pytest.raises(MissingPointer):
        load_pointer(tmp_path / "empty")


def test_fetch_writes_the_payload_and_verifies_it(tmp_path: Path) -> None:
    payload = b"energy,sigma\n0.1,2.0\n"
    d = _pointer(tmp_path, payload)
    written = fetch(d, opener=lambda url: payload)
    assert written == [d / "cross_section.csv"]
    assert (d / "cross_section.csv").read_bytes() == payload


def test_fetch_refuses_bytes_that_do_not_match_the_recorded_digest(tmp_path: Path) -> None:
    """The pointer is the only thing tying a downloaded file to the run that
    produced it. Corrupted or substituted bytes must fail loudly, and must not
    be left on disk looking like a successful fetch."""
    d = _pointer(tmp_path)
    with pytest.raises(ChecksumMismatch):
        fetch(d, opener=lambda url: b"not what was published")
    assert not (d / "cross_section.csv").exists()


def test_fetch_skips_a_file_already_present_and_correct(tmp_path: Path) -> None:
    payload = b"energy,sigma\n0.1,2.0\n"
    d = _pointer(tmp_path, payload)
    (d / "cross_section.csv").write_bytes(payload)

    def refuse(url: str) -> bytes:
        raise AssertionError(f"re-downloaded {url} despite a correct local copy")

    assert fetch(d, opener=refuse) == []


def test_fetch_replaces_a_file_whose_digest_is_wrong(tmp_path: Path) -> None:
    """A half-written or edited local copy is not a cache hit."""
    payload = b"energy,sigma\n0.1,2.0\n"
    d = _pointer(tmp_path, payload)
    (d / "cross_section.csv").write_bytes(b"stale")
    assert fetch(d, opener=lambda url: payload) == [d / "cross_section.csv"]
    assert (d / "cross_section.csv").read_bytes() == payload


def test_a_pointer_may_not_send_the_fetcher_at_the_local_filesystem(tmp_path: Path) -> None:
    """`urlopen` speaks file:// as readily as https. A committed pointer is
    reviewed, not validated, and reading a local path because a JSON field
    asked us to is never what a fetch means."""
    (tmp_path / "secret").write_bytes(b"local")
    d = tmp_path / "run"
    d.mkdir()
    (d / "artifacts.json").write_text(
        json.dumps(
            {
                "git_sha": "0" * 40,
                "url_prefix": f"file://{tmp_path}/",
                "artifacts": {
                    "secret": {"sha256": hashlib.sha256(b"local").hexdigest(), "bytes": 5}
                },
            }
        )
    )
    with pytest.raises(ArtifactStoreError, match="non-https"):
        fetch(d)


def test_the_request_identifies_itself_rather_than_defaulting_to_urllib() -> None:
    """Cloudflare answers `Python-urllib/3.x` with 403 in front of the bucket,
    so the default agent makes every fetch fail for every reader. Measured, not
    assumed: curl sent with that agent is refused and urllib sent with curl's
    is served."""
    from qscat_run.artifact_store import _user_agent

    agent = _user_agent()
    assert agent.startswith("qscat-run/")
    assert "urllib" not in agent.lower()
    # Says who is calling and where to complain, rather than posing as a browser.
    assert "github.com/VanaMartin/qscat" in agent


def test_pointer_path_is_the_documented_name(tmp_path: Path) -> None:
    assert pointer_path(tmp_path).name == "artifacts.json"


def test_every_committed_pointer_in_this_repo_is_wellformed() -> None:
    """Guards the pointers themselves: a typo in a URL prefix or a truncated
    digest is invisible until someone tries to fetch, which may be months
    later and on someone else's machine."""
    repo = Path(__file__).resolve().parents[3]
    pointers = sorted(repo.glob("validation/**/artifacts.json"))
    for p in pointers:
        ptr = load_pointer(p.parent)
        assert ptr.url_prefix.startswith("https://"), p
        assert ptr.url_prefix.endswith("/"), p
        assert len(ptr.git_sha) == 40, p
        assert ptr.artifacts, f"{p} lists no artifacts"
        for name, entry in ptr.artifacts.items():
            assert len(entry.sha256) == 64, f"{p}: {name}"
            assert entry.bytes > 0, f"{p}: {name}"
