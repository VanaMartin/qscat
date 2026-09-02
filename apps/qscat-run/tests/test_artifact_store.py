"""The read side of the artifact store.

A published sweep is a function of (config, commit, image) that costs
minutes-to-hours to recompute, so it lives at `data.qscat.org` rather than in
git. What stays in git is a KB-sized `artifacts.json` pointer, and these tests
cover what that pointer has to guarantee: that a reader gets the bytes the
maintainer published, or a clear error -- never silently different bytes.

They also guard the OTHER half of the bargain, which has no runtime code to
fail loudly: only the expensive OUTPUT moves out. The resolved config and the
manifest are inputs and provenance, they are kilobytes, and a clone with no
network that lacks either cannot say what it would have to re-run.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pytest
from qscat_run.artifact_store import (
    ArtifactPointer,
    ArtifactStoreError,
    ChecksumMismatch,
    MissingPointer,
    fetch,
    load_pointer,
)

_REPO = Path(__file__).resolve().parents[3]

#: Records a published run directory must carry IN GIT rather than fetch. The
#: resolved config is the input the sweep is a function of; the manifest says
#: what produced it. Both are KB-sized, so neither is what the store exists
#: for, and putting them behind a download is what makes an offline clone
#: unable to reproduce a published number.
_TRACKED_RECORDS = ("config.resolved.yaml", "manifest.json")


def _committed_pointers() -> list[Path]:
    return sorted(_REPO.glob("validation/**/artifacts.json"))


def _pointer(tmp_path: Path, payload: bytes = b"energy,sigma\n0.1,2.0\n") -> Path:
    d = tmp_path / "o2-ve"
    d.mkdir(parents=True)
    (d / "artifacts.json").write_text(
        json.dumps(
            {
                "git_sha": "0" * 40,
                "url_prefix": "https://data.qscat.org/o2-ve/",
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


def test_the_url_carries_the_content_digest(tmp_path: Path) -> None:
    """The digest goes in the filename, so the address is a function of the
    bytes. This is what makes the scheme hold still: an earlier version put
    the producing commit in the path and every rebase orphaned it."""
    payload = b"energy,sigma\n0.1,2.0\n"
    p = load_pointer(_pointer(tmp_path, payload))
    assert isinstance(p, ArtifactPointer)
    assert p.git_sha == "0" * 40
    digest12 = hashlib.sha256(payload).hexdigest()[:12]
    assert p.url_for("cross_section.csv") == (
        f"https://data.qscat.org/o2-ve/cross_section.{digest12}.csv"
    )


def test_the_extension_survives_so_a_browser_still_knows_what_it_got(tmp_path: Path) -> None:
    url = load_pointer(_pointer(tmp_path)).url_for("cross_section.csv")
    assert url.endswith(".csv")


def test_identical_bytes_address_identically_whatever_produced_them(tmp_path: Path) -> None:
    """Republishing a reproducible run is a no-op rather than a second copy.
    Measured on the real bucket before this change: publishing one experiment
    from a branch and then from main left 25% of it as exact duplicates,
    because the commit differed and the content did not."""
    payload = b"energy,sigma\n0.1,2.0\n"
    a = load_pointer(_pointer(tmp_path / "run-a", payload))
    b = load_pointer(_pointer(tmp_path / "run-b", payload))
    assert a.url_for("cross_section.csv") == b.url_for("cross_section.csv")


def test_a_name_without_an_extension_still_gets_its_digest(tmp_path: Path) -> None:
    d = tmp_path / "run"
    d.mkdir(parents=True)
    payload = b"x"
    sha = hashlib.sha256(payload).hexdigest()
    (d / "artifacts.json").write_text(
        json.dumps(
            {
                "git_sha": "0" * 40,
                "url_prefix": "https://data.qscat.org/run/",
                "artifacts": {"LICENSE": {"sha256": sha, "bytes": 1}},
            }
        )
    )
    assert load_pointer(d).url_for("LICENSE") == f"https://data.qscat.org/run/LICENSE.{sha[:12]}"


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


def test_every_committed_pointer_in_this_repo_is_wellformed() -> None:
    """Guards the pointers themselves: a typo in a URL prefix or a truncated
    digest is invisible until someone tries to fetch, which may be months
    later and on someone else's machine.

    Skips rather than passing silently while no results have been migrated --
    a guard with nothing to guard should say so, not report green.
    """
    pointers = _committed_pointers()
    if not pointers:
        pytest.skip("no artifacts.json committed yet -- nothing has been migrated")
    for p in pointers:
        ptr = load_pointer(p.parent)
        assert ptr.url_prefix.startswith("https://"), p
        assert ptr.url_prefix.endswith("/"), p
        assert len(ptr.git_sha) == 40, p
        assert ptr.artifacts, f"{p} lists no artifacts"
        for name, entry in ptr.artifacts.items():
            assert len(entry.sha256) == 64, f"{p}: {name}"
            assert entry.bytes > 0, f"{p}: {name}"


def test_a_published_run_keeps_its_inputs_beside_the_pointer() -> None:
    """The store holds expensive output; the inputs stay in the clone.

    A reader with no network still has to be able to say what the published
    numbers are a function of, and to re-run it. That takes two kilobyte-sized
    records -- the resolved config and the manifest -- and if either is only
    reachable over HTTPS, the offline half of the design is a claim rather
    than a property. Reading the pointer is not enough on its own: the file
    has to be THERE, and it has to not be listed as something to download.
    """
    pointers = _committed_pointers()
    if not pointers:
        pytest.skip("no artifacts.json committed yet -- nothing has been migrated")
    for p in pointers:
        directory = p.parent
        artifacts = load_pointer(directory).artifacts
        for name in _TRACKED_RECORDS:
            assert (directory / name).is_file(), (
                f"{directory} publishes artifacts but has no {name}; a clone "
                "with no network cannot say what produced them or how to re-run it"
            )
            assert name not in artifacts, (
                f"{p} lists {name} as a fetch-only artifact. It belongs in git: "
                "it is kilobytes, it is an input rather than an output, and "
                "behind a download it is unavailable exactly when it is needed"
            )


def test_the_records_a_published_run_ships_are_tracked_not_merely_present() -> None:
    """Present on a developer's disk is not the same as present in a clone.

    These directories are ignored wholesale and re-populated by `qscat-run
    fetch`, so a resolved config sitting in a working tree looks identical
    whether it was committed or downloaded five minutes ago. Only `git
    ls-files` can tell the difference, and that difference is the whole
    invariant -- the fetched copy is exactly what a fresh clone will not have.
    """
    pointers = _committed_pointers()
    if not pointers:
        pytest.skip("no artifacts.json committed yet -- nothing has been migrated")
    if not (_REPO / ".git").exists():
        # An unpacked tarball or the Docker build context, which excludes
        # `.git`. Nothing to read; the pointer-content half above still runs.
        pytest.skip("no git repository here -- cannot ask what is tracked")
    tracked = set(
        subprocess.run(
            ["git", "ls-files", "-z", "--", "validation"],
            cwd=_REPO,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.split("\0")
    )
    for p in pointers:
        directory = p.parent.relative_to(_REPO)
        for name in _TRACKED_RECORDS:
            assert f"{directory}/{name}" in tracked, (
                f"{directory}/{name} exists but is not tracked -- the .gitignore "
                "allow-list for this directory has to name it, or it reaches "
                "nobody who clones"
            )
