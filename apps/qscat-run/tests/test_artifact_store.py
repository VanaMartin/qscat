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
from urllib.parse import urlsplit

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


# --- The pointer as untrusted input -----------------------------------------
#
# `artifacts.json` is a file: reviewed when it is committed, but review is not
# validation, and a fetch acts on whatever it says. Everything below is an
# attack on that -- a name that walks out of the run directory, a digest that
# cannot be a digest, a URL that only LOOKS like the published store -- and
# every one of them has to fail before a byte is requested or written.

_ANY_DIGEST = "0" * 64


def _write_pointer(
    directory: Path,
    artifacts: dict[str, object],
    *,
    url_prefix: str = "https://data.qscat.org/o2-ve/",
    git_sha: str = "0" * 40,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "artifacts.json").write_text(
        json.dumps({"git_sha": git_sha, "url_prefix": url_prefix, "artifacts": artifacts})
    )
    return directory


def _explode(url: str) -> bytes:
    raise AssertionError(f"a rejected pointer still reached the network: {url}")


#: Names that must never become a destination. Traversal in both separators,
#: because a POSIX `Path` reads `..\..\x` as one harmless component and a
#: publisher on another platform does not; the dot names, which address a
#: directory rather than a file; an absolute path; the empty name; and forms
#: that would have to be escaped to appear in a URL, where the escaped and the
#: unescaped spelling name different things.
_REFUSED_NAMES = [
    "../escape.csv",
    "../../escape.csv",
    "sub/../../escape.csv",
    "/etc/passwd",
    "",
    ".",
    "..",
    "..\\..\\escape.csv",
    "sub\\escape.csv",
    "cross section.csv",
    "cross%2e%2e%2fsection.csv",
    "escape.csv\x00.png",
    "sub//escape.csv",
    ".hidden.csv",
]


@pytest.mark.parametrize("name", _REFUSED_NAMES)
def test_a_name_that_is_not_a_plain_file_below_the_directory_is_refused(
    tmp_path: Path, name: str
) -> None:
    """Reading the pointer is where this has to fail: `--list` derives a URL
    from every name, so a name that cannot be a destination must not survive
    long enough to become an address either."""
    d = _write_pointer(tmp_path / "run", {name: {"sha256": _ANY_DIGEST, "bytes": 1}})
    with pytest.raises(ArtifactStoreError):
        load_pointer(d)


@pytest.mark.parametrize("name", _REFUSED_NAMES)
def test_a_refused_name_downloads_nothing_and_writes_nothing(tmp_path: Path, name: str) -> None:
    """The checksum is not the boundary. It says the bytes are the published
    bytes; it says nothing about where they land, and a file written outside
    the run directory has already done its damage by the time it verifies."""
    d = _write_pointer(tmp_path / "run", {name: {"sha256": _ANY_DIGEST, "bytes": 1}})
    with pytest.raises(ArtifactStoreError):
        fetch(d, opener=_explode)
    assert sorted(p.name for p in tmp_path.rglob("*") if p.is_file()) == ["artifacts.json"]


def test_the_traversal_name_this_guards_against_would_otherwise_have_escaped(
    tmp_path: Path,
) -> None:
    """States the defect in one line: `directory / "../../escape.csv"` is a
    path OUTSIDE `directory`, and joining is where that happened -- not
    downloading, and not verifying."""
    assert (tmp_path / "run" / "../../escape.csv").resolve() == (
        tmp_path.parent / "escape.csv"
    ).resolve()


@pytest.mark.parametrize(
    "digest", ["0" * 63, "0" * 65, "g" * 64, "", "830cffb8a044", "A" * 64, 12345, None]
)
def test_a_field_that_cannot_be_a_sha256_is_refused(tmp_path: Path, digest: object) -> None:
    """A digest is what makes a URL mean one thing, and its first twelve
    characters ARE the object key -- so a short, non-hex or upper-case digest
    does not merely fail to verify later, it addresses the wrong object now."""
    d = _write_pointer(tmp_path / "run", {"cross_section.csv": {"sha256": digest, "bytes": 1}})
    with pytest.raises(ArtifactStoreError):
        load_pointer(d)


@pytest.mark.parametrize("count", [-1, 1.5, "399052", True, None])
def test_a_byte_count_that_is_not_a_count_is_refused(tmp_path: Path, count: object) -> None:
    """`True` is in the list on purpose: `bool` is a subclass of `int`, so the
    obvious `isinstance(n, int)` accepts it and a pointer would record an
    artifact of one byte."""
    d = _write_pointer(
        tmp_path / "run", {"cross_section.csv": {"sha256": _ANY_DIGEST, "bytes": count}}
    )
    with pytest.raises(ArtifactStoreError):
        load_pointer(d)


#: URL prefixes that must not be fetched from. The first two are other
#: protocols. Most of the rest begin with the eight characters `https://`,
#: which is why a string prefix check is not a check: only parsing says which
#: host the connection would actually go to.
_REFUSED_PREFIXES = [
    "http://data.qscat.org/o2-ve/",
    "file:///etc/",
    "https://evil.example/@data.qscat.org/o2-ve/",
    "https://data.qscat.org@evil.example/o2-ve/",
    "https://data.qscat.org.evil.example/o2-ve/",
    "https://notdata.qscat.org/o2-ve/",
    "https:/\\/\\evil.example/o2-ve/",
    "https://data.qscat.org:8443/o2-ve/",
    "https://data.qscat.org/o2-ve/?token=x",
    "https://data.qscat.org/o2-ve/#frag",
    "https://data.qscat.org/../o2-ve/",
    "https://data.qscat.org/o2-ve",
    "https:///o2-ve/",
]


@pytest.mark.parametrize("prefix", _REFUSED_PREFIXES)
def test_a_url_prefix_that_is_not_the_published_store_is_refused(
    tmp_path: Path, prefix: str
) -> None:
    d = _write_pointer(
        tmp_path / "run",
        {"cross_section.csv": {"sha256": _ANY_DIGEST, "bytes": 1}},
        url_prefix=prefix,
    )
    with pytest.raises(ArtifactStoreError):
        load_pointer(d)


def test_the_deceptive_prefixes_all_pass_the_check_they_replace() -> None:
    """Why parsing, and not `startswith`. Each of these reads as HTTPS and
    mentions the right hostname somewhere in the string, and each of them
    would connect somewhere else."""
    for prefix in [
        "https://evil.example/@data.qscat.org/o2-ve/",
        "https://data.qscat.org@evil.example/o2-ve/",
        "https://data.qscat.org.evil.example/o2-ve/",
    ]:
        assert prefix.startswith("https://")


def test_a_refused_prefix_downloads_nothing(tmp_path: Path) -> None:
    d = _write_pointer(
        tmp_path / "run",
        {"cross_section.csv": {"sha256": _ANY_DIGEST, "bytes": 1}},
        url_prefix="https://data.qscat.org.evil.example/o2-ve/",
    )
    with pytest.raises(ArtifactStoreError):
        fetch(d, opener=_explode)


def test_the_host_is_matched_case_insensitively(tmp_path: Path) -> None:
    """Hostnames are case-insensitive, so refusing a capitalised one would be
    a bug rather than a boundary."""
    payload = b"energy,sigma\n0.1,2.0\n"
    d = _write_pointer(
        tmp_path / "run",
        {
            "cross_section.csv": {
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
        },
        url_prefix="https://DATA.QSCAT.ORG/o2-ve/",
    )
    assert fetch(d, opener=lambda url: payload) == [d / "cross_section.csv"]


def test_a_nested_name_lands_in_its_subdirectory(tmp_path: Path) -> None:
    """Nested names are supported: a run writes its wavefunction, eigenstate
    and resonance snapshots into subdirectories of the run directory, so a
    pointer published for such a run has to be able to name them."""
    payload = b"\x93NUMPY-ish"
    d = _write_pointer(
        tmp_path / "run",
        {
            "wavefunction/psi_E0.05.npz": {
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
        },
    )
    written = fetch(d, opener=lambda url: payload)
    assert written == [d / "wavefunction" / "psi_E0.05.npz"]
    assert (d / "wavefunction" / "psi_E0.05.npz").read_bytes() == payload


def test_a_nested_name_keeps_its_directory_in_the_published_url(tmp_path: Path) -> None:
    """The key in the store mirrors the path in the run directory, and the
    digest goes before the extension of the FILE -- not before whatever
    follows the last dot in the whole name."""
    d = _write_pointer(
        tmp_path / "run", {"wavefunction/psi_E0.05.npz": {"sha256": _ANY_DIGEST, "bytes": 1}}
    )
    assert load_pointer(d).url_for("wavefunction/psi_E0.05.npz") == (
        f"https://data.qscat.org/o2-ve/wavefunction/psi_E0.05.{_ANY_DIGEST[:12]}.npz"
    )


def test_fetch_refuses_to_write_through_a_symlink_that_leaves_the_directory(
    tmp_path: Path,
) -> None:
    """A name can be a plain file and the destination still be elsewhere:
    `resolve()` follows symlinks, so the check has to be made on what the
    write would actually open, not on what the name looks like."""
    outside = tmp_path / "outside.csv"
    outside.write_bytes(b"do not touch")
    payload = b"energy,sigma\n0.1,2.0\n"
    d = _write_pointer(
        tmp_path / "run",
        {
            "cross_section.csv": {
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
        },
    )
    (d / "cross_section.csv").symlink_to(outside)
    with pytest.raises(ArtifactStoreError):
        fetch(d, opener=_explode)
    assert outside.read_bytes() == b"do not touch"


def test_one_escaping_target_stops_the_whole_fetch_before_any_of_it(tmp_path: Path) -> None:
    """Containment is proved for every wanted artifact before the first
    download, so a pointer cannot deliver half a run and then escape."""
    outside = tmp_path / "outside.npz"
    outside.write_bytes(b"do not touch")
    payload = b"energy,sigma\n0.1,2.0\n"
    sha = hashlib.sha256(payload).hexdigest()
    d = _write_pointer(
        tmp_path / "run",
        {
            "cross_section.csv": {"sha256": sha, "bytes": len(payload)},
            "cross_section.npz": {"sha256": sha, "bytes": len(payload)},
        },
    )
    (d / "cross_section.npz").symlink_to(outside)
    with pytest.raises(ArtifactStoreError):
        fetch(d, opener=lambda url: payload)
    assert not (d / "cross_section.csv").exists()
    assert outside.read_bytes() == b"do not touch"


@pytest.mark.parametrize(
    "text",
    [
        "{not json",
        "[]",
        '{"url_prefix": "https://data.qscat.org/o2-ve/", "artifacts": {}}',
        '{"git_sha": "0", "url_prefix": "https://data.qscat.org/o2-ve/", "artifacts": {}}',
        '{"git_sha": "' + "0" * 40 + '", "url_prefix": 7, "artifacts": {}}',
        '{"git_sha": "'
        + "0" * 40
        + '", "url_prefix": "https://data.qscat.org/o2-ve/", "artifacts": []}',
        '{"git_sha": "'
        + "0" * 40
        + '", "url_prefix": "https://data.qscat.org/o2-ve/", "artifacts": {"a.csv": 3}}',
    ],
)
def test_a_malformed_pointer_says_so_instead_of_raising_a_decoding_error(
    tmp_path: Path, text: str
) -> None:
    """The CLI turns `ArtifactStoreError` into a message and anything else
    into a traceback, so a hand-edited pointer has to fail as a pointer."""
    d = tmp_path / "run"
    d.mkdir()
    (d / "artifacts.json").write_text(text)
    with pytest.raises(ArtifactStoreError):
        load_pointer(d)


def test_every_committed_pointer_points_at_the_published_store() -> None:
    """The hostname policy, stated over the real pointers: this repository
    fetches from the one read-only bucket `docs/adr/0008` binds, and from
    nowhere else."""
    pointers = _committed_pointers()
    if not pointers:
        pytest.skip("no artifacts.json committed yet -- nothing has been migrated")
    for p in pointers:
        host = urlsplit(load_pointer(p.parent).url_prefix).hostname
        assert host == "data.qscat.org", f"{p} fetches from {host}"
