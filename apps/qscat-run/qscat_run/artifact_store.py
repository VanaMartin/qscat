"""Fetching published run artifacts that are too large to keep in git.

A sweep result is a function of (config, commit, container image): it can
always be recomputed, but recomputing costs minutes to hours, so the bytes
live in public object storage and only a pointer is committed.

The pointer is a KB-sized `artifacts.json` written by the maintainer at
publication time, next to the run's `config.resolved.yaml` and `manifest.json`:

    {
      "git_sha": "69742d8...",
      "url_prefix": "https://data.qscat.org/o2-ve/",
      "artifacts": {"cross_section.csv": {"sha256": "830cffb8...", "bytes": 399052}}
    }

Objects are addressed BY CONTENT: the digest goes into the filename, so
`cross_section.csv` above lives at

    https://data.qscat.org/o2-ve/cross_section.830cffb8a044.csv

The folder is a human label, not an identity -- it keeps a URL readable
enough to paste into a paper, while the hash decides what the URL means.

Addressing by content rather than by commit is what makes the scheme hold
still. An earlier version put the producing commit in the path, and every
rebase orphaned it: three times, because a branch commit is not a stable
address. It also stored the same bytes repeatedly -- publishing one
experiment from a branch and then from main left 25% of the bucket as exact
duplicates, since the commit changed and the content did not. Under content
addressing a re-run that reproduces its numbers republishes to the same key
and costs nothing, which is the common case here: the O2 sweeps came out
bit-identical on all three runs.

Two properties this module guarantees:

* **The bytes are the published bytes.** Every file is checked against the
  recorded digest, and a mismatch raises rather than leaving a plausible file
  on disk. Without this the pointer would only be a hint, and a truncated
  download would read as data.
* **A URL means one thing.** Not by policy -- by construction. Different
  content hashes differently and therefore lives elsewhere, so a link in a
  note cannot quietly come to mean something else. Correcting a number
  produces a new key and a new pointer; the old URL keeps the old value for
  as long as it is kept, which is what a cited number needs.

`git_sha` stays in the manifest, but it is a RECORD, not an address. It
answers "what produced these bytes"; `git blame` on the pointer answers "when
did this repository start citing them". Those are different questions and the
first one deserves an answer that does not depend on history staying still.

The repository still stands alone in the sense that matters: the *inputs* a
test or a claim depends on stay in git, and `config.resolved.yaml` records how
to regenerate any fetched artifact from scratch. What moves out is only the
expensive, reproducible output.

So a pointer lists only files that are NOT in git. `config.resolved.yaml` and
`manifest.json` sit beside it, tracked, and never appear among its artifacts:
they are kilobytes, they are the input and the provenance rather than the
result, and behind a download they would be missing exactly when someone
without a network wants to know what a published number came from. Nothing
here can enforce that at runtime -- a directory whose config arrived by fetch
is indistinguishable from one whose config was cloned -- so it is enforced by
a guard over the committed pointers instead.

A pointer is UNTRUSTED STRUCTURED INPUT, and this module treats it as one.
It is a JSON file that decides, on the reader's machine, which host is
contacted and which paths are written; a committed pointer is reviewed, but
review is not validation, and a pointer can also reach a reader by any route
a file reaches a reader. So every field is checked before a byte is requested
or written -- names, digests, byte counts and the URL prefix alike -- and the
destination of each artifact is resolved and proved to lie inside the
requested directory first. The checksum is NOT that boundary: it says the
bytes are the published bytes and says nothing about where they land, and a
file written outside the run directory has already done its damage by the
time it verifies.
"""

from __future__ import annotations

import hashlib
import json
import re
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from urllib.parse import SplitResult, urlsplit

__all__ = [
    "ARTIFACT_HOSTS",
    "POINTER_NAME",
    "ArtifactEntry",
    "ArtifactPointer",
    "ArtifactStoreError",
    "ChecksumMismatch",
    "MissingPointer",
    "fetch",
    "load_pointer",
    "pointer_path",
]

POINTER_NAME = "artifacts.json"

#: Hostnames a pointer may name. The bucket is published read-only at exactly
#: one address (see `docs/adr/0008-computed-artifacts-live-in-public-object-
#: storage.md`), and every pointer this tool reads is a committed file in this
#: repository, so "any HTTPS host" would not be a policy anyone chose -- it
#: would be the absence of one. Matching is EXACT and case-insensitive, never
#: by suffix: suffix matching is what makes `data.qscat.org.example.net` look
#: like the store.
#:
#: A fork that publishes its results elsewhere changes this line. That is the
#: intended cost: pointing a fetch at a new host is a reviewable edit to code,
#: not something a data file can decide. Tests never need it -- `fetch` takes
#: an `opener`, which replaces the network outright.
ARTIFACT_HOSTS = frozenset({"data.qscat.org"})

#: One path component of an artifact name. Deliberately narrow: the name is
#: BOTH a path below the run directory and a path segment of the URL, pasted
#: in unescaped, so it is restricted to characters that need no escaping in
#: either. Anything outside this -- a separator, a space, a percent sign, a
#: NUL, a leading dot -- makes the file on disk and the object in the store
#: two different names, which is exactly the confusion an escape exploits.
_NAME_COMPONENT = re.compile(r"[A-Za-z0-9][A-Za-z0-9._+~-]*")

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_SHA = re.compile(r"[0-9a-f]{40}")

#: Read timeout in seconds. Generous: some published sweeps are megabytes and
#: the reader may be far from the edge node.
_TIMEOUT = 120


def _user_agent() -> str:
    """Identify this tool by name and version.

    LOAD-BEARING, not politeness. `urllib` defaults to `Python-urllib/3.x`,
    which Cloudflare answers with **403** in front of the artifact bucket --
    measured: `curl` sent with that agent is refused too, and `urllib` sent
    with curl's is served, so it is the agent and not the object or the client.
    Left at the default, every `qscat-run fetch` fails for every reader.

    The value says what it is and where to complain, rather than impersonating
    a browser: an operator reading the logs should be able to tell who called.
    """
    try:
        from importlib import metadata

        version = metadata.version("qscat-run")
    except Exception:
        version = "unknown"
    return f"qscat-run/{version} (+https://github.com/VanaMartin/qscat)"


class ArtifactStoreError(RuntimeError):
    """Base class for artifact-store failures."""


class MissingPointer(ArtifactStoreError):
    """No `artifacts.json` in the directory, so nothing says what to fetch."""


class ChecksumMismatch(ArtifactStoreError):
    """Downloaded bytes do not match the digest the maintainer published."""


def _check_store_url(url: object, what: str) -> SplitResult:
    """Parse `url` and refuse anything that is not the published store.

    Parsed, never matched as a string. `https://evil.example/@data.qscat.org/`
    and `https://data.qscat.org@evil.example/` both begin with `https://` and
    both mention the right hostname; both connect somewhere else. Only a
    parser knows which host a URL names, and `hostname` is the field that
    answers it -- `netloc` still carries userinfo and a port.
    """
    if not isinstance(url, str):
        raise ArtifactStoreError(f"{what} must be a string, got {type(url).__name__}")
    parts = urlsplit(url)
    if parts.scheme != "https":
        # `urlopen` speaks file:// and ftp:// as readily as https, and reading
        # a local path because a JSON field asked us to is never what a fetch
        # means.
        raise ArtifactStoreError(f"refusing a non-https {what}: {url!r}")
    if parts.username is not None or parts.password is not None:
        raise ArtifactStoreError(
            f"refusing a {what} carrying credentials, which also hide the real host: {url!r}"
        )
    host = parts.hostname
    if host is None or host.lower() not in ARTIFACT_HOSTS:
        raise ArtifactStoreError(
            f"refusing a {what} pointing at {host!r} rather than "
            f"{'/'.join(sorted(ARTIFACT_HOSTS))}: {url!r}"
        )
    if parts.port is not None:
        raise ArtifactStoreError(
            f"refusing a {what} that names a port; the store is read over "
            f"HTTPS on its default port: {url!r}"
        )
    if parts.query or parts.fragment:
        raise ArtifactStoreError(
            f"refusing a {what} with a query or fragment -- the store addresses "
            f"objects by path alone: {url!r}"
        )
    return parts


def _check_url_prefix(prefix: object) -> None:
    """The prefix an artifact name is concatenated onto.

    Concatenation is why the trailing slash is required rather than repaired:
    without it `.../o2-ve` + `cross_section...` silently names a different
    folder, and a fetch that fails with 404 is the good case.
    """
    parts = _check_store_url(prefix, "artifact URL prefix")
    path = parts.path
    if not path.startswith("/") or not path.endswith("/"):
        raise ArtifactStoreError(
            f"an artifact URL prefix must end in '/'; names are appended to it: {prefix!r}"
        )
    segments = path[1:-1].split("/") if len(path) > 1 else []
    if any(segment in ("", ".", "..") for segment in segments):
        raise ArtifactStoreError(f"refusing an artifact URL prefix with a dot segment: {prefix!r}")


def _check_name(name: object) -> None:
    """The name of one artifact: a relative path to a file below the run directory.

    NESTED NAMES ARE SUPPORTED, and they have to be: a run writes its
    wavefunction snapshots, eigenstates and resonance states into
    subdirectories of its output directory, so a pointer published for such a
    run names `wavefunction/psi_E0.05.npz`. What is refused is every spelling
    that is not a relative path to a plain file below the directory -- the
    empty name, an absolute one, the dot names, a Windows separator (which a
    POSIX `Path` would read as one innocent component), and any character
    that would have to be escaped in a URL.

    Rejecting per component rather than normalising is the point: `..` never
    reaches a path join, so there is no arithmetic left to get wrong, and the
    name that goes into the URL is the name that goes onto the disk.
    """
    if not isinstance(name, str) or not name:
        raise ArtifactStoreError(f"an artifact name must be a non-empty string, got {name!r}")
    if "\\" in name:
        raise ArtifactStoreError(f"refusing an artifact name with a backslash: {name!r}")
    for component in name.split("/"):
        if not _NAME_COMPONENT.fullmatch(component):
            raise ArtifactStoreError(
                f"refusing the artifact name {name!r}: {component!r} is not a plain "
                "file or directory name. Names are relative paths below the run "
                "directory, built from letters, digits and '. _ + ~ -', and may not "
                "begin with a dot."
            )


@dataclass(frozen=True)
class ArtifactEntry:
    sha256: str
    bytes: int

    def __post_init__(self) -> None:
        if not isinstance(self.sha256, str) or not _SHA256.fullmatch(self.sha256):
            raise ArtifactStoreError(
                f"{self.sha256!r} is not a sha256: 64 lower-case hex characters are "
                "required. The first twelve of them ARE the object key, so a "
                "malformed digest does not merely fail to verify later -- it "
                "addresses the wrong object now."
            )
        # `bool` is a subclass of `int`, so the obvious isinstance check reads
        # `true` as a one-byte artifact.
        if isinstance(self.bytes, bool) or not isinstance(self.bytes, int) or self.bytes < 0:
            raise ArtifactStoreError(f"{self.bytes!r} is not a byte count")


#: Hex characters of the sha256 that go into a published filename. 12 gives a
#: ~1e-18 collision chance across any plausible number of artifacts, and keeps
#: the URL short enough to read aloud. The pointer always carries the FULL
#: digest, which is what the download is verified against -- the truncation
#: only names the object.
_URL_DIGEST_CHARS = 12


@dataclass(frozen=True)
class ArtifactPointer:
    """A validated pointer: no instance of this exists that names a bad target.

    The checks live here rather than in `load_pointer` so that the guarantee
    belongs to the TYPE. Everything downstream -- URL derivation, the fetch
    loop, `--list` -- may then read the fields without re-deciding whether
    they are safe, and a second construction path could not skip the checks.
    """

    git_sha: str
    url_prefix: str
    artifacts: dict[str, ArtifactEntry]

    def __post_init__(self) -> None:
        if not isinstance(self.git_sha, str) or not _GIT_SHA.fullmatch(self.git_sha):
            raise ArtifactStoreError(
                f"{self.git_sha!r} is not a git sha: 40 lower-case hex characters "
                "are required. It is a record of what produced the bytes, not an "
                "address for them, but a pointer that cannot name its commit was "
                "not written by the publisher."
            )
        _check_url_prefix(self.url_prefix)
        for name in self.artifacts:
            _check_name(name)

    def url_for(self, name: str) -> str:
        """The content-addressed URL for `name`.

        Derived rather than stored: the digest is already in the pointer, so a
        second copy of the address could only ever disagree with it.

        The digest goes before the extension of the FILE. Splitting on the
        last dot of the whole name would put it in the wrong place for a
        nested artifact whose DIRECTORY carries the only dot.
        """
        entry = self.artifacts[name]
        suffix = PurePosixPath(name).suffix  # "" when the file has no extension
        stem = name[: len(name) - len(suffix)]
        digest = entry.sha256[:_URL_DIGEST_CHARS]
        return f"{self.url_prefix}{stem}.{digest}{suffix}"


def pointer_path(directory: str | Path) -> Path:
    return Path(directory) / POINTER_NAME


def load_pointer(directory: str | Path) -> ArtifactPointer:
    path = pointer_path(directory)
    if not path.is_file():
        raise MissingPointer(
            f"{path} not found. A run directory that keeps its artifacts in the "
            f"store carries a {POINTER_NAME}; one that keeps them in git does not "
            "need fetching."
        )
    try:
        raw = json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ArtifactStoreError(f"{path} is not readable as JSON: {exc}") from exc
    try:
        return _pointer_from(raw)
    except ArtifactStoreError as exc:
        # Named, so the reader knows WHICH pointer to open. `fetch` may be
        # given several directories.
        raise ArtifactStoreError(f"{path}: {exc}") from exc


def _pointer_from(raw: object) -> ArtifactPointer:
    """Build a pointer out of decoded JSON, rejecting anything that is not one.

    Every failure is an `ArtifactStoreError`: the CLI turns those into a
    message and anything else into a traceback, and a hand-edited pointer is
    a bad file rather than a bug in the fetcher.
    """
    if not isinstance(raw, dict):
        raise ArtifactStoreError(f"expected a JSON object, got {type(raw).__name__}")
    missing = [key for key in ("git_sha", "url_prefix", "artifacts") if key not in raw]
    if missing:
        raise ArtifactStoreError(f"missing {', '.join(missing)}")
    artifacts = raw["artifacts"]
    if not isinstance(artifacts, dict):
        raise ArtifactStoreError(f"'artifacts' must be an object, got {type(artifacts).__name__}")
    entries: dict[str, ArtifactEntry] = {}
    for name, entry in artifacts.items():
        if not isinstance(entry, dict) or not {"sha256", "bytes"} <= set(entry):
            raise ArtifactStoreError(f"{name!r} needs an object with 'sha256' and 'bytes'")
        entries[name] = ArtifactEntry(sha256=entry["sha256"], bytes=entry["bytes"])
    return ArtifactPointer(git_sha=raw["git_sha"], url_prefix=raw["url_prefix"], artifacts=entries)


def _download(url: str) -> bytes:
    # Checked again here, on the function that actually opens a socket, so it
    # does not depend on its caller having checked. Cheap, and it is the one
    # place where getting it wrong reaches the network.
    _check_store_url(url, "URL")
    request = urllib.request.Request(url, headers={"User-Agent": _user_agent()})
    with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
        return bytes(response.read())


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _contained_target(base: Path, name: str) -> Path:
    """Where `name` will be written, proved to be inside `base`.

    `base` must already be resolved, and BOTH sides have to be: on macOS
    `/tmp` is a symlink to `/private/tmp`, so comparing a resolved target
    against an unresolved base rejects every fetch under a temporary
    directory. `Path.resolve()` is non-strict, so it works on a file that does
    not exist yet -- it normalises the path and resolves the symlinks in
    whatever part of it does exist, on macOS and Linux alike.

    That existing part is the reason this is not merely a re-check of the
    name. A name may be a perfectly plain filename and the destination still
    be somewhere else, because a file or directory already in the run
    directory is a symlink pointing out of it; the write would follow it.
    What must be contained is what the write would open, not what the name
    looks like.
    """
    target = (base / name).resolve()
    if target == base or not target.is_relative_to(base):
        raise ArtifactStoreError(f"refusing to write {name!r} to {target}, which is outside {base}")
    return target


def fetch(
    directory: str | Path,
    *,
    names: list[str] | None = None,
    opener: Callable[[str], bytes] = _download,
) -> list[Path]:
    """Download the artifacts `directory` points at, and return what was written.

    A file already present with the right digest is left alone and not
    returned -- so a second call is free, and an interrupted fetch resumes
    rather than restarting. A file present with the WRONG digest is replaced:
    a half-written or locally edited copy is not a cache hit.

    `opener` exists so the verification logic can be tested without a network.
    """
    directory = Path(directory)
    pointer = load_pointer(directory)
    wanted = names if names is not None else sorted(pointer.artifacts)

    # Every destination is proved before the first download, not one at a
    # time: a fetch that delivered half a run and then escaped would have
    # escaped. Loading the pointer has already refused every name that cannot
    # be one; this is the second half of that, the half only the filesystem
    # can answer.
    base = directory.resolve()
    targets: dict[str, Path] = {}
    for name in wanted:
        if name not in pointer.artifacts:
            raise ArtifactStoreError(
                f"{name!r} is not listed in {pointer_path(directory)}; "
                f"published names are {sorted(pointer.artifacts)}"
            )
        targets[name] = _contained_target(base, name)

    written: list[Path] = []
    for name in wanted:
        entry = pointer.artifacts[name]
        target = targets[name]
        if target.is_file() and _digest(target.read_bytes()) == entry.sha256:
            continue

        url = pointer.url_for(name)
        data = opener(url)
        got = _digest(data)
        if got != entry.sha256:
            # Deliberately NOT written to disk: a file that failed its check
            # must not be left where a later run would treat it as data.
            raise ChecksumMismatch(
                f"{url}\n  expected sha256 {entry.sha256}\n  got      sha256 {got}\n"
                "The published object should be immutable, so this means a "
                "truncated download or a corrupted pointer -- not a newer version."
            )
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        written.append(target)

    return written
