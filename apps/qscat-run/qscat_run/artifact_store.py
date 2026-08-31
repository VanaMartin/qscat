"""Fetching published run artifacts that are too large to keep in git.

A sweep result is a function of (config, commit, container image): it can
always be recomputed, but recomputing costs minutes to hours, so the bytes
live in public object storage and only a pointer is committed.

The pointer is a KB-sized `artifacts.json` written by the maintainer at
publication time, next to the run's `config.resolved.yaml` and `manifest.json`:

    {
      "git_sha": "c884f51...",
      "url_prefix": "https://data.qscat.org/main/c884f51/o2-ve/",
      "artifacts": {"cross_section.csv": {"sha256": "...", "bytes": 399052}}
    }

Reads are anonymous HTTPS -- no account, no credentials, no client library.
Writes are not possible over that hostname at all; publishing goes through the
S3 API with a maintainer token, from the private qscat-infra repository.

Two properties this module exists to guarantee:

* **The bytes are the published bytes.** Every file is checked against the
  recorded digest, and a mismatch raises rather than leaving a plausible file
  on disk. Without this the pointer would only be a hint, and a truncated
  download would read as data.
* **The URL is immutable.** Published paths are never overwritten (the
  publisher refuses, and has no `--force`), so a pointer committed today
  resolves to the same bytes indefinitely. Corrected numbers arrive as a new
  commit with a new pointer, and the old URL keeps the old values -- which is
  what a note that cites a number needs.

The repository still stands alone in the sense that matters: the *inputs* a
test or a claim depends on stay in git, and `config.resolved.yaml` records how
to regenerate any fetched artifact from scratch. What moves out is only the
expensive, reproducible output.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

__all__ = [
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


@dataclass(frozen=True)
class ArtifactEntry:
    sha256: str
    bytes: int


@dataclass(frozen=True)
class ArtifactPointer:
    git_sha: str
    url_prefix: str
    artifacts: dict[str, ArtifactEntry]

    def url_for(self, name: str) -> str:
        return f"{self.url_prefix}{name}"


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
    raw = json.loads(path.read_text())
    return ArtifactPointer(
        git_sha=raw["git_sha"],
        url_prefix=raw["url_prefix"],
        artifacts={
            name: ArtifactEntry(sha256=entry["sha256"], bytes=entry["bytes"])
            for name, entry in raw["artifacts"].items()
        },
    )


def _download(url: str) -> bytes:
    # The URL comes out of a file, and `urlopen` speaks file:// and ftp:// as
    # readily as https. A pointer is reviewed like any other committed file,
    # but "reviewed" is not "validated", and reading a local path because a
    # JSON field asked us to is never what a fetch means.
    if not url.startswith("https://"):
        raise ArtifactStoreError(f"refusing to fetch a non-https URL: {url!r}")
    request = urllib.request.Request(url, headers={"User-Agent": _user_agent()})
    with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
        return bytes(response.read())


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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

    written: list[Path] = []
    for name in wanted:
        entry = pointer.artifacts.get(name)
        if entry is None:
            raise ArtifactStoreError(
                f"{name!r} is not listed in {pointer_path(directory)}; "
                f"published names are {sorted(pointer.artifacts)}"
            )
        target = directory / name
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
