"""Read an immutable, fetched upstream main without touching a checkout's files."""

from __future__ import annotations

import fnmatch
import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote, urlsplit

from tools.agent_search.chunks import sha256

PROFILE = ".opencode/search/profiles.json"


def git(root: Path, *args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(root), *args])


def repository_id(url: str, root: Path) -> str:
    """Share SSH/HTTPS mounts of one upstream, without retaining credentials."""
    scp = re.fullmatch(r"(?:[^/@:]+@)?([^/:]+):(.+)", url)
    if "://" not in url and scp:
        host, path = scp.groups()
        host = host.lower()
    elif "://" in url:
        parts = urlsplit(url)
        if parts.scheme == "file":
            if parts.hostname not in {None, "localhost"}:
                raise ValueError("File upstream must be local")
            return "file:" + str(Path(unquote(parts.path)).resolve())
        host = (parts.hostname or "").lower()
        if not host or parts.scheme not in {"ssh", "https", "http", "git"}:
            raise ValueError("Unsupported upstream URL")
        port = parts.port
        default_port = {"ssh": 22, "https": 443, "http": 80, "git": 9418}[parts.scheme]
        if port and port != default_port:
            host += f":{port}"
        path = parts.path
    else:
        return "file:" + str((root / url).resolve())
    path = path.strip("/").removesuffix(".git")
    if host == "github.com":
        path = path.lower()
    return f"{host}/{path}"


def read_blobs(root: Path, oids: list[str]) -> dict[str, bytes]:
    """Batch object reads; paths, filters, and dirty files never supply source bytes."""
    if not oids:
        return {}
    data = subprocess.check_output(
        ["git", "-C", str(root), "cat-file", "--batch"],
        input="".join(oid + "\n" for oid in oids).encode(),
    )
    result = {}
    offset = 0
    for expected in oids:
        end = data.index(b"\n", offset)
        oid, kind, length = data[offset:end].split()
        if oid.decode() != expected or kind != b"blob":
            raise ValueError("Unexpected Git object in source snapshot")
        length = int(length)
        result[expected] = data[end + 1 : end + 1 + length]
        offset = end + length + 2
    return result


@dataclass(frozen=True)
class MainSource:
    root: Path
    remote: str
    profile_path: str
    repo_id: str
    corpus: str = "repository"

    @classmethod
    def discover(
        cls,
        root: Path,
        remote: str = "origin",
        profile: str = PROFILE,
        *,
        corpus: str = "repository",
    ) -> MainSource:
        root = Path(git(root, "rev-parse", "--show-toplevel").decode().strip()).resolve()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", remote):
            raise ValueError("Remote must be a configured Git remote name")
        path = Path(profile)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Profile must be a repository-relative committed path")
        if corpus not in {"repository", "articles"}:
            raise ValueError("Corpus must be repository or articles")
        url = git(root, "remote", "get-url", remote).decode().strip()
        return cls(root, remote, path.as_posix(), repository_id(url, root), corpus)

    @property
    def owner(self) -> dict:
        owner = {"repo_id": self.repo_id, "branch": "main", "profile_path": self.profile_path}
        if self.corpus != "repository":
            owner["corpus"] = self.corpus
        return owner

    @property
    def table(self) -> str:
        key = json.dumps(self.owner, sort_keys=True)
        prefix = "qscat_articles" if self.corpus == "articles" else "qscat_knowledge"
        return f"{prefix}_main_{sha256(key)[:16]}"

    @property
    def ref(self) -> str:
        return f"refs/remotes/{self.remote}/main"

    def fetch(self) -> Snapshot:
        url = git(self.root, "remote", "get-url", self.remote).decode().strip()
        if repository_id(url, self.root) != self.repo_id:
            raise ValueError("Configured upstream changed since mount")
        env = {**os.environ, "GIT_TERMINAL_PROMPT": "0"}
        env.setdefault("GIT_SSH_COMMAND", "ssh -oBatchMode=yes -oConnectTimeout=10")
        result = subprocess.run(
            [
                "git",
                "-C",
                str(self.root),
                "fetch",
                "--no-tags",
                self.remote,
                f"+refs/heads/main:{self.ref}",
            ],
            env=env,
            capture_output=True,
            timeout=60,
        )
        if result.returncode:
            # Git errors may contain credential-bearing URLs. Keep those on the
            # Git side of the boundary instead of returning them through MCP.
            raise RuntimeError(f"Cannot fetch {self.remote}/main (Git exit {result.returncode})")
        url = git(self.root, "remote", "get-url", self.remote).decode().strip()
        if repository_id(url, self.root) != self.repo_id:
            raise ValueError("Configured upstream changed during fetch")
        commit = git(self.root, "rev-parse", f"{self.ref}^{{commit}}").decode().strip()
        tree = git(self.root, "rev-parse", f"{commit}^{{tree}}").decode().strip()
        blobs = {}
        for entry in git(self.root, "ls-tree", "-r", "-z", commit).split(b"\0"):
            if not entry:
                continue
            metadata, path = entry.split(b"\t", 1)
            mode, kind, oid = metadata.split()
            if kind == b"blob" and mode in {b"100644", b"100755"}:
                blobs[path.decode("utf-8")] = oid.decode()
        if self.profile_path not in blobs:
            raise ValueError("Indexing profile is absent from fetched main")
        policy = read_blobs(self.root, [blobs[self.profile_path]])[blobs[self.profile_path]]
        return Snapshot(self, commit, tree, blobs, json.loads(policy))


@dataclass(frozen=True)
class Snapshot:
    source: MainSource
    commit: str
    tree: str
    blobs: dict[str, str]
    profile: dict

    def selected(self) -> dict[str, str]:
        policy = self.profile[self.source.corpus]
        articles = self.source.corpus == "articles"
        extensions = policy.get("extensions", [".md"]) if articles else policy["extensions"]
        excluded = policy.get("exclude", []) if articles else policy["exclude"]

        def matches(path, patterns):
            return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)

        return {
            path: oid
            for path, oid in sorted(self.blobs.items())
            if matches(path, policy["include"])
            and not matches(path, excluded)
            and (
                Path(path).suffix in extensions
                or matches(path, policy.get("extensionless_files", []))
            )
            and (not articles or Path(path).name != "README.md")
        }
