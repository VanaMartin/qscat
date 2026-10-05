"""Populate and reconcile a worktree-owned LanceDB code corpus.

Run from the repository root with the isolated .opencode/search environment.
Location metadata is kept outside the embedding source field. Existing vectors
are reused by exact payload/model keys, and removed rows are reconciled by ID.
"""

from __future__ import annotations

import argparse
import fcntl
import fnmatch
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

from tools.agent_search.chunks import sha256, split_units

VERSION = "qmodeling-search-v1"
CHUNKER = "python-ast-rust-treesitter-markdown-v2"


def fingerprint(value) -> str:
    return sha256(json.dumps(value, sort_keys=True, ensure_ascii=False))


def git(root: Path, *args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(root), *args])


def matches(path: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)


def selected_files(root: Path, policy: dict, include_new: list[str]) -> list[str]:
    names = {p.decode() for p in git(root, "ls-files", "-z").split(b"\0") if p}
    for name in include_new:
        path = (root / name).resolve()
        names.add(path.relative_to(root).as_posix())
    chosen = []
    for name in sorted(names):
        path = root / name
        if not path.is_file() or path.is_symlink():
            continue
        if not path.resolve().is_relative_to(root):
            raise ValueError(f"Source escapes worktree: {name}")
        if not matches(name, policy["include"]) or matches(name, policy["exclude"]):
            continue
        if path.suffix not in policy["extensions"] and not matches(
            name, policy.get("extensionless_files", [])
        ):
            continue
        chosen.append(name)
    return chosen


def source_kind(path: str) -> str:
    parts = Path(path).parts
    if "tests" in parts or Path(path).name.startswith("test_"):
        return "test"
    if path.endswith((".md", ".rst")):
        return (
            "guidance"
            if path.startswith((".claude/", ".opencode/")) or path in {"CLAUDE.md", "AGENTS.md"}
            else "documentation"
        )
    return "code" if path.endswith((".py", ".rs")) else "configuration"


def make_records(
    root: Path,
    profile: dict,
    files: list[str],
    tokenizer,
    model_key: str,
    *,
    cached_files: dict | None = None,
    cached_rows: list[dict] | None = None,
) -> tuple[list[dict], dict]:
    records = []
    hashes = {}
    worktree = sha256(str(root))
    by_path: dict[str, list[dict]] = {}
    for row in cached_rows or []:
        by_path.setdefault(row["path"], []).append(row)
    for name in files:
        raw = (root / name).read_bytes()
        file_hash = hashlib.sha256(raw).hexdigest()
        hashes[name] = file_hash
        if cached_files is not None and cached_files.get(name) == file_hash:
            records.extend(
                {k: v for k, v in r.items() if k != "vector"} for r in by_path.get(name, [])
            )
            continue
        try:
            language, chunks = split_units(
                name, raw.decode("utf-8"), tokenizer, profile["embedding"]
            )
        except (SyntaxError, UnicodeError, ValueError) as error:
            raise ValueError(f"Cannot index {name}: {error}") from error
        for chunk in chunks:
            row_id = fingerprint(
                [
                    "qscat",
                    name,
                    language,
                    chunk.kind,
                    chunk.anchor,
                    chunk.declaration,
                    chunk.content_sha256,
                    chunk.duplicate,
                ]
            )
            doc = (
                f"source={name}\nanchor={chunk.anchor}\nkind={chunk.kind}\n"
                f"declaration={chunk.declaration}\n"
                f"lines={chunk.start_line}-{chunk.end_line}\nroot={root}\n"
                f"source_sha256={file_hash}\nchunk_id={row_id}\n\n{chunk.content}"
            )
            records.append(
                {
                    "id": row_id,
                    "doc": doc,
                    "text": chunk.text,
                    "path": name,
                    "anchor": chunk.anchor,
                    "kind": chunk.kind,
                    "language": language,
                    "source_kind": source_kind(name),
                    "repo_id": "qscat",
                    "worktree_id": worktree,
                    "start_line": chunk.start_line,
                    "end_line": chunk.end_line,
                    "declaration": chunk.declaration,
                    "source_sha256": file_hash,
                    "content_sha256": chunk.content_sha256,
                    "payload_sha256": chunk.payload_sha256,
                    "embedding_key": fingerprint([model_key, chunk.payload_sha256]),
                    "tokens": chunk.tokens,
                }
            )
    if not records:
        raise ValueError("No searchable chunks selected; refusing to publish an empty corpus")
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Chunk IDs are not unique")
    return records, hashes


def embedding(profile: dict):
    import torch
    from lancedb.embeddings import get_registry

    torch.set_num_threads(2)
    settings = profile["embedding"]
    model = (
        get_registry()
        .get(settings["registry"])
        .create(
            name=settings["model"],
            device=settings["device"],
            normalize=True,
            trust_remote_code=False,
        )
    )
    transformer = model.embedding_model
    revision = getattr(transformer[0].auto_model.config, "_commit_hash", None)
    if not revision:
        raise ValueError("Cannot identify the embedding model revision")
    if transformer.max_seq_length < settings["total_input_max_tokens"]:
        raise ValueError("Profile exceeds the installed model's input limit")
    spec = {
        "model": settings["model"],
        "revision": revision,
        "normalize": True,
        "dimensions": model.ndims(),
        "max_input_tokens": transformer.max_seq_length,
        "precision": "float32",
        "sentence_transformers": version("sentence-transformers"),
        "transformers": version("transformers"),
        "tokenizers": version("tokenizers"),
    }
    if spec["dimensions"] != settings["dimensions"]:
        raise ValueError("Embedding dimensions differ from the profile")
    return model, spec


def schema(model):
    from lancedb.pydantic import LanceModel, Vector
    from pydantic import create_model

    strings = [
        "id",
        "doc",
        "path",
        "anchor",
        "kind",
        "language",
        "source_kind",
        "repo_id",
        "worktree_id",
        "source_sha256",
        "content_sha256",
        "payload_sha256",
        "embedding_key",
    ]
    integers = ["start_line", "end_line", "declaration", "tokens"]
    fields = {name: (str, ...) for name in strings}
    fields.update({name: (int, ...) for name in integers})
    fields["text"] = (str, model.SourceField())
    fields["vector"] = (Vector(model.ndims()), model.VectorField())
    return create_model("CodeChunk", __base__=LanceModel, **fields).to_arrow_schema()


def manifest_path(db_path: Path, table: str) -> Path:
    return db_path / f"{table}.manifest.json"


def read_manifest(db_path: Path, table: str) -> dict | None:
    path = manifest_path(db_path, table)
    return json.loads(path.read_text()) if path.exists() else None


def select_table(root: Path, db_path: Path, base: str = "qscat_knowledge") -> str:
    """Reuse the existing owned base table, otherwise isolate this checkout."""
    manifest = read_manifest(db_path, base)
    if manifest and manifest.get("root") == str(root) and manifest.get("format") == VERSION:
        return base
    return f"{base}_{sha256(str(root))[:16]}"


def check_owner(manifest: dict | None, root: Path) -> None:
    if manifest and (manifest.get("format") != VERSION or manifest.get("root") != str(root)):
        raise ValueError(
            "Table belongs to a different worktree or index format; use a separate table"
        )


@contextmanager
def writer_lock(db_path: Path, table: str):
    db_path.mkdir(parents=True, exist_ok=True)
    with (db_path / f"{table}.writer.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def write_manifest(db_path: Path, table: str, value: dict) -> None:
    path = manifest_path(db_path, table)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def reconcile(
    root: Path,
    db_path: Path,
    table_name: str,
    profile: dict,
    include_new: list[str],
    *,
    plan: bool = False,
    only_if_empty: bool = False,
) -> dict:
    import lancedb

    started = time.monotonic()
    with writer_lock(db_path, table_name):
        previous = read_manifest(db_path, table_name)
        check_owner(previous, root)
        additions = sorted(
            set(include_new) | set(previous.get("include_new", []) if previous else [])
        )
        files = selected_files(root, profile["repository"], additions)
        db = lancedb.connect(str(db_path))
        existing = db.open_table(table_name) if table_name in db.table_names() else None
        old = existing.to_arrow().to_pylist() if existing is not None else []
        if existing is not None and not previous and not (only_if_empty and not old):
            raise ValueError("Existing table has no ownership manifest; refusing to replace it")
        if old and any(r.get("worktree_id") != sha256(str(root)) for r in old):
            raise ValueError("Existing rows have a different worktree owner")
        if only_if_empty and old:
            return {
                "table": table_name,
                "root": str(root),
                "mode": "bootstrap",
                "state": "populated",
                "chunks": len(old),
                "new_embeddings": 0,
            }
        model, spec = embedding(profile)
        if previous and previous["embedding"] != spec:
            raise ValueError("Embedding model/configuration changed; use a new table")
        can_reuse = (
            previous
            and existing is not None
            and previous["state"] == "ready"
            and previous.get("table_version") == existing.version
            and previous["chunks"] == len(old)
            and previous["profile_sha256"] == fingerprint(profile)
            and previous["chunker"] == CHUNKER
        )
        records, hashes = make_records(
            root,
            profile,
            files,
            model.embedding_model.tokenizer,
            fingerprint(spec),
            cached_files=previous["files"] if can_reuse else None,
            cached_rows=old if can_reuse else None,
        )
        summary = {
            "table": table_name,
            "root": str(root),
            "files": len(files),
            "chunks": len(records),
            "max_tokens": max(r["tokens"] for r in records),
            "embedding": spec,
            "parsed_files": sum(
                not can_reuse or previous["files"].get(p) != digest for p, digest in hashes.items()
            ),
        }
        if plan:
            return {**summary, "mode": "plan", "seconds": round(time.monotonic() - started, 2)}
        cache = {r["embedding_key"]: r["vector"] for r in old}
        needed = {r["embedding_key"]: r["text"] for r in records if r["embedding_key"] not in cache}
        items = list(needed.items())
        for offset in range(0, len(items), 64):
            batch = items[offset : offset + 64]
            vectors = model.compute_source_embeddings([text for _, text in batch])
            if len(vectors) != len(batch) or any(v is None for v in vectors):
                raise ValueError("Incomplete embedding batch")
            cache.update({key: vector for (key, _), vector in zip(batch, vectors, strict=True)})
            print(f"embedded {min(offset + 64, len(items))}/{len(items)}", file=sys.stderr)
        for row in records:
            row["vector"] = cache[row["embedding_key"]]
        old_by_id = {r["id"]: r for r in old}
        new_ids = {r["id"] for r in records}
        deleted = len(set(old_by_id) - new_ids)
        changed = [
            r
            for r in records
            if r["id"] not in old_by_id
            or any(
                value != old_by_id[r["id"]].get(key) for key, value in r.items() if key != "vector"
            )
        ]
        # Verify a coherent source snapshot before committing derived rows.
        if selected_files(root, profile["repository"], additions) != files or any(
            hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
            for name, digest in hashes.items()
        ):
            raise ValueError("Source changed during indexing; retry catch-up")
        # Persist ownership before publishing, so interrupted first writes can be
        # recovered and an incomplete table is never reported as current.
        manifest = {
            "format": VERSION,
            "root": str(root),
            "worktree_id": sha256(str(root)),
            "table": table_name,
            "git_head": git(root, "rev-parse", "HEAD").decode().strip(),
            "files": hashes,
            "snapshot_id": fingerprint(hashes),
            "profile_sha256": fingerprint(profile),
            "chunker": CHUNKER,
            "include_new": additions,
            "embedding": spec,
            "chunks": len(records),
            "fts": "native_bm25",
            "state": "building",
            "updated_at": time.time(),
        }
        write_manifest(db_path, table_name, manifest)
        if existing is None or (only_if_empty and not old):
            # A mount may adopt an empty table left by a basic connector. Only
            # zero-row tables can have their schema replaced through this hook.
            existing = db.create_table(
                table_name,
                data=records,
                schema=schema(model),
                mode="overwrite" if existing is not None else "create",
            )
        elif changed or deleted:
            # The ownership guard makes this a complete reconciliation of one
            # dedicated worktree table, not an unscoped deletion in a shared corpus.
            (
                existing.merge_insert("id")
                .when_matched_update_all()
                .when_not_matched_insert_all()
                .when_not_matched_by_source_delete()
                .execute(records)
            )
        if (
            changed
            or deleted
            or not previous
            or previous["state"] != "ready"
            or previous.get("fts") != "native_bm25"
        ):
            existing.create_fts_index(
                "text",
                use_tantivy=False,
                replace=True,
                stem=False,
                remove_stop_words=False,
                max_token_length=200,
                writer_heap_size=8 * 1024 * 1024,
            )
        manifest.update(state="ready", table_version=existing.version)
        write_manifest(db_path, table_name, manifest)
        return {
            **summary,
            "mode": "sync",
            "new_embeddings": len(needed),
            "reused_rows": sum(r["embedding_key"] not in needed for r in records),
            "updated_rows": len(changed),
            "deleted_rows": deleted,
            "snapshot_id": manifest["snapshot_id"],
            "fts": manifest["fts"],
            "seconds": round(time.monotonic() - started, 2),
        }


def status(root: Path, db_path: Path, table: str, profile: dict) -> dict:
    import lancedb

    manifest = read_manifest(db_path, table)
    if not manifest:
        return {"table": table, "state": "absent"}
    check_owner(manifest, root)
    db = lancedb.connect(str(db_path))
    if table not in db.table_names():
        return {"table": table, "state": "absent"}
    stored = db.open_table(table)
    healthy = (
        manifest["state"] == "ready"
        and stored.version == manifest.get("table_version")
        and stored.count_rows() == manifest["chunks"]
    )
    current = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in selected_files(root, profile["repository"], manifest["include_new"])
    }
    changed = sorted(
        name
        for name in set(current) | set(manifest["files"])
        if current.get(name) != manifest["files"].get(name)
    )
    policy_changed = (
        manifest["profile_sha256"] != fingerprint(profile) or manifest["chunker"] != CHUNKER
    )
    state = "incomplete" if not healthy else "stale" if changed or policy_changed else "current"
    return {
        "table": table,
        "root": str(root),
        "state": state,
        "changed_paths": changed,
        "profile_changed": policy_changed,
        "files": len(current),
        "chunks": manifest["chunks"],
        "snapshot_id": manifest["snapshot_id"],
    }


def query(db_path: Path, table_name: str, text: str, mode: str, limit: int) -> list[dict]:
    import lancedb

    table = lancedb.connect(str(db_path)).open_table(table_name)
    search = table.search(text, query_type=mode)
    if mode in {"vector", "hybrid"}:
        search = search.metric("cosine")
    return (
        search.limit(limit)
        .select(["path", "anchor", "source_kind", "start_line", "end_line", "doc"])
        .to_list()
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "sync", "bootstrap", "status", "query"):
        command = commands.add_parser(name)
        command.add_argument("--root", type=Path, default=Path.cwd())
        command.add_argument(
            "--db",
            type=Path,
            default=Path(
                os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")
            ).expanduser(),
        )
        command.add_argument("--table", help="Explicit table; default selection is worktree-owned")
        command.add_argument("--profile", type=Path, default=Path(".opencode/search/profiles.json"))
        if name in {"plan", "sync", "bootstrap"}:
            command.add_argument("--include-new", action="append", default=[])
        if name == "query":
            command.add_argument("text")
            command.add_argument("--mode", choices=["vector", "fts", "hybrid"], default="hybrid")
            command.add_argument("--limit", type=int, default=5)
    args = parser.parse_args()
    root = args.root.resolve()
    db_path = args.db.expanduser().resolve()
    profile = json.loads((root / args.profile).read_text())
    args.table = args.table or select_table(root, db_path, profile["repository"]["table"])
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", args.table):
        parser.error("Table name must be an identifier")
    if args.command == "query":
        if not 1 <= args.limit <= 50:
            parser.error("limit must be between 1 and 50")
        freshness = status(root, db_path, args.table, profile)
        if freshness["state"] != "current":
            parser.error(f"Index is {freshness['state']}; run status/sync before querying")
        result = {
            "freshness": freshness,
            "results": query(db_path, args.table, args.text, args.mode, args.limit),
        }
    else:
        if args.command == "status":
            result = status(root, db_path, args.table, profile)
        else:
            result = reconcile(
                root,
                db_path,
                args.table,
                profile,
                args.include_new,
                plan=args.command == "plan",
                only_if_empty=args.command == "bootstrap",
            )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
