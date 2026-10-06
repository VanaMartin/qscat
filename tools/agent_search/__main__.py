"""Publish and search a shared corpus sourced exclusively from fetched upstream main."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import uuid
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

from tools.agent_search.chunks import sha256, split_units
from tools.agent_search.git_source import PROFILE, MainSource, Snapshot, read_blobs

VERSION = "qmodeling-main-search-v2"
CHUNKER = "python-ast-rust-treesitter-markdown-v2"
CHECK_INTERVAL = 300


def fingerprint(value) -> str:
    return sha256(json.dumps(value, sort_keys=True, ensure_ascii=False))


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
    snapshot: Snapshot,
    tokenizer,
    model_key: str,
    *,
    previous: dict | None = None,
    cached_rows: list[dict] | None = None,
):
    records, hashes = [], {}
    by_path: dict[str, list[dict]] = {}
    for row in cached_rows or []:
        by_path.setdefault(row["path"], []).append(row)
    selected = snapshot.selected()
    unchanged = {
        path for path, oid in selected.items() if previous and previous["blobs"].get(path) == oid
    }
    raw_blobs = read_blobs(
        snapshot.source.root,
        list(dict.fromkeys(oid for path, oid in selected.items() if path not in unchanged)),
    )
    for path, oid in selected.items():
        if path in unchanged:
            hashes[path] = previous["files"][path]
            records.extend(
                {k: v for k, v in row.items() if k != "vector"} for row in by_path.get(path, [])
            )
            continue
        raw = raw_blobs[oid]
        digest = hashes[path] = hashlib.sha256(raw).hexdigest()
        try:
            language, chunks = split_units(
                path, raw.decode("utf-8"), tokenizer, snapshot.profile["embedding"]
            )
        except (SyntaxError, UnicodeError, ValueError) as error:
            raise ValueError(f"Cannot index {path}: {error}") from error
        for chunk in chunks:
            row_id = fingerprint(
                [
                    snapshot.source.repo_id,
                    path,
                    language,
                    chunk.kind,
                    chunk.anchor,
                    chunk.declaration,
                    chunk.content_sha256,
                    chunk.duplicate,
                ]
            )
            records.append(
                {
                    "id": row_id,
                    "text": chunk.text,
                    "path": path,
                    "anchor": chunk.anchor,
                    "kind": chunk.kind,
                    "language": language,
                    "source_kind": source_kind(path),
                    "repo_id": snapshot.source.repo_id,
                    "source_blob": oid,
                    "start_line": chunk.start_line,
                    "end_line": chunk.end_line,
                    "declaration": chunk.declaration,
                    "source_sha256": digest,
                    "content_sha256": chunk.content_sha256,
                    "payload_sha256": chunk.payload_sha256,
                    "embedding_key": fingerprint([model_key, chunk.payload_sha256]),
                    "tokens": chunk.tokens,
                    "doc": (
                        f"source={path}\nanchor={chunk.anchor}\nkind={chunk.kind}\n"
                        f"declaration={chunk.declaration}\n"
                        f"lines={chunk.start_line}-{chunk.end_line}\n"
                        f"repo={snapshot.source.repo_id}\nsource_blob={oid}\n"
                        f"source_sha256={digest}\nchunk_id={row_id}\n\n{chunk.content}"
                    ),
                }
            )
    if not records:
        raise ValueError("No searchable chunks selected; refusing to publish an empty corpus")
    if len({row["id"] for row in records}) != len(records):
        raise ValueError("Chunk IDs are not unique")
    return records, hashes, len(selected) - len(unchanged)


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
        "source_blob",
        "source_sha256",
        "content_sha256",
        "payload_sha256",
        "embedding_key",
    ]
    fields = {name: (str, ...) for name in strings}
    fields.update(
        {name: (int, ...) for name in ["start_line", "end_line", "declaration", "tokens"]}
    )
    fields["text"] = (str, model.SourceField())
    fields["vector"] = (Vector(model.ndims()), model.VectorField())
    return create_model("MainCodeChunk", __base__=LanceModel, **fields).to_arrow_schema()


def read_manifest(db_path: Path, table: str) -> dict | None:
    path = db_path / f"{table}.manifest.json"
    return json.loads(path.read_text()) if path.exists() else None


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def check_owner(manifest: dict | None, source: MainSource) -> None:
    if manifest and (manifest.get("format") != VERSION or manifest.get("owner") != source.owner):
        raise ValueError("Corpus belongs to a different upstream, profile, or index format")


@contextmanager
def writer_lock(db_path: Path, table: str):
    db_path.mkdir(parents=True, exist_ok=True)
    with (db_path / f"{table}.writer.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def open_published(db_path: Path, manifest: dict):
    """Pin one immutable generation; never reread the publication pointer mid-query."""
    import lancedb

    db = lancedb.connect(str(db_path))
    name = manifest["published_table"]
    if name not in db.table_names():
        raise ValueError("Published generation is absent")
    table = db.open_table(name)
    if table.version != manifest["table_version"]:
        raise ValueError("Published generation was externally modified")
    table.checkout(manifest["table_version"])
    if table.count_rows() != manifest["chunks"]:
        raise ValueError("Published generation has an unexpected row count")
    return table


def describe(
    source: MainSource,
    db_path: Path,
    table: str,
    *,
    manifest: dict | None,
    interval: float = CHECK_INTERVAL,
) -> dict:
    check_owner(manifest, source)
    refresh_path = db_path / f"{table}.refresh.json"
    refresh = (
        json.loads(refresh_path.read_text()) if refresh_path.exists() else {"state": "pending"}
    )
    if refresh.get("owner") and refresh["owner"] != source.owner:
        raise ValueError("Refresh record belongs to a different upstream")
    integrity = "absent"
    error = None
    if manifest:
        try:
            open_published(db_path, manifest)
            integrity = "complete"
        except ValueError as failure:
            integrity, error = "incomplete", str(failure)
    checked = refresh.get("checked_at")
    recent = checked is not None and time.time() - checked <= interval
    currency = "unknown"
    if refresh.get("fetch_ok") and recent and manifest:
        currency = "current" if manifest["source_commit"] == refresh["target_commit"] else "behind"
    state = (
        integrity if integrity != "complete" else "current" if currency == "current" else "stale"
    )
    return {
        "table": table,
        "repo_id": source.repo_id,
        "source_ref": f"{source.remote}/main",
        "source_commit": manifest["source_commit"] if manifest else None,
        "state": state,
        "integrity": integrity,
        "upstream": {
            "state": currency,
            "checked_at": checked,
            "target_commit": refresh.get("target_commit"),
        },
        "refresh": refresh,
        "error": error,
        "files": len(manifest["files"]) if manifest else 0,
        "chunks": manifest["chunks"] if manifest else 0,
        "snapshot_id": manifest["snapshot_id"] if manifest else None,
        "published_table": manifest["published_table"] if manifest else None,
    }


def status(source: MainSource, db_path: Path, table: str, *, interval=CHECK_INTERVAL) -> dict:
    return describe(
        source, db_path, table, manifest=read_manifest(db_path, table), interval=interval
    )


def legacy_vectors(source: MainSource, db_path: Path, db) -> dict:
    """Import only a vector cache; legacy locations and ownership are never published."""
    cache = {}
    for path in sorted(db_path.glob("*.manifest.json")):
        old = json.loads(path.read_text())
        if old.get("format") != "qmodeling-search-v1":
            continue
        root = Path(old["root"])
        try:
            owner = MainSource.discover(root)
            stored = db.open_table(old["table"])
        except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError):
            continue
        if owner.repo_id != source.repo_id or old.get("state") != "ready":
            continue
        if stored.version != old.get("table_version") or stored.count_rows() != old["chunks"]:
            continue
        for row in stored.to_arrow().to_pylist():
            # A legacy cache contributes only self-consistent payload/model pairs.
            if sha256(row["text"]) == row["payload_sha256"] and row["embedding_key"] == fingerprint(
                [fingerprint(old["embedding"]), row["payload_sha256"]]
            ):
                cache[row["embedding_key"]] = row["vector"]
    return cache


def reconcile(snapshot: Snapshot, db_path: Path, table: str, *, plan: bool = False) -> dict:
    """Build under the caller's writer lock and atomically publish a complete generation."""
    import lancedb

    started = time.monotonic()
    previous = read_manifest(db_path, table)
    check_owner(previous, snapshot.source)
    old = []
    healthy = False
    if previous:
        try:
            old = open_published(db_path, previous).to_arrow().to_pylist()
            healthy = True
        except ValueError:
            pass
    selected = snapshot.selected()
    policy_key = fingerprint(snapshot.profile)
    if healthy and previous["profile_sha256"] == policy_key and previous["chunker"] == CHUNKER:
        if previous["blobs"] == selected:
            if not plan and previous["source_commit"] != snapshot.commit:
                updated = {
                    **previous,
                    "source_commit": snapshot.commit,
                    "source_tree": snapshot.tree,
                    "updated_at": time.time(),
                }
                write_json(db_path / f"{table}.manifest.json", updated)
            return {
                "mode": "plan" if plan else "sync",
                "table": table,
                "source_commit": snapshot.commit,
                "files": len(selected),
                "chunks": previous["chunks"],
                "parsed_files": 0,
                "new_embeddings": 0,
                "updated_rows": 0,
                "deleted_rows": 0,
                "reused_rows": len(old),
                "snapshot_id": previous["snapshot_id"],
                "seconds": time.monotonic() - started,
            }
    model, spec = embedding(snapshot.profile)
    if previous and previous["embedding"] != spec:
        raise ValueError("Embedding model/configuration changed; use a new logical table")
    can_reuse = (
        healthy and previous["profile_sha256"] == policy_key and previous["chunker"] == CHUNKER
    )
    records, hashes, parsed = make_records(
        snapshot,
        model.embedding_model.tokenizer,
        fingerprint(spec),
        previous=previous if can_reuse else None,
        cached_rows=old if can_reuse else None,
    )
    summary = {
        "table": table,
        "source_commit": snapshot.commit,
        "files": len(selected),
        "chunks": len(records),
        "max_tokens": max(r["tokens"] for r in records),
        "parsed_files": parsed,
        "embedding": spec,
    }
    if plan:
        return {**summary, "mode": "plan"}
    db = lancedb.connect(str(db_path))
    if table in db.table_names():
        raise ValueError("Logical name collides with an existing table; use a new logical table")
    cache = {row["embedding_key"]: row["vector"] for row in old}
    if not previous:
        cache.update(legacy_vectors(snapshot.source, db_path, db))
    needed = {
        row["embedding_key"]: row["text"] for row in records if row["embedding_key"] not in cache
    }
    items = list(needed.items())
    for offset in range(0, len(items), 64):
        batch = items[offset : offset + 64]
        vectors = model.compute_source_embeddings([text for _, text in batch])
        if len(vectors) != len(batch) or any(vector is None for vector in vectors):
            raise ValueError("Incomplete embedding batch")
        cache.update({key: vector for (key, _), vector in zip(batch, vectors, strict=True)})
        print(f"embedded {min(offset + 64, len(items))}/{len(items)}", file=sys.stderr)
    for row in records:
        row["vector"] = cache[row["embedding_key"]]
    old_by_id = {row["id"]: row for row in old}
    new_ids = {row["id"] for row in records}
    changed = sum(
        row["id"] not in old_by_id
        or any(
            value != old_by_id[row["id"]].get(key) for key, value in row.items() if key != "vector"
        )
        for row in records
    )
    generation = f"{table}_g_{uuid.uuid4().hex}"
    candidate = db.create_table(generation, data=records, schema=schema(model))
    candidate.create_fts_index(
        "text",
        use_tantivy=False,
        replace=True,
        stem=False,
        remove_stop_words=False,
        max_token_length=200,
        writer_heap_size=8 * 1024 * 1024,
    )
    manifest = {
        "format": VERSION,
        "owner": snapshot.source.owner,
        "table": table,
        "published_table": generation,
        "table_version": candidate.version,
        "source_commit": snapshot.commit,
        "source_tree": snapshot.tree,
        "files": hashes,
        "blobs": selected,
        "snapshot_id": fingerprint([hashes, policy_key, CHUNKER]),
        "profile_sha256": policy_key,
        "chunker": CHUNKER,
        "embedding": spec,
        "chunks": len(records),
        "fts": "native_bm25",
        "updated_at": time.time(),
    }
    open_published(db_path, manifest)
    # The complete manifest is the only publication pointer. Existing generations
    # remain readable by requests in other processes that captured the old pointer.
    write_json(db_path / f"{table}.manifest.json", manifest)
    return {
        **summary,
        "mode": "sync",
        "new_embeddings": len(needed),
        "reused_rows": sum(row["embedding_key"] not in needed for row in records),
        "updated_rows": changed,
        "deleted_rows": len(set(old_by_id) - new_ids),
        "snapshot_id": manifest["snapshot_id"],
        "seconds": round(time.monotonic() - started, 2),
    }


def ensure_main(
    source: MainSource,
    db_path: Path,
    table: str,
    *,
    force: bool = False,
    interval: float = CHECK_INTERVAL,
    plan: bool = False,
) -> dict:
    with writer_lock(db_path, table):
        previous = read_manifest(db_path, table)
        check_owner(previous, source)
        refresh_path = db_path / f"{table}.refresh.json"
        current = describe(source, db_path, table, manifest=previous, interval=interval)
        attempted = current["refresh"].get("attempted_at", 0)
        if (
            not force
            and not plan
            and time.time() - attempted < interval
            and (
                current["integrity"] == "complete"
                and current["refresh"]["state"] in {"complete", "failed"}
            )
        ):
            return {"mode": "reuse", **current}
        refresh = {
            "owner": source.owner,
            "state": "checking",
            "attempted_at": time.time(),
            "fetch_ok": False,
        }
        if not plan:
            write_json(refresh_path, refresh)
        try:
            snapshot = source.fetch()
            refresh.update(
                state="indexing",
                fetch_ok=True,
                checked_at=time.time(),
                target_commit=snapshot.commit,
            )
            if not plan:
                write_json(refresh_path, refresh)
            result = reconcile(snapshot, db_path, table, plan=plan)
            refresh.update(state="complete", completed_at=time.time())
            if not plan:
                write_json(refresh_path, refresh)
            return result
        except Exception as error:
            refresh.update(state="failed", error=str(error))
            if not plan:
                write_json(refresh_path, refresh)
            raise


def query(
    source: MainSource,
    db_path: Path,
    table: str,
    text: str,
    mode: str,
    limit: int,
    *,
    interval: float = CHECK_INTERVAL,
    manifest: dict | None = None,
) -> dict:
    if not 1 <= limit <= 50 or mode not in {"vector", "fts", "hybrid"}:
        raise ValueError("Use a supported query mode and a limit between 1 and 50")
    manifest = manifest or read_manifest(db_path, table)
    freshness = describe(source, db_path, table, manifest=manifest, interval=interval)
    if freshness["integrity"] != "complete":
        error = freshness["refresh"].get("error") or freshness["error"] or "inspect refresh status"
        raise RuntimeError(f"Main index is {freshness['integrity']}: {error}")
    stored = open_published(db_path, manifest)
    search = stored.search(text, query_type=mode)
    if mode in {"vector", "hybrid"}:
        search = search.metric("cosine")
    results = (
        search.limit(limit)
        .select(
            [
                "id",
                "path",
                "anchor",
                "source_kind",
                "start_line",
                "end_line",
                "declaration",
                "source_sha256",
                "source_blob",
                "doc",
            ]
        )
        .to_list()
    )
    return {"freshness": freshness, "results": results}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("plan", "sync", "bootstrap", "status", "query"):
        command = commands.add_parser(name)
        command.add_argument("--root", type=Path, default=Path.cwd())
        command.add_argument("--remote", default="origin")
        command.add_argument(
            "--db",
            type=Path,
            default=Path(os.environ.get("LANCEDB_URI", "~/.local/share/opencode/lancedb")),
        )
        command.add_argument("--table", help="Logical corpus name; defaults to upstream identity")
        command.add_argument(
            "--profile", default=PROFILE, help="Committed repository-relative policy"
        )
        if name == "query":
            command.add_argument("text")
            command.add_argument("--mode", choices=["vector", "fts", "hybrid"], default="hybrid")
            command.add_argument("--limit", type=int, default=5)
    args = parser.parse_args()
    source = MainSource.discover(args.root, args.remote, args.profile)
    db_path = args.db.expanduser().resolve()
    table = args.table or source.table
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", table):
        parser.error("Table name must be an identifier")
    if args.command == "status":
        result = status(source, db_path, table)
    elif args.command == "query":
        result = query(source, db_path, table, args.text, args.mode, args.limit)
    else:
        result = ensure_main(source, db_path, table, force=True, plan=args.command == "plan")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
