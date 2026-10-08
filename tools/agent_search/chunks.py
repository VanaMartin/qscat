"""Source-anchored units and tokenizer-bounded embedding payloads."""

from __future__ import annotations

import ast
import hashlib
import re
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class Unit:
    anchor: str
    kind: str
    start_line: int
    content: str


@dataclass(frozen=True)
class Chunk:
    anchor: str
    kind: str
    declaration: int
    start_line: int
    end_line: int
    content: str
    text: str
    content_sha256: str
    payload_sha256: str
    duplicate: int
    tokens: int


def sha256(text: str) -> str:
    """Fingerprint exact text rather than source coordinates."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def python_units(source: str) -> list[Unit]:
    """Keep methods independent of surrounding declarations and line shifts."""
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    units: list[Unit] = []

    def add(anchor: str, kind: str, start: int, end: int) -> None:
        content = "".join(lines[start - 1 : end]).rstrip("\n")
        if content.strip():
            units.append(Unit(anchor, kind, start, content))

    def scope(nodes: list[ast.stmt], parent: str, start: int, end: int) -> None:
        cursor = start
        for node in nodes:
            if not isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                continue
            first = min([node.lineno, *(d.lineno for d in node.decorator_list)])
            last = node.end_lineno
            assert last is not None
            add(parent or "@module", "class_context" if parent else "module", cursor, first - 1)
            name = f"{parent}.{node.name}" if parent else node.name
            if isinstance(node, ast.ClassDef):
                # Include the class signature/decorators and fields in its context;
                # methods get their own units instead of duplicating the whole class.
                scope(node.body, name, first, last)
            else:
                add(name, "method" if parent else "function", first, last)
            cursor = last + 1
        add(parent or "@module", "class_context" if parent else "module", cursor, end)

    scope(tree.body, "", 1, len(lines))
    return units


def markdown_units(source: str) -> list[Unit]:
    """Use heading ancestry, ignoring headings inside fenced examples."""
    lines = source.splitlines(keepends=True)
    units: list[Unit] = []
    parents: list[tuple[int, str]] = []
    anchor = "@document"
    start = 1
    fence: tuple[str, int] | None = None
    for number, line in enumerate(lines, 1):
        marker = re.match(r"^\s{0,3}(`{3,}|~{3,})", line)
        if marker:
            run = marker.group(1)
            if fence is None:
                fence = (run[0], len(run))
            elif run[0] == fence[0] and len(run) >= fence[1]:
                fence = None
            continue
        if fence is not None:
            continue
        heading = re.match(r"^(#{1,6})\s+(.+?)\s*#*\s*$", line)
        if not heading:
            continue
        content = "".join(lines[start - 1 : number - 1]).rstrip("\n")
        if content.strip():
            units.append(Unit(anchor, "section", start, content))
        level = len(heading.group(1))
        while parents and parents[-1][0] >= level:
            parents.pop()
        parents.append((level, heading.group(2)))
        anchor = " > ".join(title for _, title in parents)
        start = number
    content = "".join(lines[start - 1 :]).rstrip("\n")
    if content.strip():
        units.append(Unit(anchor, "section", start, content))
    return units


def rust_units(source: str) -> list[Unit]:
    """Parse Rust items and impl methods with Tree-sitter, preserving attributes."""
    import tree_sitter_rust
    from tree_sitter import Language, Parser

    raw = source.encode("utf-8")
    tree = Parser(Language(tree_sitter_rust.language())).parse(raw)
    if tree.root_node.has_error:
        raise ValueError("Rust parse failed")
    units = []
    containers = {"impl_item", "trait_item", "mod_item"}
    items = containers | {
        "function_item",
        "function_signature_item",
        "struct_item",
        "enum_item",
        "type_item",
        "const_item",
        "static_item",
        "macro_definition",
    }

    def add(anchor: str, kind: str, start: int, end: int) -> None:
        content = raw[start:end].decode("utf-8").rstrip("\n")
        if content.strip():
            units.append(Unit(anchor, kind, raw[:start].count(b"\n") + 1, content))

    def scope(node, parent: str, start: int, end: int) -> None:
        cursor = start
        pending = None
        for child in node.named_children:
            if child.type in {"attribute_item", "line_comment", "block_comment"}:
                if pending is None:
                    pending = child.start_byte
                continue
            if child.type not in items:
                pending = None
                continue
            first = child.start_byte if pending is None else pending
            pending = None
            add(parent or "@module", "item_context" if parent else "module", cursor, first)
            field = child.child_by_field_name("name") or child.child_by_field_name("type")
            label = raw[field.start_byte : field.end_byte].decode("utf-8") if field else child.type
            name = f"{parent}::{label}" if parent else label
            body = child.child_by_field_name("body")
            if child.type in containers and body:
                scope(body, name, first, child.end_byte)
            else:
                add(name, child.type.removesuffix("_item"), first, child.end_byte)
            cursor = child.end_byte
        add(parent or "@module", "item_context" if parent else "module", cursor, end)

    scope(tree.root_node, "", 0, len(raw))
    return units


def source_units(path: str, source: str) -> tuple[str, list[Unit]]:
    """Select structural parsers for Python/Markdown and a file unit otherwise."""
    suffix = path.rsplit(".", 1)[-1].lower() if "." in path else ""
    if suffix == "py":
        return "python", python_units(source)
    if suffix in {"md", "mdx"}:
        return "markdown", markdown_units(source)
    if suffix == "rs":
        return "rust", rust_units(source)
    language = {
        "rs": "rust",
        "toml": "toml",
        "yaml": "yaml",
        "yml": "yaml",
        "json": "json",
        "jsonc": "json",
        "sh": "shell",
        "rst": "rst",
    }.get(suffix, "text")
    return language, [Unit("@file", "file", 1, source)] if source.strip() else []


def split_units(
    path: str,
    source: str,
    tokenizer,
    settings: dict,
    *,
    units: list[Unit] | None = None,
    prefer_whole: bool = False,
) -> tuple[str, list[Chunk]]:
    """Count the complete payload, using line boundaries when they fit the budget."""
    if units is None:
        language, units = source_units(path, source)
    else:
        language = "markdown"
    target = settings["chunk_target_tokens"]
    maximum = settings["chunk_max_tokens"]
    total = settings["total_input_max_tokens"]
    overlap = settings["overlap_tokens"]
    if not 0 <= overlap < target <= maximum < total:
        raise ValueError("Invalid tokenizer budgets")
    declarations: Counter = Counter()
    chunks: list[Chunk] = []
    for unit in units:
        declaration = declarations[(unit.anchor, unit.kind)]
        declarations[(unit.anchor, unit.kind)] += 1
        prefix = f"{language} {unit.kind}: {unit.anchor}\n{path}\n\n"
        encoded = tokenizer(
            unit.content, add_special_tokens=False, return_offsets_mapping=True, verbose=False
        )
        offsets = encoded["offset_mapping"]
        if not offsets:
            continue
        whole = (
            prefer_whole
            and len(offsets) <= maximum
            and len(
                tokenizer(prefix + unit.content, add_special_tokens=True, verbose=False)[
                    "input_ids"
                ]
            )
            <= total
        )
        starts = [pair[0] for pair in offsets]
        duplicates: Counter = Counter()
        first = 0
        while first < len(offsets):
            last = len(offsets) if whole else min(first + target, len(offsets))
            lo = 0 if first == 0 else starts[first]
            hi = len(unit.content) if last == len(offsets) else starts[last]
            # Preserve whole source lines when doing so does not make tiny chunks.
            if last < len(offsets):
                boundary = unit.content.rfind("\n", lo, hi)
                end_token = bisect_left(starts, boundary + 1)
                if end_token - first >= target // 2:
                    last, hi = end_token, boundary + 1
            content = unit.content[lo:hi]
            count = len(
                tokenizer(prefix + content, add_special_tokens=True, verbose=False)["input_ids"]
            )
            while count > total and last > first + 1:
                last -= 1
                hi = starts[last]
                content = unit.content[lo:hi]
                count = len(
                    tokenizer(prefix + content, add_special_tokens=True, verbose=False)["input_ids"]
                )
            if count > total:
                raise ValueError(f"Anchor leaves no embedding budget: {path}: {unit.anchor}")
            content_hash = sha256(content)
            duplicate = duplicates[content_hash]
            duplicates[content_hash] += 1
            chunks.append(
                Chunk(
                    unit.anchor,
                    "note_fragment" if prefer_whole and not whole else unit.kind,
                    declaration,
                    unit.start_line + unit.content[:lo].count("\n"),
                    unit.start_line + unit.content[:hi].count("\n") - int(content.endswith("\n")),
                    content,
                    prefix + content,
                    content_hash,
                    sha256(prefix + content),
                    duplicate,
                    count,
                )
            )
            if last == len(offsets):
                break
            first = max(first + 1, last - overlap)
    return language, chunks
