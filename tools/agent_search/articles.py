"""Index tracked reference notes without promoting their prose to verified full text."""

from __future__ import annotations

import re
from pathlib import Path

from tools.agent_search.chunks import Unit, markdown_units, split_units

CHUNKER = "reference-note-blocks-v1"
STRING_FIELDS = [
    "paper_id",
    "title",
    "doi_or_url",
    "source_edition",
    "pagination",
    "note_path",
    "extraction_version",
    "verification_status",
]
LIST_FIELDS = ["printed_page", "locator"]
PAGE = re.compile(r"\bpp?\.\s*(?:R)?\d+(?:-\d+)?(?:[–—-](?:R)?\d+)?", re.IGNORECASE)


def declaration(source: str, label: str) -> str:
    """Read a possibly wrapped bold-label declaration, stopping at the next label."""
    match = re.search(
        rf"^\*\*{label}:\*\*\s*(.*?)(?=^\*\*[\w ]+:\*\*|\n\s*\n|\Z)",
        source,
        re.MULTILINE | re.DOTALL,
    )
    if not match:
        raise ValueError(f"Reference note is missing its {label} declaration")
    return " ".join(match.group(1).split())


def paper_metadata(path: str, source: str) -> dict:
    paper_id = Path(path).stem
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", paper_id):
        raise ValueError("Reference note stem must be an identifier usable as paper_id")
    title = re.search(r"^#\s+(.+)$", source, re.MULTILINE)
    if not title:
        raise ValueError("Reference note is missing its title")
    origin = declaration(source, "Source")
    pagination = declaration(source, "Pagination")
    doi = re.search(r"10\.\d{4,9}/[^\s<>`]+", origin)
    url = re.search(r"https?://[^\s<>`]+", origin)
    arxiv = re.search(r"arXiv:([\d.]+(?:v\d+)?)", origin)
    if doi:
        identifier = "https://doi.org/" + doi.group().rstrip(".,;")
    elif url:
        identifier = url.group().rstrip(".,;")
    elif arxiv:
        identifier = "https://arxiv.org/abs/" + arxiv.group(1)
    else:
        raise ValueError("Reference note has no DOI or stable source URL")
    return {
        "paper_id": paper_id,
        "title": title.group(1),
        "doi_or_url": identifier,
        # These are the note's edition declarations, not a guessed journal/preprint
        # mapping. Preserve update/front-matter caveats for the source resolver.
        "source_edition": origin,
        "pagination": pagination,
        "note_path": path,
        "extraction_version": CHUNKER,
    }


def note_units(source: str) -> list[Unit]:
    """Keep paragraphs, table rows and cited multi-line equations together when possible."""
    result = []
    for section in markdown_units(source):
        lines = section.content.splitlines(keepends=True)
        start = 0
        fence = None

        def emit(end, lines=lines, section=section):
            nonlocal start
            content = "".join(lines[start:end]).rstrip("\n")
            if content.strip():
                result.append(
                    Unit(section.anchor, "note_block", section.start_line + start, content)
                )
            start = end

        for i, line in enumerate(lines):
            marker = re.match(r"^\s{0,3}(`{3,}|~{3,})", line)
            if marker:
                emit(i)
                run = marker.group(1)
                if fence is None:
                    fence = (run[0], len(run))
                elif run[0] == fence[0] and len(run) >= fence[1]:
                    fence = None
                start = i + 1
            elif fence:
                # A trailing page locator closes an equation, including its
                # preceding continuation lines. Oversized equations are flagged.
                if PAGE.search(line):
                    emit(i + 1)
            elif not line.strip():
                emit(i)
                start = i + 1
            elif line.lstrip().startswith("|"):
                emit(i)
                emit(i + 1)
            elif re.match(r"^[-*]\s|^\d+\.\s", line):
                emit(i)
        emit(len(lines))
    return result


def split_note(path: str, source: str, tokenizer, settings: dict):
    metadata = paper_metadata(path, source)
    language, chunks = split_units(
        path,
        source,
        tokenizer,
        settings,
        units=note_units(source),
        prefer_whole=True,
    )
    return language, chunks, metadata


def citation_metadata(content: str, kind: str) -> dict:
    """Return literal locator clauses present in this excerpt; never invent page offsets."""
    locators = []
    pages = []
    # Restrict clauses to one Markdown cell/line and stop at the next citation.
    # These are discovery locators, not an assertion that every clause refers to
    # the owning paper: notes can compare several sources in one paragraph.
    for line in content.splitlines():
        matches = list(PAGE.finditer(line))
        for i, match in enumerate(matches):
            pages.append(re.sub(r"^pp?\.\s*", "", match.group(), flags=re.IGNORECASE))
            end = matches[i + 1].start() if i + 1 < len(matches) else len(line)
            clause = line[match.start() : end].split("|", 1)[0].strip().rstrip(".;")
            locators.append(clause)
    return {
        "printed_page": list(dict.fromkeys(pages)),
        "locator": list(dict.fromkeys(locators)),
        "extraction_page": None,
        "verification_status": "note_fragment" if kind == "note_fragment" else "tracked_note",
    }


def filtered_hybrid(
    table, text: str, predicate: str, columns: list[str], limit: int, *, vector
) -> list[dict]:
    """Prefilter both candidate sets before reciprocal-rank fusion on one pinned table."""
    from lancedb.rerankers import RRFReranker

    # LanceDB 0.21.2's hybrid builder passes _postfilter as the child builder's
    # prefilter argument, reversing its meaning. Public child queries avoid that
    # bug: post-filtering a global top-k loses relevant hits in a smaller paper.
    vector = (
        table.search(vector, query_type="vector")
        .metric("cosine")
        .where(predicate, prefilter=True)
        .limit(limit)
        .select(columns)
        .with_row_id(True)
        .to_arrow()
    )
    lexical = (
        table.search(text, query_type="fts")
        .where(predicate, prefilter=True)
        .limit(limit)
        .select(columns)
        .with_row_id(True)
        .to_arrow()
    )
    fused = RRFReranker().rerank_hybrid(text, vector, lexical)
    return fused.slice(0, limit).drop(["_rowid"]).to_pylist()
