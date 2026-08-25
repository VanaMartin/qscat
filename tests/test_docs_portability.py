"""The docs/physics notes must render in two places.

They are the working notes of the repository -- read as files in a clone --
and pages on the published site. GitHub's Markdown renderer handles
`$...$` and `$$...$$`, but NOT conf.py-defined MathJax macros, and it shows
MyST directives such as `{dropdown}` as visible junk.

So the notes are restricted to constructs both renderers understand. This
test enforces that; see docs/superpowers/specs/2026-08-19-docs-latex-and-
theory-ia-design.md, decisions D1 and D2.

Site-first pages -- the molecule guides and the section index pages, which
nobody reads raw -- are exempt, by an explicit allowlist rather than a path
rule, so adding an eighth is a visible edit here and not a silent opt-out.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PHYSICS_DIR = REPO_ROOT / "docs" / "physics"
CONF_PY = REPO_ROOT / "docs" / "conf.py"

# Exempt from the note rules: read on the site, never raw.
SITE_FIRST_PAGES = frozenset(
    {
        "README.md",
        "discretisation.md",
        "solvers.md",
        "engine.md",
        "dissociation.md",
        "time-dependent.md",
        "resonances.md",
        "open-directions.md",
    }
)

# ``` fenced blocks: their contents are not prose and are skipped.
_FENCE_RE = re.compile(r"^```.*?^```", re.DOTALL | re.MULTILINE)
_DISPLAY_MATH_RE = re.compile(r"\$\$.*?\$\$", re.DOTALL)

# A MyST directive opener: ```{name} or :::{name}, any name. Scanned
# unfenced-stripped, because the directive IS the fence. The scan is
# deliberately name-agnostic -- an allowlist of the six sphinx-design
# directives that happened to exist when this was written let {figure},
# {doc} and {toctree} through, and those render as literal source on
# github.com exactly like {dropdown} does.
_DIRECTIVE_RE = re.compile(r"^(?:```+|:::+)\{([a-z][a-z0-9_-]*)\}", re.MULTILINE)

# A MyST inline role: {name}`content`. The lookbehind is what keeps maths
# out of it -- `V^{-*}_{dk}` inside a code span ends in `}` followed by a
# backtick and would otherwise read as a role named `dk`. A real role
# opener follows whitespace or punctuation, never an identifier character.
_ROLE_RE = re.compile(r"(?<![\w^_\\])\{([a-z][a-z0-9_-]*)\}`")

# A markdown ATX heading.
_HEADING_RE = re.compile(r"^#{1,6} +(.*)$", re.MULTILINE)


def _strip_fences(text: str) -> str:
    return _FENCE_RE.sub("", text)


def find_myst_directives(text: str) -> list[str]:
    """Return every MyST/sphinx-design directive opener found in `text`."""
    return [m.group(0) for m in _DIRECTIVE_RE.finditer(text)]


def find_myst_roles(text: str) -> list[str]:
    """Return every MyST inline role (``{doc}`x```, ``{ref}`x```) in `text`."""
    return [f"{{{m.group(1)}}}" for m in _ROLE_RE.finditer(_strip_fences(text))]


def find_math_in_headings(text: str) -> list[str]:
    """Return every heading containing a `$`.

    A heading also renders in the sidebar and in `toctree` entries, where
    MathJax does not run, so `$\\sigma$` in a title shows as literal source
    there. Headings stay plain unicode.
    """
    return [m.group(1) for m in _HEADING_RE.finditer(_strip_fences(text)) if "$" in m.group(1)]


def find_bare_math_environments(text: str) -> list[str]:
    """Return `\\begin{align|gather|equation}` used OUTSIDE a `$$` block.

    GitHub renders these only when wrapped in `$$`; bare, they show as
    literal LaTeX source.
    """
    outside = _DISPLAY_MATH_RE.sub("", _strip_fences(text))
    pattern = re.compile(r"\\begin\{(align\*?|gather\*?|equation\*?)\}")
    return [m.group(0) for m in pattern.finditer(outside)]


def find_macro_definitions(text: str) -> list[str]:
    """Return every LaTeX macro definition; none may appear in a note."""
    pattern = re.compile(r"\\(newcommand|renewcommand|def)\b")
    return [m.group(0) for m in pattern.finditer(text)]


def _notes() -> list[Path]:
    return sorted(p for p in PHYSICS_DIR.glob("*.md") if p.name not in SITE_FIRST_PAGES)


# --- detector unit tests -------------------------------------------------


def test_find_myst_directives_flags_a_dropdown():
    text = "Some prose.\n\n```{dropdown} Parameters\ncontent\n```\n"
    assert find_myst_directives(text) == ["```{dropdown}"]


def test_find_myst_directives_flags_colon_fence_form():
    text = ":::{grid} 2\n:::\n"
    assert find_myst_directives(text) == [":::{grid}"]


def test_find_myst_directives_flags_a_figure():
    """The name-agnostic scan is the point: {figure} is not sphinx-design."""
    text = "```{figure} figures/x.png\n:width: 90%\n\nA caption.\n```\n"
    assert find_myst_directives(text) == ["```{figure}"]


def test_find_myst_directives_ignores_a_language_code_fence():
    text = "```python\nx = 1\n```\n"
    assert find_myst_directives(text) == []


def test_find_myst_directives_passes_a_plain_note():
    text = "Prose with $\\sigma$ and a table.\n\n| a | b |\n|---|---|\n"
    assert find_myst_directives(text) == []


def test_find_myst_roles_flags_doc_and_ref():
    text = "See {doc}`h2plus-dr` and {ref}`some-label` for more.\n"
    assert find_myst_roles(text) == ["{doc}", "{ref}"]


def test_find_myst_roles_ignores_a_subscripted_identifier():
    """`V^{-*}_{dk}` ends in `}` + a backtick; it is not a role."""
    text = "the conjugated `V^{-*}_{dk}` becomes `V^+_{dk}` here\n"
    assert find_myst_roles(text) == []


def test_find_myst_roles_ignores_display_maths():
    text = "$$\n\\begin{aligned} a &= b \\end{aligned}\n$$\n\n`code`\n"
    assert find_myst_roles(text) == []


def test_find_math_in_headings_flags_a_dollar_title():
    text = "## 4. PRA 77's $V_d$ is not qscat's `Vd`\n\nbody\n"
    assert find_math_in_headings(text) == ["4. PRA 77's $V_d$ is not qscat's `Vd`"]


def test_find_math_in_headings_passes_a_unicode_title():
    text = "## 8.2 The \u03c6\u207b gate\n\nbody with $\\phi^{-}$ maths\n"
    assert find_math_in_headings(text) == []


def test_find_math_in_headings_ignores_a_comment_in_a_fence():
    text = "```bash\n# a $VAR shell comment\n```\n"
    assert find_math_in_headings(text) == []


def test_find_bare_math_environments_flags_a_bare_align():
    text = "\\begin{align} a &= b \\end{align}\n"
    assert find_bare_math_environments(text) == ["\\begin{align}"]


def test_find_bare_math_environments_allows_aligned_inside_display_math():
    text = "$$\n\\begin{aligned} a &= b \\end{aligned}\n$$\n"
    assert find_bare_math_environments(text) == []


def test_find_bare_math_environments_allows_align_inside_display_math():
    text = "$$\n\\begin{align} a &= b \\end{align}\n$$\n"
    assert find_bare_math_environments(text) == []


def test_find_macro_definitions_flags_newcommand():
    text = "\\newcommand{\\Psip}{\\Psi^{(+)}}\n"
    assert find_macro_definitions(text) == ["\\newcommand"]


def test_find_macro_definitions_passes_plain_latex():
    text = "$\\Psi^{(+)}$ and $\\sigma_\\mathrm{DA}$\n"
    assert find_macro_definitions(text) == []


# --- tree scan -----------------------------------------------------------


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_myst_directives(note: Path):
    found = find_myst_directives(note.read_text())
    assert not found, (
        f"{note.name} uses MyST directives {found}, which github.com renders as "
        f"visible junk. Use a plain markdown image with the caption as ordinary "
        f"prose. Directives belong in docs/molecules/ or a section index page "
        f"(see SITE_FIRST_PAGES)."
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_myst_roles(note: Path):
    found = sorted(set(find_myst_roles(note.read_text())))
    assert not found, (
        f"{note.name} uses MyST roles {found}; github.com shows the role name "
        f"and the backticks literally. Link with a relative markdown link -- "
        f"[the note](other-note.md) -- which resolves in both renderers."
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_math_in_headings(note: Path):
    found = find_math_in_headings(note.read_text())
    assert not found, (
        f"{note.name} has maths in the headings {found}. A heading also renders "
        f"in the sidebar and in toctree entries, where MathJax does not run. "
        f"Write the symbol as plain unicode there and keep the maths in the body."
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_bare_math_environments(note: Path):
    found = find_bare_math_environments(note.read_text())
    assert not found, (
        f"{note.name} uses {found} outside a $$ block; github.com will show it "
        f"as literal LaTeX. Wrap it: $$ \\begin{{aligned}} ... \\end{{aligned}} $$"
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_defines_no_macros(note: Path):
    found = find_macro_definitions(note.read_text())
    assert not found, (
        f"{note.name} defines a LaTeX macro {found}; macros do not render on github.com."
    )


# Every key sphinx.ext.mathjax reads. Guarding only one of them is a hole:
# Sphinx 9.1 defaults to MathJax 4 (MATHJAX_URL is .../mathjax@4/...), so
# `mathjax4_config` is the LIVE route and `mathjax3_config` is inert against
# that default -- a macro added via the live key would sail past a check that
# watched the dead one.
_MATHJAX_CONFIG_KEYS = (
    "mathjax_config",
    "mathjax2_config",
    "mathjax3_config",
    "mathjax4_config",
    "mathjax_config_path",
)


def test_conf_py_defines_no_mathjax_macros():
    """The other half of 'no custom macros' (spec D1).

    Checks for an actual assignment, not a bare substring match -- conf.py's
    own comment explaining that no such config is defined legitimately
    contains these words too.
    """
    text = CONF_PY.read_text()
    defined = [
        key for key in _MATHJAX_CONFIG_KEYS if re.search(rf"^\s*{key}\s*=", text, re.MULTILINE)
    ]
    assert not defined, (
        f"conf.py defines {defined}. Macros defined through any of these render "
        f"on the site but not on github.com, and the notes are read in both."
    )


def test_the_scan_actually_covers_the_notes():
    """Guard against the glob silently matching nothing."""
    assert len(_notes()) >= 20
