"""The project has four public locations, and exactly one way to install it.

The four:

    https://qscat.org                  the landing page (the project's identity)
    https://vanamartin.github.io/qscat/ the API reference and theory notes
    https://github.com/VanaMartin/qscat source, issues, changelog
    https://data.qscat.org             the published computed artifacts

They are easy to conflate -- for a while the package metadata named GitHub as
the homepage because there was nothing else to name, and prose called
`qscat.com` the project's domain. `qscat.com` is the maintainer's MAILBOX and
stays one; it is not, and never was, the website. So the rule this file
enforces is asymmetric: `@qscat.com` is fine anywhere, a bare `qscat.com` is
not.

`qscat` is repo-only and is not published to PyPI (and will not be until the
citation article is out), so `pip install qscat` names a distribution that does
not exist. The install route is a clone plus `uv sync --all-packages`, with
`--extra plot` for the figure helpers. That defect had six sites when it was
first found -- two documentation pages, four error messages and a docstring in
the shipped library -- and nothing asserted on any of them, which is why this
file exists.

The scan deliberately skips `docs/superpowers/`: those plans and specs are a
dated record of what was decided at the time, and rewriting history to match
the present would destroy the thing they are for.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

LANDING_PAGE = "https://qscat.org"
DOCUMENTATION = "https://vanamartin.github.io/qscat/"
REPOSITORY = "https://github.com/VanaMartin/qscat"
ARTIFACT_STORE = "https://data.qscat.org"

QSCAT_PYPROJECT = REPO_ROOT / "libs/qscat/pyproject.toml"
CITATION = REPO_ROOT / "CITATION.cff"

# The mailbox stays on qscat.com even though the website moved to qscat.org --
# a domain change is not a mailbox change. Both metadata files carry a comment
# saying so; this pins the address itself, so a well-meaning "fix" to match the
# website fails here rather than silently changing where mail goes.
CONTACT_EMAIL = "martin@qscat.com"


def _scanned_files() -> list[Path]:
    """Every file whose prose speaks for the project to an outside reader.

    The shipped library and app sources (their error messages and docstrings
    reach a user at the moment an import fails), the published documentation,
    and the READMEs and metadata a visitor lands on first.
    """
    paths: list[Path] = []
    for pattern in ("libs/qscat/qscat/**/*.py", "apps/qscat-run/qscat_run/**/*.py"):
        paths.extend(REPO_ROOT.glob(pattern))
    paths.extend(
        p
        for p in REPO_ROOT.glob("docs/**/*.md")
        if "superpowers" not in p.relative_to(REPO_ROOT).parts
    )
    paths.extend(
        REPO_ROOT / name
        for name in (
            "README.md",
            "CONTRIBUTING.md",
            "CITATION.cff",
            "libs/qscat/README.md",
            "libs/qscat/pyproject.toml",
            "apps/qscat-run/README.md",
            "apps/qscat-run/pyproject.toml",
        )
    )
    return sorted(set(paths))


# `pip install qscat`, `pip install "qscat[plot]"`, `pip install 'qscat[mumps]'`
# -- any quoting, any extra. A bare `pip install` of something else (the CI
# wheel smoke test installs `dist/*.whl`) is not this project's install
# guidance and is not matched.
_PIP_INSTALL_QSCAT_RE = re.compile(r"""pip install\s+["']?qscat""")

# A bare `qscat.com`: the website claim. `martin@qscat.com` is the mailbox and
# is explicitly allowed, so the match requires a non-`@` character (or the
# start of the line) in front.
_QSCAT_COM_WEBSITE_RE = re.compile(r"(?<!@)\bqscat\.com\b")


def find_pypi_install_guidance(text: str) -> list[str]:
    """Return every `pip install qscat...` in `text`; none may exist."""
    return [m.group(0) for m in _PIP_INSTALL_QSCAT_RE.finditer(text)]


def find_qscat_com_as_website(text: str) -> list[str]:
    """Return every `qscat.com` in `text` that is not part of an email address."""
    return [m.group(0) for m in _QSCAT_COM_WEBSITE_RE.finditer(text)]


# --- detector unit tests -------------------------------------------------


def test_find_pypi_install_guidance_flags_the_plain_form():
    assert find_pypi_install_guidance("pip install qscat  # core") == ["pip install qscat"]


@pytest.mark.parametrize("quote", ['"', "'"])
def test_find_pypi_install_guidance_flags_a_quoted_extra(quote: str):
    text = f"the plot extra (`pip install {quote}qscat[plot]{quote}`)"
    assert find_pypi_install_guidance(text) == [f"pip install {quote}qscat"]


def test_find_pypi_install_guidance_ignores_an_unrelated_pip_install():
    """The CI wheel smoke test pip-installs a built artifact; that is not this."""
    assert find_pypi_install_guidance("pip install dist/*.whl") == []


def test_find_pypi_install_guidance_passes_the_uv_route():
    assert find_pypi_install_guidance("uv sync --all-packages --extra plot") == []


def test_find_qscat_com_as_website_flags_a_bare_domain():
    text = "the project's domain is `qscat.com`."
    assert find_qscat_com_as_website(text) == ["qscat.com"]


def test_find_qscat_com_as_website_allows_the_mailbox():
    assert find_qscat_com_as_website("reported to martin@qscat.com") == []


def test_find_qscat_com_as_website_allows_qscat_org():
    text = "https://qscat.org and https://data.qscat.org"
    assert find_qscat_com_as_website(text) == []


# --- metadata ------------------------------------------------------------


def test_package_urls_name_all_four_locations_distinctly():
    urls = tomllib.loads(QSCAT_PYPROJECT.read_text())["project"]["urls"]
    assert urls["Homepage"] == LANDING_PAGE
    assert urls["Documentation"] == DOCUMENTATION
    assert urls["Repository"] == REPOSITORY
    assert urls["Issues"].startswith(REPOSITORY)
    assert urls["Changelog"].startswith(REPOSITORY)


def test_package_contact_email_stays_on_the_mailbox_domain():
    project = tomllib.loads(QSCAT_PYPROJECT.read_text())["project"]
    assert [a["email"] for a in project["authors"]] == [CONTACT_EMAIL]
    assert [m["email"] for m in project["maintainers"]] == [CONTACT_EMAIL]


def test_citation_separates_the_landing_page_from_the_source():
    """In CFF, `url` is the landing page and `repository-code` is the source."""
    cff = yaml.safe_load(CITATION.read_text())
    assert cff["url"] == LANDING_PAGE
    assert cff["repository-code"] == REPOSITORY
    assert cff["preferred-citation"]["url"] == LANDING_PAGE
    assert cff["preferred-citation"]["repository-code"] == REPOSITORY


def test_citation_contact_email_stays_on_the_mailbox_domain():
    cff = yaml.safe_load(CITATION.read_text())
    assert [a["email"] for a in cff["authors"]] == [CONTACT_EMAIL]
    assert [a["email"] for a in cff["preferred-citation"]["authors"]] == [CONTACT_EMAIL]


# --- tree scan -----------------------------------------------------------


@pytest.mark.parametrize("path", _scanned_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_file_recommends_no_pypi_package(path: Path):
    found = find_pypi_install_guidance(path.read_text())
    assert not found, (
        f"{path.relative_to(REPO_ROOT)} tells a reader to {found}, but qscat is "
        f"repo-only and no such distribution exists. The install route is "
        f"`git clone {REPOSITORY}` then `uv sync --all-packages` "
        f"(`--extra plot` for the figure helpers)."
    )


@pytest.mark.parametrize("path", _scanned_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_file_does_not_call_qscat_com_the_website(path: Path):
    found = find_qscat_com_as_website(path.read_text())
    assert not found, (
        f"{path.relative_to(REPO_ROOT)} names {found} outside an email address. "
        f"qscat.com is the maintainer's mailbox; the project's website is "
        f"{LANDING_PAGE}."
    )


def test_the_scan_actually_covers_the_tree():
    """Guard against the globs silently matching nothing."""
    scanned = _scanned_files()
    assert len(scanned) >= 50
    names = {p.relative_to(REPO_ROOT).as_posix() for p in scanned}
    assert "docs/getting-started.md" in names
    assert "libs/qscat/qscat/viz/plot.py" in names


# --- published guidance --------------------------------------------------


def test_getting_started_teaches_the_clone_plus_uv_route():
    text = (REPO_ROOT / "docs/getting-started.md").read_text()
    assert f"git clone {REPOSITORY}" in text
    assert "uv sync --all-packages" in text
    assert "--extra plot" in text


def test_the_front_pages_link_the_landing_page_and_the_artifact_store():
    for page in ("README.md", "docs/index.md"):
        text = (REPO_ROOT / page).read_text()
        assert LANDING_PAGE in text, f"{page} does not link {LANDING_PAGE}"
        assert ARTIFACT_STORE in text, f"{page} does not link {ARTIFACT_STORE}"
