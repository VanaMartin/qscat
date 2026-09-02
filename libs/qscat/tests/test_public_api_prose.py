"""The installed prose must describe the callables that actually ship.

`qscat.core`'s module docstring IS the installed API inventory: a reader who
has only the package consults it before reading any signature. That prose can
fall behind the code silently -- nothing executes it -- and the reader then
calls a keyword the shipped callable no longer has and gets a `TypeError`.
The checks below tie the inventory to `inspect.signature`, so a signature
change that leaves the prose stale fails here instead.

Scoped to the dissociative-recombination entry points, where prose and
signature have drifted apart before: the free function
`qscat.core.dr_cross_section` is sigma-only, while the object-API method
`ScatteringProblem.dr_cross_section` keeps its flags -- a distinction the
inventory has to state correctly.
"""

from __future__ import annotations

import importlib.util
import inspect
import re
from collections.abc import Callable
from typing import Any

import qscat.core
from qscat.core import dr_cross_section, dr_solve

# A keyword-shaped switch, the form every detail-bearing flag in this package
# takes (`return_wavefunction`, `store_amplitude`, ...).
_SWITCH = re.compile(r"\b(?:return|store)_[a-z_]+\b")


def _api_bullets(doc: str) -> list[str]:
    """The `Public API:` bullets of a module docstring, one entry per bullet."""
    _, marker, inventory = doc.partition("Public API:")
    assert marker, "qscat.core's docstring must carry a `Public API:` inventory"
    return re.split(r"\n {2}- ", inventory)[1:]


def _dr_prose() -> str:
    """The inventory bullets that describe the DR entry points."""
    bullets = _api_bullets(qscat.core.__doc__ or "")
    dr = [b for b in bullets if "dr_cross_section" in b or "dr_solve" in b]
    assert dr, "qscat.core's inventory must describe the DR entry points"
    return "\n".join(dr)


def _keyword_only(fn: Callable[..., Any]) -> set[str]:
    return {
        name
        for name, p in inspect.signature(fn).parameters.items()
        if p.kind is inspect.Parameter.KEYWORD_ONLY
    }


def test_dr_prose_names_every_keyword_dr_cross_section_accepts() -> None:
    prose = _dr_prose()
    missing = sorted(k for k in _keyword_only(dr_cross_section) if k not in prose)
    assert not missing, (
        f"qscat.core's inventory does not mention {missing}, which "
        f"dr_cross_section accepts; the prose has fallen behind the signature"
    )


def test_dr_prose_names_no_switch_the_shipped_callables_lack() -> None:
    prose = _dr_prose()
    shipped = _keyword_only(dr_cross_section) | _keyword_only(dr_solve)
    phantom = sorted(set(_SWITCH.findall(prose)) - shipped)
    assert not phantom, (
        f"qscat.core's inventory advertises {phantom}, which neither "
        f"dr_cross_section nor dr_solve accepts; a reader following it gets a "
        f"TypeError"
    )


def test_dr_cross_section_is_sigma_only_and_the_prose_routes_details_to_dr_solve() -> None:
    switches = sorted(k for k in _keyword_only(dr_cross_section) if _SWITCH.fullmatch(k))
    assert not switches, (
        f"dr_cross_section grew the switches {switches}; it is documented as "
        f"sigma-only, so either the signature or the inventory is wrong"
    )
    prose = _dr_prose()
    for name in ("dr_solve", "DrResult"):
        assert name in prose, (
            f"the inventory must name `{name}` as the route to wavefunctions, "
            f"amplitudes and other detailed results"
        )


def test_core_prose_promises_no_deprecation_machinery_the_package_dropped() -> None:
    # The package ships no deprecation shims (ADR 0004 point 2: pre-1.0,
    # removal is immediate), so no installed docstring may promise a
    # deprecation cycle a caller could rely on.
    assert importlib.util.find_spec("qscat._deprecation") is None
    assert "deprecat" not in (qscat.core.__doc__ or "").lower(), (
        "qscat.core's docstring promises a deprecation cycle, but the package "
        "carries no deprecation machinery"
    )
