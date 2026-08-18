"""Tests for `plot_resonance_levels`' shift panel, which used to lie by default.

The right-hand panel differences two level series and labels the result a shift
in meV. It paired them by sorted index, which is correct only when both sets are
complete and ordered alike -- and that assumption fails exactly where the
physics is interesting. On H2+ two BO levels 20 uHa apart correspond to exact
poles 154 uHa apart, so index pairing crosses them and reports two shifts
belonging to neither level.

These tests pin the fix: an explicit `pairing` is honoured, and it is honoured
in the caller's own order rather than a sorted one.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("matplotlib")

from qscat.core import plot_resonance_levels  # noqa: E402

# Two series in the crossing geometry, deliberately NOT in matching order:
# `exact[0]` belongs with `bo[1]` and vice versa.
_BO = np.array([-0.1000 - 1e-6j, -0.1000200 - 1e-6j], dtype=np.complex128)
_EXACT = np.array([-0.1000210 - 2e-6j, -0.0999960 - 2e-6j], dtype=np.complex128)


def _render(tmp_path, name, *, bo=_BO, exact=_EXACT, **kw):
    path = tmp_path / f"{name}.png"
    plot_resonance_levels({"BO": bo, "exact": exact}, path=path, baseline="BO", **kw)
    assert path.exists() and path.stat().st_size > 0
    return path.read_bytes()


def test_an_explicit_pairing_changes_the_shift_panel(tmp_path):
    """The whole point: the same inputs must render differently once paired.

    Index pairing differences `exact[0]-BO[0]` and `exact[1]-BO[1]`; the correct
    pairing differences `exact[0]-BO[1]` and `exact[1]-BO[0]`. If `pairing` were
    ignored, these two renders would be byte-identical.
    """
    default = _render(tmp_path, "default")
    paired = _render(tmp_path, "paired", pairing={"exact": [(0, 1), (1, 0)]})
    assert default != paired


def test_a_pairing_that_reproduces_the_default_renders_the_same(tmp_path):
    """A sanity anchor: the mechanism is not simply perturbing every render.

    Inputs must be pre-sorted for this comparison, because supplying a
    `pairing` also suppresses the internal sort -- the pairing indexes the
    caller's own ordering, so re-sorting would silently re-map every pair it
    names.
    """
    bo, exact = np.sort_complex(_BO), np.sort_complex(_EXACT)
    default = _render(tmp_path, "d2", bo=bo, exact=exact)
    same = _render(tmp_path, "s2", bo=bo, exact=exact, pairing={"exact": [(0, 0), (1, 1)]})
    assert default == same


def test_a_pairing_suppresses_the_internal_sort(tmp_path):
    """Unsorted input renders differently with and without a pairing.

    Without `pairing` the series are sorted before anything is drawn; with it
    they are used as given. A caller who computed indices against their own
    array must get those indices back.
    """
    identity = {"exact": [(0, 0), (1, 1)]}
    unsorted = _render(tmp_path, "u1", pairing=identity)
    presorted = _render(
        tmp_path, "u2", bo=np.sort_complex(_BO), exact=np.sort_complex(_EXACT), pairing=identity
    )
    assert unsorted != presorted


def test_a_series_absent_from_the_pairing_is_skipped_not_index_paired(tmp_path):
    """Silence must mean 'no shift shown', never 'fall back to index pairing'.

    A quiet fallback is how the wrong pairing survived: it produced a plausible
    panel with no indication that nothing had been asserted about the matching.
    """
    with pytest.warns(UserWarning, match="No artists with labels"):
        empty = _render(tmp_path, "empty", pairing={})
    default = _render(tmp_path, "d3")
    assert empty != default


def test_plot_still_works_without_a_baseline(tmp_path):
    path = tmp_path / "single.png"
    plot_resonance_levels({"BO": _BO}, path=path)
    assert path.exists()


def test_an_out_of_range_pairing_index_is_an_error_not_a_silent_drop(tmp_path):
    with pytest.raises(IndexError):
        _render(tmp_path, "bad", pairing={"exact": [(0, 5)]})
