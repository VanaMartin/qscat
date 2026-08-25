"""Task 3 (sub-project #9): tests for `validation.n2.ti_curve` and the
generic `projects.n2_2d_cross_section.cross_section_plot.plot_cross_sections`.

Both solving tests run on the PRODUCTION working grid; what is kept small is
the number of solves, not the deck:

- `test_v4_dense_curve_matches_houfek_at_anchors` needs the REAL,
  physically-converged `working_tgrid()` (it is checking numbers against
  Houfek's golden data), so instead of shrinking the grid it shrinks the
  *energy count*: `E_grid` is just the 3 distinct anchor energies
  (0.02, 0.1, 0.2 Ha), not a dense sweep -- 3 sparse-LU solves on the
  working grid (measured 26.6 s), not the minutes a full dense curve costs.
- `test_compute_ti_curve_shape_and_physical` checks only
  shape/dtype/realness/non-negativity, which does not depend on convergence,
  and is cheaper (measured 9.5 s) only because it asks for fewer channels --
  it too builds `_WORKING_TGRID` (see its comment below). Do not read its
  name as a promise of a toy deck.

V4 tolerance: gated anchors (clear of their own vibrational threshold, see
`validation.n2.reference.ANCHOR_MARGIN_HA`) use `GATED_RTOL` -- the same
*tight*, differential-oracle bound `validation/n2/exact2d.py`'s anchors are
already gated at (this module and `exact2d.py` share `_build_system`
verbatim, so the same solver/grid/vib-state combination is exercised here,
just through `ti_curve.compute_ti_curve`'s driver instead of `exact2d`'s).
The two DOCUMENTED-LIMITED anchors -- elastic (E=0.2, v'=0) and
near-threshold (E=0.02, v'=1) -- use the loose cross-model
`reference.ANCHOR_FACTOR` band instead, per the brief: this module is not
re-deriving the GATED/DOCUMENTED-LIMITED split, it reuses
`exact2d.compute_exact2d_results()`'s already-computed `.gated` classification
verbatim (same rationale `exact2d.py` itself gives for reusing `cross_section.py`'s
classification).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_plot import plot_cross_sections
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from validation.n2 import reference
from validation.n2.exact2d import GATED_RTOL, compute_exact2d_results
from validation.n2.ti_curve import compute_ti_curve, houfek_reference

import pytest

# Both solving tests below consume the production WORKING_GRID (N=26857;
# measured 26.6 s + 9.5 s), and the Houfek anchor gate's GATED_RTOL is
# calibrated to that grid's round-off floor, so the deck cannot shrink
# without the gate passing while measuring nothing -- the same
# docs/adr/0005 point 7 case as validation/n2/test_anchors.py, and the same
# resolution: the whole module is `slow`. It still runs in the Docker `test`
# image and on a `validate:n2` (or `validate:all`) pull-request label.
pytestmark = pytest.mark.slow

# These are `convergence.WORKING_GRID`'s parameters -- this is the PRODUCTION
# deck (N=26857), not a reduced one, despite being built by hand here. The
# only nominal difference, nuc_r_max 22.0 vs 20.0, changes the length of the
# ECS tail elements and not the nuclear point count (see WORKING_GRID's own
# comment on that axis), so both spellings give the identical 107 x 251 grid.
# Built explicitly rather than via `working_tgrid()` so this module states the
# parameters it exercises; if it is ever made cheap, shrink it deliberately
# and rename it, do not assume it is already small.
_WORKING_TGRID = TensorGrid(
    [
        n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)


def test_compute_ti_curve_shape_and_physical() -> None:
    E_grid = np.array([0.05, 0.1, 0.15])
    vprimes = [0, 1]
    sigma = compute_ti_curve(E_grid, vprimes, tgrid=_WORKING_TGRID)

    assert sigma.shape == (len(E_grid), len(vprimes))
    assert np.all(np.isfinite(sigma))
    assert np.all(sigma >= 0.0)
    assert sigma.dtype == np.float64


def test_v4_dense_curve_matches_houfek_at_anchors() -> None:
    anchor_rows = reference.anchors()  # [(e_row, channel, sigma_houfek), ...]
    gated_by_key = {(r.energy_ha, r.channel): r.gated for r in compute_exact2d_results()}

    e_rows = sorted({e for e, _ch, _s in anchor_rows})
    E_grid = np.asarray(e_rows, dtype=np.float64)
    assert len(E_grid) == 3  # 0.02, 0.1, 0.2 Ha -- documents the expected anchor spread

    vprimes = [0, 1, 2, 3]
    sigma = compute_ti_curve(E_grid, vprimes)  # working_tgrid() default

    for e_row, channel, sigma_houfek in anchor_rows:
        i = int(np.argmin(np.abs(E_grid - e_row)))
        j = vprimes.index(channel)
        sigma_computed = float(sigma[i, j])
        ratio = sigma_computed / sigma_houfek if sigma_houfek != 0 else float("inf")

        if gated_by_key[(e_row, channel)]:
            dev = abs(ratio - 1.0)
            assert dev < GATED_RTOL, (
                f"anchor (E={e_row}, v'={channel}): ratio {ratio:.6f} deviates "
                f"{dev:.3e} from 1, outside GATED_RTOL={GATED_RTOL}"
            )
        else:
            assert 1.0 / reference.ANCHOR_FACTOR <= ratio <= reference.ANCHOR_FACTOR, (
                f"documented-limited anchor (E={e_row}, v'={channel}): ratio "
                f"{ratio:.6f} outside the loose factor-{reference.ANCHOR_FACTOR} band"
            )


def test_houfek_reference_slices_requested_channels() -> None:
    vprimes = [0, 2]
    e_ref, sigma_ref = houfek_reference(vprimes)
    assert e_ref.ndim == 1
    assert sigma_ref.shape == (e_ref.shape[0], len(vprimes))
    assert np.all(sigma_ref >= 0.0)


def test_plot_cross_sections_writes_nonempty_png(tmp_path: Path) -> None:
    E_grid = np.linspace(0.01, 0.2, 10)
    sigma = np.stack([np.full(10, 1.0), np.full(10, 0.1)], axis=1)
    e_ref = np.linspace(0.0, 0.2, 20)
    sigma_ref = np.stack([np.full(20, 1.1), np.full(20, 0.09)], axis=1)

    out_path = tmp_path / "test_plot.png"
    plot_cross_sections(
        E_grid,
        sigma,
        channels=[0, 1],
        reference=(e_ref, sigma_ref),
        thresholds={1: 0.05},
        title="test plot",
        path=out_path,
    )

    assert out_path.exists()
    assert out_path.stat().st_size > 0
