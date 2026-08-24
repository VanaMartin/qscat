"""GATE: the exact-2D NO sigma_DA against the ONE published NO DA curve.

Vana 2017 thesis, **Fig. 3.14, p. 46** (`reference/literature/vana-2017-thesis.md`)
is the only published dissociative-attachment cross section for the NO-like
model. Its bottom panel plots `VE: 0 -> DA` on a LINEAR axis in units of
1e-9 a0^2 over 0.170-0.200 Ha, comparing a time-independent 2-D reference
against three time-dependent S-matrix routes and the LCP (the LCP curve is
drawn scaled by 1e-5 to share the axis -- i.e. the LCP over-predicts by five
orders, which is itself the figure's message).

Read off that panel: the curve is zero below threshold, jumps to its peak of
**~1.35e-9 a0^2** at the first sampled energy above it (~0.1720 Ha), and then
decays smoothly and monotonically, passing ~1e-9 near 0.1730 and falling
below 1e-10 by ~0.1875.

THIS GATE EXISTS BECAUSE THE ORACLE WAS WRONG. `da_cross_section` used to
extract sigma from the post-form volume-integral T-matrix
`<phi_e F_K | V_DR | Psi+>`, which is formally exact but numerically unusable
here: NO's sigma_DA is ~1e-9 a0^2 while the T-matrix integrand sums to ~2.6,
so the answer is the residue of a ~1e6-fold cancellation. Any part of the
integrand that has not decayed at the edge of the integration region then
dominates. On the shipped decks it did, twice over -- the electronic real
region (16 bohr) and the nuclear one (`R_inf` = 9.0 bohr, where NO's Morse
`v0` is still -1.0e-5 Ha) -- and the published sigma came out 1e4 to 1e7 times
too large. F2 escaped because its sigma_DA is O(1) a0^2 and its cancellation
is only ~3-fold. See `docs/physics/diatomic-ve-cross-sections.md`.

The replacement is the outgoing-flux extraction (the same one
`qscat.core.lcp.lcp_da_cross_section` uses), which needs no cancellation at
all and is invariant under both of those box sizes -- that invariance is
gated by `test_no_da_is_invariant_under_the_electronic_box` below.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from qscat.core.dissociation import da_cross_section
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS

# The read-off values live in ONE place, `data/vana-2017-fig3.14-no-da.dat`
# (which also carries the provenance and the reading-uncertainty caveat), so
# this gate and the committed figure's overlay
# (`validation/diatomic/da_figure.py`) can never quote the panel differently.
_REFERENCE = np.loadtxt(Path(__file__).parent / "data/vana-2017-fig3.14-no-da.dat")

# Per-point tolerance, deliberately generous and deliberately NOT uniform: the
# panel is a ~4 cm LINEAR axis, so its peak is readable to ~10% while its tail,
# at a few percent of full scale, is only good to a factor of ~2.5. What the
# gate has to catch is the four-to-seven-order error that was there before, and
# these bands do that with enormous margin while claiming no precision the
# figure cannot support.
_FACTORS: tuple[float, ...] = (1.4, 1.4, 1.6, 2.5, 2.5, 3.0)

THESIS_FIG_314: tuple[tuple[float, float, float], ...] = tuple(
    (float(e), float(sigma), factor)
    for (e, sigma), factor in zip(_REFERENCE, _FACTORS, strict=True)
)


def _no_grid(e_r_max: float = 16.0) -> TensorGrid:
    """NO's eMoScat nuclear deck x an electronic grid of the given extent."""
    cfg = CONFIGS["NO"]
    return TensorGrid(
        [
            electronic_grid(r_max=e_r_max, order=cfg.e_order, n_complex=cfg.e_n_complex),
            segmented_grid(
                cfg.nuc_real,
                cfg.nuc_complex,
                angle_deg=cfg.nuc_angle,
                quadrature=cfg.nuc_quad,
            ),
        ]
    )


def _sigma_da(tg: TensorGrid, energies: np.ndarray) -> np.ndarray:
    cfg = CONFIGS["NO"]
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    return np.asarray(da_cross_section(tg, cfg.model, eps, chi, 0, energies), dtype=np.float64)[
        :, 0
    ]


@pytest.mark.slow
def test_no_da_matches_the_published_thesis_curve() -> None:
    """sigma_DA(NO) lands on Vana 2017 Fig. 3.14 at every sampled energy."""
    energies = np.array([e for e, _, _ in THESIS_FIG_314])
    sigma = _sigma_da(_no_grid(), energies)
    for (e, published, factor), got in zip(THESIS_FIG_314, sigma, strict=True):
        assert published / factor <= got <= published * factor, (
            f"sigma_DA(NO, E={e}) = {got:.4e} a0^2 is outside a factor {factor} "
            f"of the Fig. 3.14 value {published:.2e} a0^2"
        )


@pytest.mark.slow
def test_no_da_is_invariant_under_the_electronic_box() -> None:
    """The extraction must not depend on where the electronic real region ends.

    This is the property the old volume-integral extraction did NOT have: over
    exactly this 16 -> 48 bohr change it moved sigma_DA by four orders of
    magnitude. Measured with the flux extraction: 4-digit agreement.
    """
    energies = np.array([0.175, 0.185, 0.200])
    small = _sigma_da(_no_grid(e_r_max=16.0), energies)
    large = _sigma_da(_no_grid(e_r_max=48.0), energies)
    assert np.allclose(small, large, rtol=2e-3, atol=0.0), (
        f"sigma_DA moved with the electronic box: {small} vs {large}"
    )
