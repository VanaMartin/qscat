"""sigma_DA refinement stability check on eMoScat's per-molecule nuclear
deck (F2).

An earlier quadrature ladder (docs/physics/diatomic-ve-cross-sections.md,
"The discretisation must be per-molecule") showed sigma_DA does NOT converge on the shared
N2-style nuclear grid by raising quadrature alone -- the K_R~58 dissociation
wave in the [2.7, 10.7] bohr outer region is under-resolved by the grid's
ELEMENT SIZE, not its quadrature order. eMoScat's per-molecule
deck (`MoleculeConfig.da_grid`, sub-project Task 6) fixes the element size
directly (~0.2 bohr/element there). This benchmark is the refinement
convergence evidence on that deck.

DEVIATION FROM THE ORIGINAL PLAN: a true h-refinement (doubling the OUTER
real nuclear segment's element count, ~197k unknowns) was attempted here and
did not complete on this laptop (SuperLU on scipy; no output, no partial
result -- OOM-killed before any print flushed). Retrying the same heavy solve
would not produce new information, so the refined variant below instead
raises the nuclear QUADRATURE order on the SAME eMoScat elements
(`nuc_quad + 2`, ~132x1126 ~= 149k unknowns) -- lighter than the h-refined
element-doubling (~197k) while still probing resolution beyond the eMoScat
deck. A genuine h-refinement sweep (doubled elements, or MUMPS instead of
SuperLU for the factorization) needs Docker + the `qscat[mumps]` extra; see
`docs/physics/mumps-sparse-backend.md`. Run via:

    uv run python -m benchmarks.da_grid_stability

MEASURE -- this script asserts nothing; it reports sigma_DA + the relative
difference between the two resolutions so the numbers can go in the docs.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS, MoleculeConfig

# How much to raise the nuclear quadrature order for the "refined" variant
# (see module docstring: a lighter substitute for h-refining the element
# size, which OOM'd on this laptop).
_QUAD_BUMP = 2


def _da_grid_with_quad(cfg: MoleculeConfig, quadrature: int) -> TensorGrid:
    return TensorGrid(
        [
            electronic_grid(r_max=cfg.e_r_max, order=cfg.e_order, n_complex=cfg.e_n_complex),
            segmented_grid(
                cfg.nuc_real, cfg.nuc_complex, angle_deg=cfg.nuc_angle, quadrature=quadrature
            ),
        ]
    )


def stability(cfg: MoleculeConfig, E: npt.NDArray[np.float64]) -> dict[str, object]:
    """sigma_DA(E) on `cfg.da_grid()` (the eMoScat deck, quadrature
    `cfg.nuc_quad`) vs the same elements at `cfg.nuc_quad + _QUAD_BUMP` +
    their relative difference -- the convergence evidence for the docs (see
    module docstring for why this replaces a full h-refinement here).
    """
    tg_base = cfg.da_grid()
    tg_fine = _da_grid_with_quad(cfg, cfg.nuc_quad + _QUAD_BUMP)

    eps_b, chi_b = vibrational_states(tg_base.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    eps_f, chi_f = vibrational_states(tg_fine.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)

    sigma_base = da_cross_section(tg_base, cfg.model, eps_b, chi_b, 0, E)[:, 0]
    sigma_fine = da_cross_section(tg_fine, cfg.model, eps_f, chi_f, 0, E)[:, 0]
    rel_diff = np.abs(sigma_fine - sigma_base) / np.maximum(sigma_base, 1e-30)

    return {
        "E": E,
        "n_nuclear_base": tg_base.grids[1].n,
        "n_nuclear_fine": tg_fine.grids[1].n,
        "sigma_base": sigma_base,
        "sigma_fine": sigma_fine,
        "rel_diff": rel_diff,
    }


def main() -> None:
    cfg = CONFIGS["F2"]
    E = np.array([0.02, 0.03, 0.04])
    result = stability(cfg, E)
    print(
        f"F2 sigma_DA(E) quadrature-refinement stability: eMoScat deck "
        f"(nuclear n={result['n_nuclear_base']}, quad={cfg.nuc_quad}) vs "
        f"quad+{_QUAD_BUMP} (nuclear n={result['n_nuclear_fine']})"
    )
    sigma_base = result["sigma_base"]
    sigma_fine = result["sigma_fine"]
    rel_diff = result["rel_diff"]
    assert isinstance(sigma_base, np.ndarray)
    assert isinstance(sigma_fine, np.ndarray)
    assert isinstance(rel_diff, np.ndarray)
    for e, sb, sf, rel in zip(E, sigma_base, sigma_fine, rel_diff, strict=True):
        print(f"  E={e:.3f} Ha  base={sb:.6e}  fine={sf:.6e}  rel_diff={rel:.2%}")


if __name__ == "__main__":
    main()
