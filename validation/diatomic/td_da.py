"""Time-dependent dissociative-attachment (DA) three-way validation for F2/NO.

`qscat.core.time_dependent.td_da_cross_sections_all` (the nuclear-axis
`Flux`/`Dirac`/`TannorWeeks` extractors, `qscat.core.td_extractors`) vs the
exact-2D TI `da_cross_section` oracle (`validation.diatomic.da_curves`) --
the SAME differential-test pattern `n2_2d_td_cross_section` established for
VE, generalized to DA. See `docs/physics/td-da.md` for the full write-up.

**KEY LESSON**: TD-DA reuses the SAME grid + incident wavepacket shape as
the electronic-axis Tannor-Weeks VE method -- a LARGE electronic LAUNCH box
(`r0` well inside `r_max`) -- but NOT the TI `MoleculeConfig.da_grid()`'s
electronic box: that box (`e_r_max=16`) is sized for the DRIVEN solve, not
for launching a TD wavepacket, and was found (controller, SP2 Task 2) to
diverge ~1e6x when the incident packet was placed off-box. `td_launch_grid`
below builds a dedicated launch-box electronic grid (`r_max=25`, matching
`libs/qscat/tests/test_td_extractors.py`'s validated `@slow` gate config)
paired with the UNCHANGED eMoScat fine nuclear deck (`da_grid().grids[1]`,
which already resolves the fast K_R dissociation wave -- a coarse nuclear
grid reads sigma~0, the surface never sees the outgoing flux).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid
from qscat.core.time_dependent import td_da_cross_sections_all
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid

from validation.diatomic.config import CONFIGS, MoleculeConfig

__all__ = ["td_launch_grid", "WP_IN", "WP_OUT_NUCLEAR", "compute_td_da_three_way"]

# The TD launch-box electronic grid -- SAME as
# `test_td_extractors.py::test_nuclear_flux_da_converges_to_ti_oracle`'s
# validated `@slow` gate config: large enough (r_max=25) that the incident
# r0=12 launches cleanly, well inside the real region (unlike `da_grid()`'s
# electronic box, sized for the driven TI solve, not for launching a TD
# wavepacket).
_LAUNCH_R_MAX = 25.0
_LAUNCH_ORDER = 6
_LAUNCH_N_COMPLEX = 3
_LAUNCH_ANGLE_DEG = 40.0

# eMoScat's F2 nuclear test packet (controller-validated, SP2): a narrow
# (wide-K) packet placed inward of the surface/position analysis points,
# resolving the K_R~72 dissociation flux wave.
WP_IN = {"r0": 12.0, "p0": -0.5, "sigma": 3.0}
WP_OUT_NUCLEAR = {"r0_out": 8.0, "p0_out": 72.0, "sigma_out": 0.07}


def td_launch_grid(cfg: MoleculeConfig) -> TensorGrid:
    """Launch-box electronic (r_max=25) x `cfg`'s eMoScat fine nuclear deck."""
    elec = electronic_grid(
        r_max=_LAUNCH_R_MAX,
        order=_LAUNCH_ORDER,
        n_complex=_LAUNCH_N_COMPLEX,
        angle_deg=_LAUNCH_ANGLE_DEG,
    )
    nuc = cfg.da_grid().grids[1]
    return TensorGrid([elec, nuc])


def _real_index_near(tg: TensorGrid, R: float) -> int:
    """Nearest REAL-region (unscaled) nuclear DVR index to `R` (bohr)."""
    return tg.grids[1].real_index_near(R)


def compute_td_da_three_way(
    cfg: MoleculeConfig,
    E: npt.ArrayLike,
    *,
    n_steps: int = 1800,
    dt: float = 1.0,
    surface_R: float = 6.0,
    position_R: float = 6.0,
    wp_in: dict[str, float] | None = None,
    wp_out: dict[str, float] | None = None,
    n_channels: int = 1,
) -> dict[str, npt.NDArray[np.float64]]:
    """`{"flow":, "delta":, "tw":}` sigma_DA(E) (bohr^2, per anion channel)
    from ONE shared propagation on `cfg`'s TD launch grid (`td_launch_grid`).

    `surface_R`/`position_R` (bohr) are converted to the nearest real-region
    nuclear DVR index for the `Flux`/`Dirac` extractors; `wp_out` (nuclear
    R-axis Gaussian test packet) drives `TannorWeeks`. Defaults are the
    controller-validated eMoScat F2 config (module docstring).
    """
    wp_in = wp_in if wp_in is not None else WP_IN
    wp_out = wp_out if wp_out is not None else WP_OUT_NUCLEAR
    tg = td_launch_grid(cfg)
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    surface = _real_index_near(tg, surface_R)
    position = _real_index_near(tg, position_R)
    return td_da_cross_sections_all(
        tg,
        cfg.model,
        eps,
        chi,
        0,
        E,
        dt=dt,
        n_steps=n_steps,
        wp_in=wp_in,
        surface=surface,
        position=position,
        wp_out=wp_out,
        n_channels=n_channels,
    )


def main() -> None:
    """Run the F2/NO three-way TD-DA validation and print a summary.

    NOT invoked by the test suite -- a full run per molecule is ~10 min
    (see `docs/physics/td-da.md`); run manually via `uv run python -m
    validation.diatomic.td_da`.
    """
    from qscat.core.dissociation import da_cross_section

    for name in ("F2", "NO"):
        cfg = CONFIGS[name]
        e_probe = np.array([0.03, 0.04])
        tg = td_launch_grid(cfg)
        eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
        sigma_ti = np.ravel(da_cross_section(tg, cfg.model, eps, chi, 0, e_probe))
        sigma_td = compute_td_da_three_way(cfg, e_probe)
        print(f"{name}: sigma_ti = {sigma_ti}")
        for method, sigma in sigma_td.items():
            ratio = np.ravel(sigma) / sigma_ti
            print(f"{name} {method}: sigma = {np.ravel(sigma)}  ratio = {ratio}")


if __name__ == "__main__":
    main()
