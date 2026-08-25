"""Per-molecule eMoScat DECK definitions for NO and F2.

Since the exact-2D VE/DA/LCP *cross-section curves* now run through
`apps/qscat-run` (config-driven, `qscat_run.presets` carries the same decks),
this module's remaining job is narrow: it is the source of the fine per-molecule
NUCLEAR deck (`da_grid`) that the **discretisation tuner** reads as its
reproduce-and-beat reference (`validation/tuning/{calibrate,test_emoscat_decks,
test_resonance_aware}.py`), plus the electronic-grid parameters the tuner needs
to place F2's anion bound state.

The decks here are byte-identical to `qscat_run.presets`' F2/NO decks
(`_f2_nuc_grid`/`_no_nuc_grid` + the electronic `r_max=16/order=8/n_complex=6`);
`test_da_grid.py::test_diatomic_decks_match_presets` locks that invariant so the
two copies (which layering keeps separate -- `qscat_run` must not import
`validation`, and the tuner should not reach into the app's grid internals)
never drift.
"""

from __future__ import annotations

from dataclasses import dataclass

from qscat.core.grids import electronic_grid, segmented_grid
from qscat.dvr import TensorGrid
from qscat.model import F2, NO, DiatomicResonanceModel

__all__ = ["CONFIGS", "MoleculeConfig"]


@dataclass(frozen=True)
class MoleculeConfig:
    """One molecule's model + the eMoScat FEM-DVR-ECS grid parameters.

    `e_r_max`/`e_order`/`e_n_complex` size the electronic grid; the
    `nuc_*` fields are the eMoScat per-molecule nuclear deck (transcribed
    verbatim from `reference/eMoScat/input/{NO,F2}/grids.txt`, 2nd/nuclear
    declaration) that `da_grid()` builds via `segmented_grid`.
    """

    name: str
    model: DiatomicResonanceModel
    # electronic FEM-DVR-ECS grid
    e_r_max: float
    e_order: int
    e_n_complex: int
    n_vib: int  # number of neutral vibrational states to resolve
    # eMoScat per-molecule nuclear deck: (n_elements, endpoint) segment pairs.
    nuc_real: tuple[tuple[int, float], ...]
    nuc_complex: tuple[tuple[int, float], ...]
    nuc_angle: float
    nuc_quad: int

    def da_grid(self) -> TensorGrid:
        """Electronic (r_max=e_r_max) x the eMoScat fine nuclear deck.

        The fine nuclear grid the fast K_R~58 dissociation wave needs (the
        shared N2-style grid under-resolves it); the tuner's reproduce-and-beat
        reference deck.
        """
        return TensorGrid(
            [
                electronic_grid(r_max=self.e_r_max, order=self.e_order, n_complex=self.e_n_complex),
                segmented_grid(
                    self.nuc_real,
                    self.nuc_complex,
                    angle_deg=self.nuc_angle,
                    quadrature=self.nuc_quad,
                ),
            ]
        )


CONFIGS: dict[str, MoleculeConfig] = {
    "NO": MoleculeConfig(
        name="NO",
        model=NO,
        e_r_max=16.0,
        e_order=8,
        e_n_complex=6,
        n_vib=4,
        nuc_real=((1, 1.0), (1, 1.6), (37, 9.0)),
        nuc_complex=((1, 9.25), (1, 10.0), (1, 12.0), (4, 42.0)),
        nuc_angle=45.0,
        nuc_quad=14,
    ),
    "F2": MoleculeConfig(
        name="F2",
        model=F2,
        e_r_max=16.0,
        e_order=8,
        e_n_complex=6,
        n_vib=4,
        nuc_real=((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)),
        nuc_complex=(
            (1, 10.8),
            (1, 11.0),
            (1, 11.5),
            (1, 12.5),
            (1, 14.0),
            (1, 18.0),
            (4, 30.0),
            (2, 101.0),
        ),
        nuc_angle=35.0,
        nuc_quad=14,
    ),
}
