"""The resonance curve each accessible partial wave supports on its own.

The coupled model's diagonal blocks are single-channel problems, one per `l`,
differing in TWO ways: the centrifugal barrier `l(l+1)/2r^2`, and -- once the
well is anisotropic -- the angular projection `V_ll` of the interaction. This
walks the resonance of each block separately and sets it beside the shipped
`l = 1` model, which answers a question the coupled cross section cannot: is
there anything resonant in the higher waves at all, or are they inert channels
that the `l = 1` resonance merely leaks into and back out of?

At `s = 0` the projection is `l`-independent, so the spread across `l` is the
centrifugal barrier ALONE, and `l = 1` is the shipped model exactly. Turning
the anisotropy on adds the projection difference on top.

Run: `python -m validation.coupled.per_ell_curves`
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.core.grids import electronic_grid
from qscat.dvr import kinetic_sparse
from qscat.ecs import match_angle_stable
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import DiagonalChannelModel
from projects.no_coupled_channels.renormalised import RenormalisedTwoCentreWell, ScaleTable

__all__ = ["ELLS", "RESULTS", "channel_poles", "main"]

RESULTS = Path("validation/coupled/results")
ANGLES = (44.0, 52.0)
# The screen's deck, used at BOTH its angles here -- unlike a bound state, a
# resonance needs two rotated spectra to be told from a discretised continuum.
DECK = {"r_max": 16.0, "order": 8, "n_complex": 8}
# Up to l = 5 because the question is whether ANY higher wave binds, and the
# centrifugal barrier only grows with l: l = 5 contributing nothing settles
# l > 5 as well.
ELLS = (1, 2, 3, 4, 5)


def _spectrum(well: object, ell: int, R: float, angle: float) -> npt.NDArray[np.complex128]:
    g = electronic_grid(angle_deg=angle, **DECK)
    r = np.asarray(g.points, dtype=np.complex128)
    surface = DiagonalChannelModel(well, ell).surface(r, R)  # type: ignore[arg-type]
    H = (kinetic_sparse(g, 1.0) + sp.diags(surface, format="csr")).toarray()
    return np.linalg.eigvals(H)


def channel_poles(
    well: object, ell: int, R: float, *, resid_max: float = 1e-5
) -> list[tuple[complex, float]]:
    """Every GENUINE angle-stable pole of the `l` block, over a wide window.

    No seed. A seed is a disambiguator between poles, and using one here would
    presuppose the answer -- the question is whether the `l` block has a
    resonance AT ALL, so the window spans the whole physically reachable range
    and everything angle-stable in it is returned.

    The residual is what separates a pole from an artefact, not proximity to a
    guess: this deck carries a spurious near-threshold state at
    `E_res ~ +0.001`, `Gamma ~ 0.006` whose residual is 7e-4, against 1e-8 for
    the genuine `l = 1` resonance -- four orders apart, so the cut is not
    delicate.
    """
    v0 = complex(np.asarray(NO.v0(R)).ravel()[0])
    a = _spectrum(well, ell, R, ANGLES[0])
    b = _spectrum(well, ell, R, ANGLES[1])
    window = (v0.real - 0.5, v0.real + 4.0, -0.6, 0.05)
    try:
        energies, residuals, _ = match_angle_stable(a, b, window, atol=1e-3)
    except Exception:
        return []
    out = [
        (complex(e - v0), float(r))
        for e, r in zip(energies, residuals, strict=True)
        if r < resid_max and (e - v0).imag < 1e-9
    ]
    return sorted(out, key=lambda z: z[0].real)


def main(results: Path = RESULTS) -> dict[str, object]:
    d = np.load(results / "f_table_s03.npz")
    table = ScaleTable(R=np.asarray(d["R"], dtype=float), f=np.asarray(d["f"], dtype=float))
    ones = ScaleTable(R=np.array([0.0, 60.0]), f=np.array([1.0, 1.0]))

    wells = {
        "s=0 (shipped)": RenormalisedTwoCentreWell(base=NO, table=ones, s=0.0, kappa=0.5),
        "s=0.3 renormalised": RenormalisedTwoCentreWell(base=NO, table=table, s=0.3, kappa=0.5),
        "s=0.3 bare": TwoCentreWell(base=NO, s=0.3, kappa=0.5),
    }
    R_list = [1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.4, 4.2]
    report: dict[str, object] = {"ells": list(ELLS), "R": R_list, "curves": {}}

    for name, well in wells.items():
        print(f"\n=== {name} ===")
        rows = {}
        for R in R_list:
            per_l = {}
            bits = []
            for ell in ELLS:
                poles = channel_poles(well, ell, R)
                per_l[ell] = [{"re": e.real, "gamma": -2 * e.imag, "residual": r} for e, r in poles]
                if poles:
                    e, _ = poles[0]
                    bits.append(f"l={ell}: {e.real:+.5f}/{-2 * e.imag:.5f}")
                else:
                    bits.append(f"l={ell}: none")
            print(f"  R={R:5.2f}   " + "   ".join(bits))
            rows[str(R)] = per_l
        report["curves"][name] = rows  # type: ignore[index]

    results.mkdir(parents=True, exist_ok=True)
    (results / "per_ell_curves.json").write_text(json.dumps(report, indent=1))
    print(f"\nwrote {results / 'per_ell_curves.json'}")
    return report


if __name__ == "__main__":
    main()
