"""Vibrational excitation in the nonlocal resonance model.

PRA 77 Sec. III B splits the VE T-matrix into a background and a resonant
part (Eq. 28):

    T^VE_{vi->vf} = T^bg + T^res
    T^res = <chi_vf| V^{-*}_{dk_f} |Psi_d^+>                     Eq. (31)
    T^bg  = <chi_vf phi^-_{k_f}| V_int |chi_vi J_{k_i}>
            - <chi_vf| V^{-*}_{dk_f} J_{dk_i} |chi_vi>           Eq. (37)
    J_{dk}(R) = Int dr phi_d*(r;R) J_k(r)                        Eq. (38)
    sigma = 4 pi^3 |T^VE|^2 / k_i^2

`Psi_d^+` is the SAME Eq. (52) solution `dissociation.nrm_da_cross_section`
computes, so the resonant part is a contraction over already-gated machinery;
the background part is what is new here, and it is what PRA 77 shows a bare
LCP curve is missing (largest for the broadest resonance).

CONJUGATION -- three citations, three separate jobs; do not let one stand in
for another:

- **Eq. (34)** is the CONDITION selecting the implementable branch: for the
  radial case with a real discrete state, `phi_k^- = (phi_k^+)^*`.
- **Eq. (35)** is why that condition drops the conjugation on the COUPLING
  FACTOR itself, `V^{-*}_{dk} -> V^+_{dk}` (non-conjugated) -- it follows
  from `H_el` Hermiticity in the exact, unscaled theory, combined with
  Eq. (34). The paper's own emphasis: "we can use the matrix element
  `Vdk^+` but *without complex conjugation*."
- **p. 012710-6** is the separate reason the SCALAR PRODUCT itself (the
  c-product pairing of any two DVR coefficient vectors,
  `qscat.linalg.c_product`) carries no conjugation -- under exterior
  complex scaling `P H_el P` is complex symmetric, not Hermitian.

Eq. (34) governs BOTH conjugated objects Eq. (37) prints: the coupling
factor `V^{-*}_{dk_f}` (via Eq. 35, as above) AND the bra `<phi^-_{k_f}|`
in term 1 -- conjugating that bra turns it into a pairing against
`phi^+_{k_f}` (`scattering_state`, the OUTGOING solution `coupling.v_dk_plus`
also builds), not `phi^-` (`scattering_state_minus`, which Eq. (37) has no
actual use for despite printing a `-` label). Both conjugations happen in
the same unscaled theory that Eq. (34)/(35) are stated in; `p. 012710-6`'s
"no conjugation" rule is the separate, later step of forming the
complex-scaled c-product of the (already de-conjugated) result.

Eq. (38) is printed with `phi_d*(r;R)`, conjugating the discrete state. That
is a no-op here: Eq. (34)/(35) already restrict this model to a real
`phi_d`, so `phi_d* = phi_d` and `j_dk` below conjugates nothing.

WEIGHTS. DVR vectors are coefficients (`c_j = f(x_j) sqrt(w_j)`), so a
c-product of two coefficient vectors already IS the quadrature integral and a
function multiplying them is elementwise. No explicit weights appear here.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid
from qscat.linalg import c_product

from .coupling import v_dk_plus
from .discrete_state import DiscreteState, electronic_hamiltonian
from .dissociation import solve_nuclear
from .ingredients import NrmIngredients, nrm_ingredients
from .nonlocal_potential import continue_to_tail, nonlocal_operator
from .scattering import incident_coefficients, scattering_state

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["j_dk", "nrm_ve_cross_section", "t_background", "t_resonant"]


def j_dk(
    elec_grid: FemDvrEcsGrid,
    phi_d: DiscreteState,
    R_values: npt.NDArray[np.float64],
    energy: float,
    ell: int,
) -> npt.NDArray[np.complex128]:
    """`J_dk(R)` of Eq. (38): the overlap of `phi_d(.;R)` with the incident wave.

    Parameters
    ----------
    elec_grid : FemDvrEcsGrid
        The electronic radial grid.
    phi_d : DiscreteState
        The discrete-state choice; consumed only through `phi_d.phi_d(R)`.
    R_values : ndarray
        Nuclear coordinates at which to evaluate the overlap.
    energy : float
        The real electron energy `E = k^2/2` (hartree), positive.
    ell : int
        Partial wave.

    Returns
    -------
    ndarray
        Complex `J_dk(R)`, shape `(R_values.size,)`.
    """
    if energy <= 0.0:
        raise ValueError(f"energy must be positive, got {energy}")
    k = float(np.sqrt(2.0 * energy))
    inc = incident_coefficients(elec_grid, k, ell)
    R = np.asarray(R_values, dtype=np.float64)
    out = np.empty(R.size, dtype=np.complex128)
    for j in range(R.size):
        out[j] = c_product(phi_d.phi_d(float(R[j])), inc)
    return out


def t_resonant(
    chi_f: npt.NDArray[np.complex128],
    v_dk_f: npt.NDArray[np.complex128],
    psi_d: npt.NDArray[np.complex128],
) -> complex:
    """`T^res` of Eq. (31): `<chi_vf| V^{-*}_{dk_f} |Psi_d^+>`.

    All three arguments live on the full nuclear grid. `chi_f` and `psi_d`
    are DVR coefficient vectors; `v_dk_f` is the coupling as a function of
    `R` (`nonlocal_potential.continue_to_tail`'s output), so it multiplies
    elementwise -- no explicit weights, per the module docstring.

    The paper's `V^{-*}_{dk_f}` is not conjugated here: Eq. (34) is the
    CONDITION (radial case, real discrete state, so `phi_k^- = (phi_k^+)^*`);
    Eq. (35) is why that condition drops the conjugation off the COUPLING
    FACTOR, collapsing it to `V^+_{dk_f}` by `H_el` Hermiticity ("we can use
    the matrix element `Vdk^+` but *without complex conjugation*", p.
    012710-4); p. 012710-6 is the separate reason the SCALAR PRODUCT itself
    (`qscat.linalg.c_product`) carries no conjugation, since `P H_el P` is
    complex symmetric, not Hermitian, under exterior complex scaling. The
    paper warns this distinction "becomes important when the background
    terms ... are added to the resonant T matrix, since `V^±_dk` are in
    general complex even when the discrete state is real" (p. 012710-4) --
    Eq. (37)'s background term must use the same non-conjugated convention.

    Parameters
    ----------
    chi_f : ndarray
        Final vibrational state, DVR coefficients on the nuclear grid.
    v_dk_f : ndarray
        `V^+_{dk_f}(R)`, the coupling continued to the tail, on the same
        nuclear grid.
    psi_d : ndarray
        `Psi_d^+(R)` of Eq. (52) (the same solution the DA path computes),
        DVR coefficients on the same nuclear grid.

    Returns
    -------
    complex
        `T^res_{vi->vf}`.
    """
    if not (chi_f.size == v_dk_f.size == psi_d.size):
        raise ValueError(
            "chi_f, v_dk_f and psi_d must have the same length, got "
            f"{chi_f.size}, {v_dk_f.size}, {psi_d.size}"
        )
    return complex(c_product(chi_f, v_dk_f * psi_d))


def _t_background_term2(
    chi_i: npt.NDArray[np.complex128],
    chi_f: npt.NDArray[np.complex128],
    v_dk_f: npt.NDArray[np.complex128],
    j_dk_i: npt.NDArray[np.complex128],
) -> complex:
    """The second term of Eq. (37): `<chi_vf| V_dk_f J_dk_i |chi_vi>`.

    Separated so it can be checked in isolation -- it is a pure contraction,
    while the first term needs a scattering solve at every nuclear node.
    """
    return complex(c_product(chi_f, v_dk_f * j_dk_i * chi_i))


def t_background(
    elec_grid: FemDvrEcsGrid,
    nuclear_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    R_values: npt.NDArray[np.float64],
    chi_i: npt.NDArray[np.complex128],
    chi_f: npt.NDArray[np.complex128],
    v_dk_f: npt.NDArray[np.complex128],
    e_kin_i: float,
    e_kin_f: float,
) -> complex:
    """`T^bg` of Eq. (37), the non-resonant background T-matrix.

    `T^bg = <chi_vf phi^-_{k_f}| V_int |chi_vi J_{k_i}> - <chi_vf| V_dk_f J_dk_i |chi_vi>`

    Eq. (37) is printed in the UNSCALED theory, where the bra `<phi^-_{k_f}|`
    genuinely means complex conjugation. Eq. (34) gives `(phi_k^-)^* =
    phi_k^+` for the radial case with a real discrete state, so conjugating
    that bra turns it into a (non-conjugating, c-product) pairing against
    `phi^+_{k_f}` -- the SAME reading `coupling.v_dk_plus` applies to the
    structurally identical Eq. (35), whose `V^{-*}_{dk} -> V^+_{dk}`
    collapse rests on the same `H_el` Hermiticity argument at the same real
    energy. `p. 012710-6`'s "no conjugation" rule is about the
    complex-scaled SCALAR PRODUCT itself and is a separate step from what
    Eq. (34) does to the bra's `-`/`+` label -- exactly as for `V_dk` in
    `t_resonant`'s docstring. So term 1's P-space state is built with
    `scattering_state` (OUTGOING boundary, i.e. `phi^+`) at the FINAL
    kinetic energy `e_kin_f`, c-producted (no further conjugation) against
    `V_int * J_{k_i}` -- not `scattering_state_minus`, which has no
    consumer in this module. The second term is a contraction; `v_dk_f` is
    expected on the FULL nuclear grid, already continued into the tail, and
    is the SAME non-conjugated `coupling.v_dk_plus` output `t_resonant`
    consumes, per the same Eq. (34)/(35) reasoning applied to `V^{-*}_{dk_f}`
    there. The paper warns the conjugation convention "becomes important
    when the background terms ... are added to the resonant T matrix, since
    `V^±_dk` are in general complex even when the discrete state is real"
    (p. 012710-4).

    PRA 77 shows this term is what a bare LCP curve omits -- non-negligible for
    elastic and even some inelastic channels, and largest for the broadest
    resonance (p. 012710-1, 012710-9-10).

    Parameters
    ----------
    elec_grid, nuclear_grid : FemDvrEcsGrid
        The electronic and nuclear radial grids.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice; consumed only through `phi_d.phi_d(R)`.
    R_values : ndarray
        Real nuclear nodes at which the electronic scattering solve of
        term 1 is evaluated -- must match `nuclear_grid`'s real nodes (the
        `continue_to_tail` contract).
    chi_i, chi_f : ndarray
        Initial/final vibrational states, DVR coefficients on the full
        nuclear grid.
    v_dk_f : ndarray
        `V^+_{dk_f}(R)` at the final energy, already continued to the full
        nuclear grid (`nonlocal_potential.continue_to_tail`'s output).
    e_kin_i, e_kin_f : float
        Initial/final electron kinetic energies (hartree); both positive.

    Returns
    -------
    complex
        `T^bg_{vi->vf}`.
    """
    if e_kin_i <= 0.0 or e_kin_f <= 0.0:
        raise ValueError(f"both energies must be positive, got {e_kin_i} and {e_kin_f}")
    R = np.asarray(R_values, dtype=np.float64)
    k_i = float(np.sqrt(2.0 * e_kin_i))
    inc_i = incident_coefficients(elec_grid, k_i, model.ell)
    ident = np.eye(elec_grid.n, dtype=np.complex128)

    inner = np.empty(R.size, dtype=np.complex128)
    for j in range(R.size):
        r_j = float(R[j])
        d = phi_d.phi_d(r_j)
        h_el = electronic_hamiltonian(elec_grid, model, r_j)
        p = ident - np.outer(d, d)  # Eq. (57)-(58), bilinear -- matches coupling.v_dk_plus
        php = p @ h_el @ p
        # phi^+ at the FINAL energy, per Eq. (34): conjugating Eq. (37)'s bra
        # <phi^-_{k_f}| turns it into a non-conjugating pairing with phi^+_{k_f}.
        phi_plus_f = scattering_state(php, elec_grid, e_kin_f, model.ell)
        v_int_j = np.asarray(model.v_int(elec_grid.points, r_j), dtype=np.complex128)
        inner[j] = c_product(phi_plus_f, v_int_j * inc_i)

    inner_full = continue_to_tail(inner, R, nuclear_grid)
    term1 = complex(c_product(chi_f, chi_i * inner_full))

    j_i = continue_to_tail(j_dk(elec_grid, phi_d, R, e_kin_i, model.ell), R, nuclear_grid)
    term2 = _t_background_term2(chi_i, chi_f, v_dk_f, j_i)
    return term1 - term2


def _psi_d_for_energy(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    e_kin: float,
    ing: NrmIngredients,
    n_states: int | None,
) -> npt.NDArray[np.complex128]:
    """`Psi_d^+` at one incident energy -- the Eq. (52) solve, shared with DA.

    Exposed (module-private) so `nrm_ve_cross_section` and
    `dissociation.nrm_da_cross_section` are provably solving the SAME
    equation rather than two independently-typed copies of it.

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice under test.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
    v_init : int
        Initial vibrational level.
    e_kin : float
        Incident electron kinetic energy (hartree), positive.
    ing : NrmIngredients
        Precomputed ingredients (`ingredients.nrm_ingredients`).
    n_states : int, optional
        Truncate the sum over projected electronic states in `F(E)`. `None`
        uses all.

    Returns
    -------
    ndarray
        `Psi_d^+(R)` as DVR coefficients on `nuclear_grid`.
    """
    e_total = float(e_kin) + float(eps[v_init])
    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuclear_grid)
    v_dk_i = continue_to_tail(
        v_dk_plus(elec_grid, model, phi_d, ing.R, float(e_kin)), ing.R, nuclear_grid
    )
    f = nonlocal_operator(ing, nuclear_grid, model, e_total, n_states=n_states)
    return solve_nuclear(nuclear_grid, model.mu, v_d_full, f, v_dk_i * chi[v_init], e_total)


def nrm_ve_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    elec_grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: float | npt.ArrayLike,
    *,
    ingredients: NrmIngredients | None = None,
    n_states: int | None = None,
    include_background: bool = True,
) -> npt.NDArray[np.float64]:
    """`sigma_{v_init->v'}(E)` in the nonlocal resonance model (bohr^2).

    `sigma = 4 pi^3 |T^res + T^bg|^2 / k_i^2`, Eq. (28)/(31)/(37), on
    `qscat.core.driven`'s own normalization so the exact 2-D and nonlocal
    curves compare directly rather than through two conventions that happen
    to agree.

    `T^res` (Eq. 31) is a contraction of the SAME `Psi_d^+` solution
    `dissociation.nrm_da_cross_section` computes (`_psi_d_for_energy`, built
    once per energy and reused across every `vprimes` entry). `T^bg`
    (Eq. 37) is the non-resonant background PRA 77 shows a bare LCP curve is
    missing -- largest for the broadest resonance.

    `include_background=False` drops `T^bg`, giving the paper's "nonlocal"
    curve as against its "nonlocal + bg" one (Figs. 4 and 8 plot both). That
    difference is the method's own argument, not a debugging switch.

    Parameters
    ----------
    nuclear_grid, elec_grid : FemDvrEcsGrid
        The nuclear and electronic radial grids.
    model : ResonanceModel
        The molecule.
    phi_d : DiscreteState
        The discrete-state choice under test.
    eps, chi : ndarray
        Neutral vibrational energies and states (`qscat.core.vibrational`).
    v_init : int
        Initial vibrational level.
    vprimes : list of int
        Final vibrational levels.
    E : float or array
        Incident electron kinetic energy or energies (hartree).
    ingredients : NrmIngredients, optional
        Precomputed ingredients; built here if omitted. Pass them in when
        sweeping energies or comparing discrete-state choices -- they are
        energy-independent and dominate the cost.
    n_states : int, optional
        Truncate the sum over projected electronic states in `F(E)`. `None`
        uses all.
    include_background : bool, default True
        Add `T^bg` (Eq. 37) to `T^res` before squaring. `True` is the
        paper's "nonlocal + background" curve; `False` is its bare
        "nonlocal" curve.

    Returns
    -------
    ndarray
        `sigma_{v_init->v'}` per energy; scalar `E` returns shape
        `(len(vprimes),)`, array `E` returns `(len(E), len(vprimes))` --
        `driven.ve_cross_section`'s convention. A closed channel
        (`E_tot - eps_vf <= 0`, or `E <= 0`) contributes `0.0`.
    """
    real = nuclear_grid.points.imag == 0.0
    R_desc = np.sort(nuclear_grid.points[real].real)[::-1]
    ing = ingredients or nrm_ingredients(elec_grid, model, phi_d, R_desc)

    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros((e_arr.size, len(vprimes)), dtype=np.float64)
    for ie, e_kin in enumerate(e_arr):
        if float(e_kin) <= 0.0:
            continue
        e_total = float(e_kin) + float(eps[v_init])
        psi_d = _psi_d_for_energy(
            nuclear_grid, elec_grid, model, phi_d, eps, chi, v_init, float(e_kin), ing, n_states
        )
        for jv, vp in enumerate(vprimes):
            excess = e_total - float(eps[vp])
            if excess <= 0.0:
                continue  # closed channel
            v_dk_f = continue_to_tail(
                v_dk_plus(elec_grid, model, phi_d, ing.R, excess), ing.R, nuclear_grid
            )
            t = t_resonant(chi[vp], v_dk_f, psi_d)
            if include_background:
                t += t_background(
                    elec_grid,
                    nuclear_grid,
                    model,
                    phi_d,
                    ing.R,
                    chi[v_init],
                    chi[vp],
                    v_dk_f,
                    float(e_kin),
                    excess,
                )
            out[ie, jv] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e_kin))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
