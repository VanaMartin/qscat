"""The two-centre anisotropic electron-molecule well, and the coupled
channel potential it generates.

    V_int(r, R) = -(1+kappa)/2 lam(R) exp(-alpha |r_vec - d zhat|^2)
                  -(1-kappa)/2 lam(R) exp(-alpha |r_vec + d zhat|^2),
    d = s R / 2.

Two knobs, both geometric. `s` moves the wells from the molecular centre
(s = 0, where the sum collapses to the shipped isotropic
-lam(R) exp(-alpha r^2) for ANY kappa) out onto the nuclei (s = 1). `kappa`
is the amplitude asymmetry, and it IS homonuclear-versus-heteronuclear: at
kappa = 0 the well is symmetric, only even Legendre components survive, and
within Lambda = 1 the l = 1 resonance can reach only l = 3. Turning kappa on
opens l = 1 <-> l = 2, the coupling a homonuclear molecule forbids.

The channel matrix is one angular quadrature,

    V_{ll'}(r, R) = int_{-1}^{1} Theta_{l,Lambda}(x) V_int(r, x, R)
                    Theta_{l',Lambda}(x) dx,

evaluated by Gauss-Legendre in x = cos(theta). That route rather than the
closed-form Legendre expansion because it is complex-safe for free (the nodes
are real; r and R carry the ECS phase), because at s = 0 it returns the
Kronecker delta to round-off rather than to a tolerance, and because it
survives a change of well shape. The closed form is kept as
`v_lambda_closed_form`, the differential oracle for the tests.

See docs/physics/coupled-partial-waves.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.model import DiatomicResonanceModel
from scipy.special import eval_legendre, ive

from projects.no_coupled_channels.angular import theta_lm

__all__ = ["TwoCentreWell", "v_lambda_closed_form"]


@dataclass(frozen=True)
class TwoCentreWell:
    """The shipped model's Gaussian well, moved onto the two nuclei.

    `base` supplies lam(R), alpha_c, v0(R), mu and the charge; nothing is
    restated here, so the s = 0 embedding is structural and cannot drift.
    """

    base: DiatomicResonanceModel
    s: float = 0.0
    kappa: float = 0.0
    Lambda: int = 1
    n_nodes: int = 64

    def offset(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Well offset `d = s R / 2` from the molecular centre."""
        return np.asarray(0.5 * self.s * np.asarray(R, dtype=np.complex128))

    def v_int_angular(
        self, r: npt.ArrayLike, x: float, R: npt.ArrayLike
    ) -> npt.NDArray[np.complex128]:
        """`V_int` at a single real `x = cos(theta)`.

        The Gaussian is evaluated at the SHIFTED argument in one expression.
        Splitting it into exp(-a(r^2+d^2)) * exp(2 a r d x) overflows: the
        second factor reaches e^96 against a first factor of e^-265 on NO's
        electronic grid.
        """
        rr = np.asarray(r, dtype=np.complex128)
        d = self.offset(R)
        a = self.base.alpha_c
        lam = self.base.lam(R)
        rho_a = rr**2 + d**2 - 2.0 * rr * d * x
        rho_b = rr**2 + d**2 + 2.0 * rr * d * x
        out = (
            -0.5
            * lam
            * ((1.0 + self.kappa) * np.exp(-a * rho_a) + (1.0 - self.kappa) * np.exp(-a * rho_b))
        )
        return np.asarray(out, dtype=np.complex128)

    def _quadrature(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        return np.polynomial.legendre.leggauss(self.n_nodes)

    def _project(
        self, coeff: npt.NDArray[np.float64], r: npt.ArrayLike, R: npt.ArrayLike
    ) -> npt.NDArray[np.complex128]:
        """Sum `coeff[i] * V_int(r, x_i, R)` over the quadrature nodes."""
        x, _ = self._quadrature()
        acc: npt.NDArray[np.complex128] | None = None
        for c, xi in zip(coeff, x, strict=True):
            term = c * self.v_int_angular(r, float(xi), R)
            acc = term if acc is None else acc + term
        assert acc is not None  # n_nodes >= 1
        return np.asarray(acc, dtype=np.complex128)

    def v_block(
        self, l: int, lp: int, r: npt.ArrayLike, R: npt.ArrayLike
    ) -> npt.NDArray[np.complex128]:
        """The coupled channel potential `V_{ll'}(r, R)` for `Lambda`."""
        x, w = self._quadrature()
        coeff = w * theta_lm(l, self.Lambda, x) * theta_lm(lp, self.Lambda, x)
        return self._project(coeff, r, R)

    def v_lambda(self, lam: int, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """The Legendre component `v_lambda(r, R)` -- a DIAGNOSTIC, not the
        production path: `v_block` integrates the potential directly."""
        x, w = self._quadrature()
        coeff = (2 * lam + 1) / 2 * w * eval_legendre(lam, x)
        return self._project(np.asarray(coeff, dtype=np.float64), r, R)


def v_lambda_closed_form(
    well: TwoCentreWell, lam: int, r: npt.ArrayLike, R: npt.ArrayLike
) -> npt.NDArray[np.complex128]:
    """`v_lambda` in closed form -- the differential oracle, REAL argument only.

    Expanding exp(z cos theta) = sum_lambda (2 lambda + 1) i_lambda(z)
    P_lambda(cos theta) with z = 2 alpha r d gives

        v_lambda = -(2 lambda + 1) [lam_A + (-1)^lambda lam_B]
                   exp(-alpha (r - d)^2) itilde_lambda(z),

    where `itilde_lambda(z) = exp(-z) i_lambda(z) = sqrt(pi/2z) ive(lambda+1/2, z)`
    is the EXPONENTIALLY SCALED modified spherical Bessel function. The scaling
    is not optional: the unscaled identity carries exp(-alpha(r^2+d^2))
    i_lambda(z), whose two factors overflow and underflow separately.

    The bracket makes the kappa = 0 symmetry manifest -- there lam_A = lam_B and
    every odd lambda vanishes by inspection.

    `z = 0` (at `r = 0` or `s = 0`) is the isotropic limit and is returned
    directly rather than through the singular `sqrt(pi/2z)`.
    """
    a = well.base.alpha_c
    rr = np.asarray(r, dtype=np.float64)
    d = np.real(well.offset(R))
    lam_r = np.real(well.base.lam(R))
    lam_a = 0.5 * (1.0 + well.kappa) * lam_r
    lam_b = 0.5 * (1.0 - well.kappa) * lam_r
    z = 2.0 * a * rr * d
    bracket = lam_a + (-1) ** lam * lam_b

    isotropic = -bracket * np.exp(-a * rr**2) if lam == 0 else np.zeros_like(z)
    safe = np.where(z > 0.0, z, 1.0)
    scaled = np.sqrt(np.pi / (2.0 * safe)) * ive(lam + 0.5, safe)
    general = -(2 * lam + 1) * bracket * np.exp(-a * (rr - d) ** 2) * scaled
    return np.asarray(np.where(z > 0.0, general, isotropic), dtype=np.complex128)
