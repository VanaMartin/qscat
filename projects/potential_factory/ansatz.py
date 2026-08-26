"""The analytic, ECS-safe surface family the factory emits.

`FlexibleDiatomicModel` is a strict superset of `qscat.model.DiatomicResonanceModel`:
an Expanded-Morse-Oscillator neutral curve (Le Roy; one beta == Morse), a Gaussian
electron-molecule well whose depth `lam(R)` AND range `alpha(R)` are smooth
functions of `R`, and an optional repulsive Gaussian shell `shell(R) exp(-alpha_b
(r - r_b)^2)`. `from_diatomic` embeds the published models exactly, so N2/NO/F2
are points of this parameter space (the round-trip oracle relies on that).

Every term is entire in `r`; in `R` the only singularities are the poles of
`y_p` at `|R| = R_e`, which an ECS tail pivoted at `R_0 > R_e` never reaches.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd
from qscat.model import DiatomicResonanceModel

__all__ = [
    "FlexibleDiatomicModel",
    "SmoothR",
    "TailR",
    "from_diatomic",
    "pack",
    "params",
    "unpack",
    "with_params",
    "y_p",
]


def y_p(R: npt.ArrayLike, R_e: float, p: int) -> npt.NDArray[np.complex128]:
    """Le Roy's dimensionless radial variable `(R^p - R_e^p) / (R^p + R_e^p)`."""
    Rc = np.asarray(R, dtype=np.complex128)
    out = (Rc**p - R_e**p) / (Rc**p + R_e**p)
    return np.asarray(out, dtype=np.complex128)


@dataclass(frozen=True)
class SmoothR:
    """`f(R) = f_inf + f_0 / (1 + exp(f_1 (R - R_f))) * (1 + sum_i coeffs[i] y_p(R)^(i+1))`.

    With `coeffs == ()` this is exactly Houfek's sigmoid
    `lambda_inf + lambda_0 / (1 + exp(lambda_1 (R - R_lambda)))`.
    """

    f_inf: float
    f_0: float
    f_1: float
    R_f: float
    coeffs: tuple[float, ...] = ()
    R_e: float = 2.0
    p: int = 3

    def __call__(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        sig = self.f_0 / (1.0 + np.exp(self.f_1 * (Rc - self.R_f)))
        poly = np.ones_like(Rc)
        if self.coeffs:
            y = y_p(Rc, self.R_e, self.p)
            poly = poly + sum(c * y ** (i + 1) for i, c in enumerate(self.coeffs))
        return np.asarray(self.f_inf + sig * poly, dtype=np.complex128)


@dataclass(frozen=True)
class TailR:
    """`f(R) = f_inf + (1 - y_q(R)) * sum_k coeffs[k] y_p(R)^k`.

    The long-range-correct alternative to `SmoothR`: every term carries the
    factor `1 - y_q(R) ~ 2 (R_e/R)^q` as `R -> inf`, so `f(inf) == f_inf`
    EXACTLY and the approach is the power law `R^-q` -- with `q = 4` the
    ion--atom polarisation form the anion curve must follow (the
    `polarisation_tail` of `target.py`). Inside, `P(y_p)` is a plain
    polynomial in Le Roy's bounded variable `y_p in (-1, 1)`, so the fit is
    linear in `coeffs` and well conditioned -- no sigmoid inflection to run
    off to large `R`, which is how the `SmoothR` form held `-EA` at one
    node and missed it by 0.2 eV beyond (measured on O2). `SmoothR` stays
    for the published models, whose `lam(R)` IS a sigmoid.
    """

    f_inf: float
    coeffs: tuple[float, ...]
    R_e: float = 2.0
    p: int = 3
    q: int = 4

    def __call__(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        poly = np.zeros_like(Rc)
        if self.coeffs:
            y = y_p(Rc, self.R_e, self.p)
            poly = poly + sum(c * y**k for k, c in enumerate(self.coeffs))
        tail = 1.0 - y_p(Rc, self.R_e, self.q)
        return np.asarray(self.f_inf + tail * poly, dtype=np.complex128)


@dataclass(frozen=True)
class FlexibleDiatomicModel:
    """EMO neutral curve + Gaussian well with `lam(R)`, `alpha(R)` + optional shell."""

    mu: float
    ell: int
    D_e: float
    R_e: float
    betas: tuple[float, ...]
    p: int
    lam: SmoothR | TailR
    alpha: SmoothR | TailR
    shell: SmoothR | None
    alpha_b: float
    r_b: float
    charge: int = 0

    # -- neutral curve -------------------------------------------------------
    def beta_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        y = y_p(Rc, self.R_e, self.p)
        out = sum(b * y**i for i, b in enumerate(self.betas))
        return np.asarray(out + 0.0 * Rc, dtype=np.complex128)

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """EMO: `D_e [ (1 - exp(-beta(R)(R - R_e)))^2 - 1 ]`, minimum -D_e at R_e, 0 at infinity."""
        Rc = np.asarray(R, dtype=np.complex128)
        e = np.exp(-self.beta_R(Rc) * (Rc - self.R_e))
        return np.asarray(self.D_e * ((1.0 - e) ** 2 - 1.0), dtype=np.complex128)

    # -- interaction ---------------------------------------------------------
    def lam_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        return self.lam(R)

    def alpha_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        return self.alpha(R)

    def shell_R(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        Rc = np.asarray(R, dtype=np.complex128)
        if self.shell is None:
            return np.zeros_like(Rc)
        return self.shell(Rc)

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        Rc = np.asarray(R, dtype=np.complex128)
        well = -self.lam_R(Rc) * np.exp(-self.alpha_R(Rc) * rr**2)
        if self.shell is None:
            return np.asarray(well, dtype=np.complex128)
        barrier = self.shell_R(Rc) * np.exp(-self.alpha_b * (rr - self.r_b) ** 2)
        return np.asarray(well + barrier, dtype=np.complex128)

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        rr = np.asarray(r, dtype=np.complex128)
        out = self.v0(R) + self.ell * (self.ell + 1) / (2.0 * rr**2) + self.v_int(rr, R)
        return np.asarray(out, dtype=np.complex128)

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        return potential_nd(tgrid, self.v_int)

    def with_shell(
        self, shell: SmoothR | None, alpha_b: float, r_b: float
    ) -> FlexibleDiatomicModel:
        return replace(self, shell=shell, alpha_b=alpha_b, r_b=r_b)


def from_diatomic(model: DiatomicResonanceModel) -> FlexibleDiatomicModel:
    """Embed a published `DiatomicResonanceModel` exactly (to round-off)."""
    lam0 = (model.lambda_c - model.lambda_inf) * (
        1.0 + np.exp(model.lambda_1 * (model.R_c - model.R_lambda))
    )
    return FlexibleDiatomicModel(
        mu=model.mu,
        ell=model.ell,
        D_e=model.D0,
        R_e=model.R0,
        betas=(model.alpha0,),
        p=3,
        lam=SmoothR(
            f_inf=model.lambda_inf,
            f_0=float(lam0),
            f_1=model.lambda_1,
            R_f=model.R_lambda,
            R_e=model.R0,
        ),
        alpha=SmoothR(f_inf=model.alpha_c, f_0=0.0, f_1=1.0, R_f=0.0, R_e=model.R0),
        shell=None,
        alpha_b=1.0,
        r_b=0.0,
        charge=model.charge,
    )


def _smooth_params(prefix: str, s: SmoothR | TailR) -> dict[str, float]:
    out = {f"{prefix}.f_inf": s.f_inf}
    if isinstance(s, SmoothR):
        out[f"{prefix}.f_0"] = s.f_0
        out[f"{prefix}.f_1"] = s.f_1
        out[f"{prefix}.R_f"] = s.R_f
    for i, c in enumerate(s.coeffs):
        out[f"{prefix}.c{i}"] = c
    return out


def _smooth_update(s: SmoothR | TailR, prefix: str, upd: Mapping[str, float]) -> SmoothR | TailR:
    kw = (
        {"f_inf": s.f_inf, "f_0": s.f_0, "f_1": s.f_1, "R_f": s.R_f}
        if isinstance(s, SmoothR)
        else {"f_inf": s.f_inf}
    )
    coeffs = list(s.coeffs)
    for key, val in upd.items():
        if not key.startswith(prefix + "."):
            continue
        field = key[len(prefix) + 1 :]
        if field in kw:
            kw[field] = float(val)
        elif field.startswith("c") and field[1:].isdigit() and int(field[1:]) < len(coeffs):
            coeffs[int(field[1:])] = float(val)
        else:
            raise KeyError(key)
    return replace(s, coeffs=tuple(coeffs), **kw)


def params(model: FlexibleDiatomicModel) -> dict[str, float]:
    """Flat, ordered name -> value map of every FITTABLE parameter."""
    out: dict[str, float] = {"D_e": model.D_e, "R_e": model.R_e}
    for i, b in enumerate(model.betas):
        out[f"beta{i}"] = b
    out.update(_smooth_params("lam", model.lam))
    out.update(_smooth_params("alpha", model.alpha))
    if model.shell is not None:
        out.update(_smooth_params("shell", model.shell))
        out["alpha_b"] = model.alpha_b
        out["r_b"] = model.r_b
    return out


def with_params(
    model: FlexibleDiatomicModel, updates: Mapping[str, float]
) -> FlexibleDiatomicModel:
    """A copy with the named parameters replaced. Unknown names raise KeyError."""
    known = params(model)
    for key in updates:
        if key not in known:
            raise KeyError(key)
    betas = list(model.betas)
    for i in range(len(betas)):
        if f"beta{i}" in updates:
            betas[i] = float(updates[f"beta{i}"])
    return replace(
        model,
        D_e=float(updates.get("D_e", model.D_e)),
        R_e=float(updates.get("R_e", model.R_e)),
        betas=tuple(betas),
        lam=_smooth_update(model.lam, "lam", updates),
        alpha=_smooth_update(model.alpha, "alpha", updates),
        shell=None if model.shell is None else _smooth_update(model.shell, "shell", updates),
        alpha_b=float(updates.get("alpha_b", model.alpha_b)),
        r_b=float(updates.get("r_b", model.r_b)),
    )


def pack(model: FlexibleDiatomicModel, names: Sequence[str]) -> npt.NDArray[np.float64]:
    p = params(model)
    return np.array([p[n] for n in names], dtype=np.float64)


def unpack(
    model: FlexibleDiatomicModel, names: Sequence[str], x: npt.ArrayLike
) -> FlexibleDiatomicModel:
    xs = np.asarray(x, dtype=np.float64)
    return with_params(model, dict(zip(names, xs.tolist(), strict=True)))
