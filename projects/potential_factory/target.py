"""What a molecule IS to the factory: tiered, provenance-carrying target data."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from scipy.interpolate import CubicSpline, RegularGridInterpolator

__all__ = ["Curve", "NeutralTarget", "ResonanceTarget", "CouplingTarget", "Provenance", "Target"]

FArr = npt.NDArray[np.float64]


class Curve:
    """A real curve on one coordinate: a cubic spline through a table, or a callable."""

    def __init__(self, fn: Callable[[FArr], FArr], x: FArr | None) -> None:
        self._fn = fn
        self.x = x

    @classmethod
    def from_table(cls, x: npt.ArrayLike, y: npt.ArrayLike) -> Curve:
        xs = np.asarray(x, dtype=np.float64)
        ys = np.asarray(y, dtype=np.float64)
        order = np.argsort(xs)
        spl = CubicSpline(xs[order], ys[order], extrapolate=False)
        return cls(lambda t: np.asarray(spl(t), dtype=np.float64), xs[order])

    @classmethod
    def from_callable(cls, fn: Callable[[FArr], FArr]) -> Curve:
        return cls(
            lambda t: np.asarray(fn(np.asarray(t, dtype=np.float64)), dtype=np.float64), None
        )

    def __call__(self, x: npt.ArrayLike) -> FArr:
        return self._fn(np.asarray(x, dtype=np.float64))


@dataclass(frozen=True)
class NeutralTarget:
    curve: Curve | None
    constants: dict[str, float]
    R_range: tuple[float, float]


@dataclass(frozen=True)
class ResonanceTarget:
    v_ion: Curve
    gamma: Curve
    ea: float
    R_range: tuple[float, float]


@dataclass(frozen=True)
class CouplingTarget:
    gamma_tilde: Callable[[npt.ArrayLike, npt.ArrayLike], FArr]
    eps_window: tuple[float, float]
    R_range: tuple[float, float]
    alpha_exponent: float
    # The table's own (ascending) axes, when built `from_table` -- `None` for
    # an analytic form (`from_alt_houfek`) that has no nodes to evaluate
    # exactly on. Interpolating log(Gamma_tilde) linearly in log(eps) is
    # exact ON these nodes (the interpolant passes through its own data) but
    # has real curvature error BETWEEN them, since `Gamma_tilde ~ eps^a
    # exp(-B(R) eps)`'s `exp(-B(R) eps)` factor is not linear in log(eps) --
    # a fitter that queries `gamma_tilde` off-node inherits that error as
    # spurious residual. Querying ON these nodes avoids it entirely.
    eps_nodes: FArr | None = None
    R_nodes: FArr | None = None

    @classmethod
    def from_alt_houfek(
        cls,
        *,
        a0: float,
        a1: float,
        a2: float,
        b0: float,
        b1: float,
        alpha: float,
        R_range: tuple[float, float],
        eps_window: tuple[float, float] = (0.002, 0.22),
    ) -> CouplingTarget:
        """Alt & Houfek, PRA 103, 032829 (2021) Eq. (25)-(27): 2pi eps^a A(R) exp(-B(R) eps)."""

        def g(eps: npt.ArrayLike, R: npt.ArrayLike) -> FArr:
            e = np.asarray(eps, dtype=np.float64)
            r = np.asarray(R, dtype=np.float64)
            A = (a0 + a1 * r) * np.exp(a2 * r)
            B = b0 + b1 * r
            return np.asarray(2.0 * np.pi * e**alpha * A * np.exp(-B * e), dtype=np.float64)

        return cls(g, eps_window, R_range, alpha)

    @classmethod
    def from_table(
        cls,
        eps: npt.ArrayLike,
        R: npt.ArrayLike,
        table: npt.ArrayLike,
        *,
        alpha: float,
    ) -> CouplingTarget:
        """Interpolate in `(log eps, R)` on `log(Gamma_tilde)`.

        A plain bilinear interpolant on `(eps, R)` badly distorts the
        near-threshold Wigner power law `Gamma_tilde ~ eps^(l+1/2)` when the
        eps table is coarse (linear-in-eps interpolation across a bin where
        the true value changes by orders of magnitude), as measured.
        Interpolating the LOG of the table on a
        LOG-eps axis makes the interpolant piecewise power-law in eps
        instead, which tracks the threshold exponent correctly between
        table nodes. `table` is clamped to `>= 1e-300` before taking the log
        so an exact zero doesn't diverge; NaN outside the table is
        preserved (exp(NaN) is NaN).
        """
        e = np.asarray(eps, dtype=np.float64)
        r = np.asarray(R, dtype=np.float64)
        tab = np.asarray(table, dtype=np.float64)  # shape (e.size, r.size)
        log_tab = np.log(np.clip(tab, 1e-300, None))
        interp = RegularGridInterpolator(
            (np.log(e), r), log_tab, bounds_error=False, fill_value=np.nan
        )

        def g(eps_q: npt.ArrayLike, R_q: npt.ArrayLike) -> FArr:
            eq = np.asarray(eps_q, dtype=np.float64)
            rq = np.asarray(R_q, dtype=np.float64)
            eb, rb = np.broadcast_arrays(eq, rq)
            pts = np.stack([np.log(eb).ravel(), rb.ravel()], axis=-1)
            out = np.exp(interp(pts)).reshape(eb.shape)
            return np.asarray(out, dtype=np.float64)

        return cls(
            g,
            (float(e.min()), float(e.max())),
            (float(r.min()), float(r.max())),
            alpha,
            eps_nodes=e,
            R_nodes=r,
        )


@dataclass(frozen=True)
class Provenance:
    source: str
    locator: str


@dataclass(frozen=True)
class Target:
    name: str
    mu: float
    ell: int
    charge: int
    coordinates: tuple[str, ...]
    neutral: NeutralTarget | None
    resonance: ResonanceTarget | None
    coupling: CouplingTarget | None
    eigenphase: object | None = None  # reserved: T2 tables (loader only, no loss in v1)
    provenance: dict[str, Provenance] = field(default_factory=dict)
