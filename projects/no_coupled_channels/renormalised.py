"""The two-centre well with `lam` renormalised to preserve the anion curve.

The bare `TwoCentreWell` splits `lam` between the two centres without
rescaling it, so the deeper well keeps only `(1+kappa)/2` of the original
depth and the anion unbinds beyond `R ~ 0.7/s` for ANY `s > 0`. The
dissociation limit is then wrong, and -- because this model sits within 5 % of
not binding at all -- the resonance position is badly wrong long before that,
by hundreds of meV right through the Franck-Condon region.

`RenormalisedTwoCentreWell` multiplies `lam(R)` by a per-`R` factor `f(R)`
solved so the two-centre model reproduces the SHIPPED model's `E_res(R)`. The
curve, the crossing and the asymptote are then pinned by construction, and
`Gamma(R)` is left free -- which is the physics under study: the anisotropy
should change how the electron COUPLES, not where the anion state sits.

`f` is supplied as tabulated `(R, f)` pairs from
`validation.coupled.renormalise`, not recomputed here, because solving it
needs the resonance machinery and `projects/` must not import `validation/`.
Outside the tabulated range `f` is held at its end values: below the table
the well separation `d = sR/2` is small and `f -> 1`; above it, `f` has
already reached the analytic `2/(1+kappa)`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from qscat.model import DiatomicResonanceModel

from projects.no_coupled_channels.anisotropy import TwoCentreWell

__all__ = ["RenormalisedTwoCentreWell", "ScaleTable"]


@dataclass(frozen=True)
class ScaleTable:
    """Tabulated `f(R)`, interpolated linearly and clamped outside the table.

    Clamped rather than extrapolated on purpose: `f` is a bounded, monotone
    quantity between 1 and `2/(1+kappa)`, and a linear extrapolation off
    either end leaves that range immediately -- which would put the well at a
    depth no measurement supports, in exactly the region (large `R`) where
    nothing checks it.

    **The table MUST end at or before the nuclear ECS radius.** Interpolation
    is on `Re R`, which is only legitimate where `f` is CONSTANT: past the ECS
    radius the grid points are complex, and a factor that varies with `Re R`
    there is not an analytic function of `R`, so it destroys the analytic
    continuation the whole exterior-complex-scaling construction rests on. Held
    constant instead, `f` is trivially analytic and the tail is undisturbed.
    `for_ecs_grid` builds a table with that property; the invariant is checked
    in `__post_init__` only in the sense that a caller passing a longer table
    is passing one whose tail values will never be reached by a real `R`.
    """

    R: npt.NDArray[np.float64]
    f: npt.NDArray[np.float64]

    def __call__(self, R: npt.ArrayLike) -> npt.NDArray[np.float64]:
        r = np.real(np.asarray(R, dtype=np.complex128))
        return np.interp(r, self.R, self.f, left=self.f[0], right=self.f[-1])

    @classmethod
    def for_ecs_grid(
        cls, R: npt.ArrayLike, f: npt.ArrayLike, *, r_ecs: float, asymptote: float
    ) -> ScaleTable:
        """Truncate a solved `f(R)` at the ECS radius, ending on `asymptote`.

        Keeps every solved point strictly inside the real region, then brings
        the last interval onto the analytic large-separation value so the
        constant held over the whole rotated tail is the RIGHT constant rather
        than wherever the last solved point happened to land.
        """
        Ra = np.asarray(R, dtype=float)
        fa = np.asarray(f, dtype=float)
        keep = Ra < r_ecs
        return cls(
            R=np.append(Ra[keep], r_ecs),
            f=np.append(fa[keep], asymptote),
        )


@dataclass(frozen=True)
class _ScaledBase:
    """`base` with `lam` scaled by a table; every other attribute delegated."""

    base: DiatomicResonanceModel
    table: ScaleTable

    def __getattr__(self, name: str) -> object:
        return getattr(self.base, name)

    def lam(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        return np.asarray(self.base.lam(R)) * self.table(R)


@dataclass(frozen=True)
class RenormalisedTwoCentreWell:
    """`TwoCentreWell` whose `lam` is scaled to pin the shipped anion curve.

    Composes rather than subclasses `TwoCentreWell`, so the angular quadrature,
    the `s = 0` embedding and the Legendre selection rules are inherited
    unchanged and cannot drift: at `s = 0` the table is all ones and this IS
    the shipped model, structurally.
    """

    base: DiatomicResonanceModel
    table: ScaleTable
    s: float = 0.0
    kappa: float = 0.0
    Lambda: int = 1
    n_nodes: int = 64
    _well: TwoCentreWell = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "_well",
            TwoCentreWell(
                base=_ScaledBase(self.base, self.table),
                s=self.s,
                kappa=self.kappa,
                Lambda=self.Lambda,
                n_nodes=self.n_nodes,
            ),
        )

    def __getattr__(self, name: str) -> object:
        return getattr(self._well, name)
