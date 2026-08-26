"""The H2+ ionic resonance-model form.

`IonicResonanceModel` implements the narrowed `ResonanceModel` protocol for a
CHARGED target (H2+, `charge=-1`): a Morse ion-core curve `v0(R)`, a
sigma-capture electron-molecule interaction `v_int(r, R)` (a `tanh`-gated,
Gaussian-in-r form distinct from `DiatomicResonanceModel`'s sigmoid+Gaussian
form), and a `surface` that adds the `charge/r` electron-core Coulomb tail on
top of the neutral-diatomic terms. Ported (formula only) from the extracted
H2+ dissociative-recombination model -- see `docs/physics/h2plus-dr.md`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd

__all__ = ["IonicResonanceModel"]


@dataclass(frozen=True)
class IonicResonanceModel:
    """The H2+ Morse + sigma-capture + Coulomb-tail resonance-model form,
    implementing `ResonanceModel`. Carries no grid state -- purely the
    potential-surface parameters; `qscat.model.library` holds the H2+
    instance (`H2P`).
    """

    mu: float
    ell: int
    charge: int
    V0: float
    R0: float
    alpha: float
    a1: float
    a2: float
    a3: float
    a4: float

    # Hvizdos et al., Phys. Rev. A 97, 022704 (2018), Sec. II: the nuclear ECS
    # angle must stay below pi/8 or the quartic `a3 * R**4` term in `v_int`
    # diverges under the rotation (4*theta < pi/2). The electronic bound is
    # pi/4 (from exp(-r^2/3), 2*theta < pi/2). Neutral diatomics have no such
    # bound -- their Morse + Gaussian forms are entire -- so the `ResonanceModel`
    # protocol does NOT declare this attribute at all; it is an ionic-only
    # extension that consumers read with `getattr(model, ..., None)`
    # (`qscat.core.lcp.levels._check_angle_bound`), treating its absence as "no bound".
    max_nuclear_ecs_angle_deg: float = 22.5

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Ion-core Morse potential (Hartree). Minimum -V0 at R0.

        `R` is converted to `complex128`, never a real dtype: it may be
        complex (ECS-rotated tail points), and coercing to real would
        silently discard Im(R) and corrupt the analytic continuation the
        exterior-complex-scaling method relies on.
        """
        Rc = np.asarray(R, dtype=np.complex128)
        a = self.alpha
        out = self.V0 * (np.exp(-2 * a * (Rc - self.R0)) - 2 * np.exp(-a * (Rc - self.R0)))
        return np.asarray(out, dtype=np.complex128)

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Sigma-capture electron-molecule interaction potential (Hartree).

        `v_int(r,R) = -a1 * (1 - tanh Q(R)) * S(R) * (exp(-r^2/3) / r)`,
        `Q(R) = (a2 - R - a3 R^4) / 7`, `S(R) = tanh(R/a4)^4`.

        `r`/`R` are converted to `complex128`, never a real dtype: they may
        be complex (ECS-rotated tail points), and coercing to real would
        silently discard Im(r)/Im(R) and corrupt the analytic continuation
        the exterior-complex-scaling method relies on.
        """
        rr = np.asarray(r, dtype=np.complex128)
        Rc = np.asarray(R, dtype=np.complex128)
        Q = (self.a2 - Rc - self.a3 * Rc**4) / 7.0
        S = np.tanh(Rc / self.a4) ** 4
        out = -self.a1 * (1 - np.tanh(Q)) * S * (np.exp(-(rr**2) / 3.0) / rr)
        return np.asarray(out, dtype=np.complex128)

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """The full surface `v0(R) + v_int(r, R) + ell(ell+1)/(2 r^2) +
        charge/r` -- the `charge/r` term is the electron-core Coulomb tail
        (`charge=-1` for H2+, the `-1/r` attraction).

        Must not coerce to a real dtype: `r`/`R` are complex on the ECS
        tails.
        """
        rr = np.asarray(r, dtype=np.complex128)
        out = (
            self.v0(R)
            + self.v_int(rr, R)
            + self.ell * (self.ell + 1) / (2.0 * rr**2)
            + self.charge / rr
        )
        return np.asarray(out, dtype=np.complex128)

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """`H_2D` on `tgrid` (axis 0 = electronic r, axis 1 = nuclear R)."""
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        """`V_int` evaluated on the tensor grid, flattened (C order)."""
        return potential_nd(tgrid, self.v_int)
