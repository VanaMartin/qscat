"""The `ResonanceModel` protocol and the shared diatomic resonance-model form.

`ResonanceModel` is the ENTIRE structural contract `qscat.core`'s
model-independent solvers depend on -- a driven-equation or time-dependent
solver type-annotates against this protocol, never against the concrete
`DiatomicResonanceModel` dataclass below (or any other model type that may
join later, e.g. the parked angular coupled-channel model -- see
`docs/physics/angular-coupled-channels.md`). `qscat.core` must never import
this module's concrete class, only the protocol.

`DiatomicResonanceModel` is the shared Morse + sigmoid + Gaussian-in-r form
eMoScat uses for every diatomic in its resonance model (N2, NO, F2 -- see
`qscat.model.library`): a Morse neutral curve `v0(R)`, a sigmoid interaction
strength `lambda(R)`, and a Gaussian-in-r electron-molecule interaction
`V_int(r, R) = -lambda(R) exp(-alpha_c r^2)`. Ported verbatim (formula only,
not shared code) from `projects/n2_resonance/potential.py`'s
`v0`/`lam`/`v_int` and `projects/n2_2d_cross_section/hamiltonian2d.py`'s
`potential_2d`/`interaction_diag`/`build_h2d`, parameterized by dataclass
fields instead of a `config.json`/`PARAMS` dict lookup. Consistency with
those (still-live, not-yet-rewired) N2 project modules is NOT guaranteed by
construction -- it is guaranteed BY TEST: `libs/qscat/tests/test_model.py`
cross-checks this module against them elementwise to 1e-14.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd

__all__ = ["ResonanceModel", "DiatomicResonanceModel"]


@runtime_checkable
class ResonanceModel(Protocol):
    """Structural type: anything with this shape is a model `qscat.core` can
    drive, regardless of functional form.

    `mu`/`ell`/`charge` are per-molecule constants (nuclear reduced mass,
    fixed partial wave, Coulomb charge of the residual channel -- 0 for a
    neutral target, -1 for a cation like H2+). `v0`/`v_int`/`surface` are the
    potential-energy surface, evaluated pointwise (scalars or broadcastable
    arrays; `r`/`R` may be COMPLEX on an exterior-complex-scaling tail, and an
    implementation must not coerce them to a real dtype). `hamiltonian`/`interaction_diag`
    assemble those surfaces onto a `qscat.dvr.TensorGrid`: the sparse `H_2D`
    and the flattened `diag(V_int)` a time-dependent free-reference subtracts
    off `H_2D`.

    `mu`/`ell` are declared as READ-ONLY properties (not plain mutable
    attributes) so that a frozen dataclass like `DiatomicResonanceModel` --
    whose fields are read-only by construction -- satisfies this protocol
    structurally under mypy. A plain `mu: float` annotation here would
    require a SETTABLE attribute (mypy's default for Protocol attributes),
    which a frozen dataclass field can never be; `qscat.core.driven`'s
    `model: ResonanceModel` parameter is the first real static consumer of
    this protocol (earlier code only used `isinstance` at runtime, which
    does not distinguish read-only from settable), so this mismatch was
    latent until Task 4.
    """

    @property
    def mu(self) -> float: ...
    @property
    def ell(self) -> int: ...
    @property
    def charge(self) -> int:
        """The Coulomb charge z for the channel functions; 0 for neutral
        targets, -1 for a singly-charged cation like H2+."""
        ...

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Neutral molecule's own potential (Hartree); a CHANNEL potential --
        survives as r -> infinity."""
        ...

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """`V_int(r, R) = -lambda(R) exp(-alpha_c r^2)` -- the perturbation
        ALONE, excluding `v0(R)` and the centrifugal term."""
        ...

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """The full surface `v0(R) + ell(ell+1)/(2 r^2) + v_int(r, R)` that
        goes into `hamiltonian`."""
        ...

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """`H_2D` on `tgrid` (axis 0 electronic `r`, axis 1 nuclear `R`),
        sparse CSR, complex symmetric (never Hermitian) whenever an axis has
        an ECS tail."""
        ...

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        """`V_int` evaluated on `tgrid`, flattened (C order) -- what a
        time-dependent free-reference subtracts off `hamiltonian(tgrid)`."""
        ...


@dataclass(frozen=True)
class DiatomicResonanceModel:
    """The shared Morse + sigmoid + Gaussian-in-r resonance-model form,
    implementing `ResonanceModel`. Carries no grid state -- purely the
    potential-surface parameters; `qscat.model.library` holds the
    per-molecule instances.
    """

    mu: float
    ell: int
    D0: float
    alpha0: float
    R0: float
    lambda_inf: float
    lambda_1: float
    R_lambda: float
    lambda_c: float
    R_c: float
    alpha_c: float
    charge: int = 0

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Neutral Morse potential (Hartree). Minimum -D0 at R0.

        `R` is converted to `complex128`, never a real dtype: it may be
        complex (ECS-rotated tail points), and coercing to real would
        silently discard Im(R) and corrupt the analytic continuation the
        exterior-complex-scaling method relies on.
        """
        Rc = np.asarray(R, dtype=np.complex128)
        a = self.alpha0
        out = self.D0 * (np.exp(-2 * a * (Rc - self.R0)) - 2 * np.exp(-a * (Rc - self.R0)))
        return np.asarray(out, dtype=np.complex128)

    def lam(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Sigmoid interaction strength; `lam(R_c) == lambda_c`."""
        Rc = np.asarray(R, dtype=np.complex128)
        l1, Rl = self.lambda_1, self.R_lambda
        lam0 = (self.lambda_c - self.lambda_inf) * (1 + np.exp(l1 * (self.R_c - Rl)))
        out = self.lambda_inf + lam0 / (1 + np.exp(l1 * (Rc - Rl)))
        return np.asarray(out, dtype=np.complex128)

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Electron-molecule interaction potential (Hartree).

        `r` is converted to `complex128`, never a real dtype: it may be
        complex (ECS-rotated tail points), and coercing to real would
        silently discard Im(r) and corrupt the analytic continuation the
        exterior-complex-scaling method relies on.
        """
        rr = np.asarray(r, dtype=np.complex128)
        out = -self.lam(R) * np.exp(-self.alpha_c * rr**2)
        return np.asarray(out, dtype=np.complex128)

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """The full surface `v0(R) + ell(ell+1)/(2 r^2) + v_int(r, R)`.

        Must not coerce to a real dtype: `r`/`R` are complex on the ECS
        tails.
        """
        rr = np.asarray(r, dtype=np.complex128)
        out = self.v0(R) + self.ell * (self.ell + 1) / (2.0 * rr**2) + self.v_int(rr, R)
        return np.asarray(out, dtype=np.complex128)

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """`H_2D` on `tgrid` (axis 0 = electronic r, axis 1 = nuclear R)."""
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        """`V_int` evaluated on the tensor grid, flattened (C order)."""
        return potential_nd(tgrid, self.v_int)
