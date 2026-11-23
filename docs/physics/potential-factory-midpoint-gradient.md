# Differentiate the reported two-angle midpoint

For two fixed FEM-DVR-ECS grids, let the selected simple eigenvalues be
`Ea(lambda, alpha)` and `Eb(lambda, alpha)`. The factory reports
`z = (Ea + Eb)/2`; its derivative is therefore

    dz/dp = (dEa/dp + dEb/dp)/2.

For a complex-symmetric Hamiltonian, normalize each right eigenvector separately
by its bilinear c-product. Then `dE/dV_i = psi_i**2`, without conjugation. For
`V(r) = -lambda exp(-alpha r**2) + ell(ell+1)/(2r**2) + extra(r)`, with fixed grids
and parameter-independent `extra`, each angle contributes

    dE/dlambda = sum_i psi_i**2 (-exp(-alpha r_i**2)),
    dE/dalpha  = sum_i psi_i**2 (lambda r_i**2 exp(-alpha r_i**2)).

Use that angle's own complex coordinates before averaging the two scalar
derivatives. Averaging eigenvectors or using angle-A coordinates for angle B is
not equivalent. Coordinates are bohr; lambda is Hartree and alpha is bohr^-2.

`Pole.angle_energies` retains the literal selected pair. Legacy energy/residual
construction remains valid as a value holder, but `well_pole_gradient` rejects
missing or inconsistent pair provenance. It follows each matrix's own selected
eigenvalue, not whichever eigenvalue is nearest the midpoint. The existing
single-matrix `pole_sensitivity` API retains its original meaning.

The shared reduced gradient supplies both Newton and joint-polish Jacobians.
This establishes estimator consistency, not a physical error bound. Degenerate,
self-orthogonal or changing-identity branches do not gain a derivative certificate
from this formula. Candidate admission and continuity controls remain unchanged;
forward pole convergence and parameter identifiability remain separate gates
before a molecular fit. No global potential or QSCAT promotion is implied.
