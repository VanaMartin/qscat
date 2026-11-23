# A p-wave radial boundary oracle for synthetic Gaussian wells

This toy comparator integrates a scalar radial ODE, independently of FEM-DVR
assembly and ECS eigensolvers. It diagnoses boundary errors; it is not a molecular
potential fitter. Atomic units and electron mass one are used throughout.

For `ell=1`, `V=-lambda exp(-alpha r²)+1/r²`, write the regular radial solution
as `u=r² f`. Its equation is

    f'' + 4 f'/r + 2(E + lambda exp(-alpha r²)) f = 0.

A regular origin start uses `f=1+c2 r²+c4 r⁴`, where
`c2=-(E+lambda)/5` and `c4=((E+lambda)²/5+lambda alpha)/14`.
The implementation starts at `1e-4 bohr`; the omitted terms and integration
tolerance must be negligible in the regime being checked. `regular_p_wave`
returns `(f, f')` at a declared real matching radius, not a normalized probability.
For a zero well, an independent analytic control is `f=3 j1(kr)/(kr)`.

With `k=sqrt(2E)` on the branch with positive real part, unnormalized free radial
waves are `h±(r)=exp(±ikr)(1±i/(kr))`. Their derivatives are taken analytically,
and their Wronskian `h+ h-' - h+' h- = -2ik` is nonzero. Normalization factors
cancel in boundary matching. The outgoing wave grows on the real axis for a
resonance; ordinary real-axis square-integrability is not imposed.

`p_wave_boundary` returns the outgoing value/derivative, or the finite-endpoint
solution `b(r)=h+(r)h-(Z)-h-(r)h+(Z)` whose value at the complex Dirichlet endpoint
`Z` vanishes. Matching the regular solution to `b` with a normalized Wronskian
allows a root finder to compare the infinite outgoing condition with a specific
finite ECS endpoint without another ECS grid calculation. The finite-endpoint
combination need not be purely outgoing.

`gaussian_p_wave_mismatch` neglects the Gaussian beyond its matching radius but
retains the centrifugal term exactly. Its residual is dimensionless and scaled
by both solution-vector norms; it is not a bound on energy error. Tail omission,
origin truncation, root conditioning and integration convergence must be assessed
separately. These routines support neither general partial waves nor Coulomb or
arbitrary molecular tails. Finite-endpoint agreement tests shared Hamiltonian and
boundary conventions, not the accuracy of either as a physical CO model.
