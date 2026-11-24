# Real nuclear box with electronic ECS: a projected toy Hamiltonian

This project-level construction describes an explicitly imposed real nuclear
Dirichlet box, not an outgoing nuclear resonance or physical molecular scattering
certificate. Electronic exterior complex scaling (ECS) is allowed. No extension
of a potential beyond its real nuclear input domain is implied.

In atomic units let the real orthonormal DVR coefficient columns of `B` solve
`(T_R + diag(V0)) B = B diag(eps)`. Project only the nuclear coordinate:

```
H = T_radial ⊗ I + I ⊗ diag(eps)
    + blockdiag_i [B.T @ diag(U(r_i, R)) @ B].
```

Flatten electronic index first, nuclear state second (C order). `T_radial` may
include the electronic centrifugal diagonal, which must be added exactly once.
The nuclear eigenvalues already include `V0`; do not add it again. Since `B`
contains normalized DVR **coefficients**, no further quadrature weights enter
the projection. `B.T`, not a conjugated interaction, preserves complex symmetry
when the electronic contour and `U` are complex. The basis itself is real.

`projects.potential_factory.nuclear_projection.projected_hamiltonian` accepts
the radial matrix, neutral eigenvalues, basis columns, and sampled interaction
matrix `(n_radial, n_nuclear_DVR)`; it returns a complex CSR matrix. It is an
assembler, not an eigensolver, pole selector or convergence classifier.

Analytic separable controls require pairwise sums of radial and nuclear levels.
A separate dense full-DVR matrix projected with `I ⊗ B` checks ordering, weights,
complex transpose and couplings. A complete orthogonal nuclear basis must be
isospectral to that full finite-DVR matrix. These tests establish assembly, not
convergence of a truncated basis or an ECS contour.

Compare eigenvalue spectra as unordered multisets with a one-to-one assignment,
not by lexicographically sorting complex numbers. Equal-real-part pairs can swap
order after last-bit perturbations across LAPACK implementations. Keep the same
error tolerances and check multiplicity; independent nearest-neighbor matches
alone can incorrectly reuse one reference eigenvalue twice.

For any later pole study, expose **channel** momenta relative to `eps_n`, not
`sqrt(2 E)` relative to a fixed-nuclei electronic zero. Open outgoing and closed
decaying branches need explicit continuous sheet choices. A small eigen-residual
or nearest-shift match alone cannot identify a resonance. Check independent
identity/localization witnesses and electronic angle, order, tail and nuclear
basis convergence. A BO curve from another representation can supply a seed,
but is not a comparator measuring nonadiabatic error of this interaction.

No nuclear ECS-angle test applies to an imposed real closed box. Do not bypass
the two-coordinate outgoing checks of `qscat.core.exact_resonance_states` by
passing an identical nuclear partner. This toy does not promote a new `qscat`
API; optimization requires profiling and promotion requires independent review.

## An oracle-identified electronic reference is a single-target problem

A product reference may need one electronic vector at an independently known
outgoing pole, rather than an unlabelled set of nearby continuum eigenpairs.
Requesting one shift-invert pair is appropriate only with explicit identification:
agreement with that pole in position and full width, a small eigen-residual,
nonvanishing c-self-norm and localization. It does not replace the multi-pair
search used to distinguish candidates in the coupled problem.

Before replacing a stopped multi-pair reference prescription, test the targeted
pair against a dense all-eigenpair solve on the stopped grid: require a unique
nearest dense pole, energy agreement and bilinear vector overlap. The dense
comparison shares discretization; an independent radial outgoing oracle supplies
the separate boundary/energy check. Keep the original failure and unchanged
physical accuracy gates. Convergence of one unexported partial ARPACK pair alone
cannot justify this replacement.
