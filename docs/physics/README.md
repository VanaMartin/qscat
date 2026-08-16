# Theory notes

One note per method: the derivation, the equations, the unit conventions
(atomic units throughout), the validation evidence, and the literature it
comes from. These are the working notes behind the implementation, so they
record limitations and negative results as well as what works.

## Discretisation

- {doc}`femdvr-ecs` — the FEM-DVR grid with an exterior-complex-scaled tail,
  and the four analytic benchmarks that pin it down.
- {doc}`nd-tensor-hamiltonian` — the N-dimensional sparse tensor Hamiltonian.
- {doc}`discretisation-tuning` — deriving a grid from the potential instead of
  hand-picking element lengths, and where the 1-D probes are not sufficient.
- {doc}`mumps-sparse-backend` — the complex-symmetric MUMPS backend.
- {doc}`ti-energy-sweep-reuse` — reusing the symbolic factorization across an
  energy sweep.
- {doc}`shift-invert-eigensolver` — the eigenpairs nearest a complex shift, for
  resonances that sit in the interior of the spectrum. Validated in 1-D only.

## The scattering engine

- {doc}`qscat-core-scattering` — the model-independent engine and the
  model/engine split it enforces.
- {doc}`n2-resonance` — locating the resonance pole.
- {doc}`n2-cross-section` — the one-dimensional time-independent route.
- {doc}`n2-2d-cross-section` — the exact two-dimensional driven solve, gated
  against independent published data.

## Time-dependent routes

- {doc}`n2-td-cross-section` — wavepacket propagation in one dimension.
- {doc}`n2-2d-td-cross-section` — the exact two-dimensional time-dependent
  route, including why order-1 Crank–Nicolson was not enough.
- {doc}`td-extractors` — three energy extractors sharing one propagation.
- {doc}`td-da` — the dissociative-attachment generalization.

## Molecules and approximations

- {doc}`diatomic-ve-cross-sections` — NO and F₂, and the local-complex-potential
  approximation measured against the exact oracle.
- {doc}`h2plus-dr` — dissociative recombination for an ionic target.
- {doc}`lcp-resonance-levels` — Born–Oppenheimer quasi-bound levels in the
  complex curve.
- {doc}`exact-2d-resonances` — the same levels without the approximation: poles
  of the full 2-D S-matrix, and what the Born–Oppenheimer error actually
  measures on N₂.

## Open directions

- {doc}`angular-coupled-channels` — the parked angular extension.
- {doc}`optimization-targets` — where the remaining hot paths are.
