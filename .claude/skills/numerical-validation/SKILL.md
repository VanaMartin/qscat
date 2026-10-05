---
name: numerical-validation
description: Use when validating quantum/numerical code where exact equality does not apply — analytic benchmarks, convergence studies, conservation checks, and differential testing against references.
---

# numerical-validation

## Overview

Floating-point quantum-mechanics code is almost never checked with bare
`==`. Correctness is established by comparing against known analytic
results, watching error shrink at the expected rate as resolution improves,
checking physical invariants the math guarantees, and/or matching an
independent reference implementation. Use one or more of the four techniques
below depending on what the method is; most methods need at least two.

Write validation checks before or alongside the implementation. Use
`superpowers:test-driven-development` when available; its absence does not remove
the requirement to establish the expected result independently of the code.

## When to Use

- Writing or reviewing tests for any physics/numerical routine (special
  functions, DVR grids, ECS contours, time evolution, linear algebra).
- Step 3 ("Validate") of `qm-method-lifecycle`.
- Deciding whether a Rust kernel (`python-to-rust-kernel`) matches its Python
  oracle.
- A PR touches numerical code and it's unclear what "correct" even means for it.

## The Four Techniques

### (a) Analytic benchmarks

Compare against a closed-form result: harmonic oscillator eigenvalues
(`E_n = n + 1/2` in atomic units), hydrogen energy levels
(`E_n = -1/(2n^2)` Hartree), Coulomb phase shifts, etc. Use
`pytest.approx(expected, rel=..., abs=...)` or
`np.testing.assert_allclose(actual, expected, rtol=..., atol=...)` — never a
bare `==` on floats.

```python
import numpy as np

def test_harmonic_oscillator_ground_state():
    E0 = solve_ho_eigenvalue(n=0)
    np.testing.assert_allclose(E0, 0.5, rtol=1e-8, atol=0.0)
```

### (b) Convergence studies

Refine the discretization (grid spacing, basis size, DVR order, time step, or
domain extent) across at least three resolutions. Track the actual observable
and, where available, its error against an independent result. State the expected
trend or convergence rate and the regime in which it applies.

Strict monotonic decrease is appropriate only when the method and refinement
family justify it. Resonance tracking, competing error sources, and roundoff can
produce non-monotone sequences or a plateau. Explain these with resolved states,
separate refinement knobs, and an error floor; do not turn an unexplained plateau
into a pass. Report the sequence and confirm the final observable meets the stated
error budget. A single high-resolution pass or one small successive difference
does not establish convergence.

### (c) Conservation checks

Choose invariants from the Hamiltonian, boundary conditions, and inner product:

- **Closed Hermitian systems:** test physical norm preservation and propagator
  unitarity in the appropriate metric. For a time-independent Hamiltonian, also
  check energy conservation when applicable. Quantify integration error over the
  relevant propagation interval; short-time agreement can hide accumulated drift.
- **ECS, absorbing boundaries, or other non-Hermitian systems:** ordinary norm
  preservation and unitarity are not general requirements. Check the applicable
  continuity/flux balance, real-region observables, stability, and convergence
  against an independent propagation or scattering oracle. Loss can be physical;
  apparent growth needs investigation but is not classified solely by its sign.
  Monotone decay requires additional dissipativity assumptions and is not a
  universal ECS invariant.

Preserve the bilinear, non-conjugated c-product for ECS complex-symmetric algebra;
it is not an ordinary probability norm. State explicitly which product and region
each measurement uses. Apply `assert_allclose` only to an invariant justified by
these assumptions, with both `rtol` and `atol` supplied.

### (d) Differential testing

State what is independent about the oracle and what the comparison establishes.
Two routes sharing a grid or ingredient can validate their agreement while sharing
the same discretization or normalization error. Separate model approximation,
discretization, and reference uncertainty in the error budget.

Possible comparisons:
- vs. `reference/` (`reference/eMoScat`, `reference/libXcuda`) — read-only
  oracles; never edit them, only read outputs/algorithms to compare against.
- vs. `mpmath` high-precision arithmetic for special functions or anywhere
  double precision itself is in question.
- vs. the retained Python implementation, when checking a Rust kernel (see
  `native/qscat-kernels/tests/test_l2_norm.py` for the exact pattern: same
  inputs into both implementations, assert the results agree within a stated
  tolerance).

## Tolerance Conventions

- Always state `rtol`/`atol` explicitly — don't rely on a tool's default.
- Never compare floats with bare `==`, except for kernel outputs that are
  provably bit-exact for a specific trivial case (e.g. an exact integer-valued
  input), and even then prefer `assert_allclose` for consistency.
- Looser tolerances (e.g. `rtol=1e-6`) are fine for iterative/convergent
  methods; tight tolerances (`rtol=1e-12` or better) are expected for
  differential tests between two implementations of the same deterministic
  arithmetic (e.g. Python vs. Rust on the same inputs).
- See `qscat-conventions` for the repo's default tolerance values.

## Common Mistakes

- Asserting only a single resolution converges "close enough" instead of
  running a convergence study.
- Comparing against `reference/` output but not pinning down what tolerance
  is acceptable given the reference's own precision.
- Applying closed-system norm or unitarity checks to ECS/open systems, or skipping
  the balance and stability checks justified by the actual problem.
- Writing validation after the implementation instead of alongside it, or deriving
  both the expected value and implementation from the same unverified assumption.
