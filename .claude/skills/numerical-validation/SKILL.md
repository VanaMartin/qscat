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

**REQUIRED BACKGROUND:** `superpowers:test-driven-development` — write these
checks as tests before/alongside the implementation, not after the fact.

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
    np.testing.assert_allclose(E0, 0.5, rtol=1e-8)
```

### (b) Convergence studies

Refine the discretization (grid spacing, number of basis functions, DVR
points) across a few resolutions and assert the error decreases monotonically
toward the known convergence rate (e.g. spectral methods should show
super-algebraic decay; finite differences should show the expected power
law). A single high-resolution pass proves nothing about correctness of the
discretization — you need at least 3 resolutions to see the trend.

```python
errors = [abs(compute(n) - exact) for n in (8, 16, 32, 64)]
assert all(e2 < e1 for e1, e2 in zip(errors, errors[1:])), errors
```

### (c) Conservation checks

For time evolution, assert the invariants the physics guarantees: norm
preservation (`||psi(t)|| == ||psi(0)||`) and unitarity of the propagator
(`U @ U.conj().T ≈ I`). A propagator that drifts in norm over long
integration is a correctness bug even if short-time results look fine.

```python
np.testing.assert_allclose(np.linalg.norm(psi_t), 1.0, atol=1e-10)
```

### (d) Differential testing

Compare against an independent implementation:
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
- Skipping conservation checks for time-evolution code because "the energy
  looked right" — energy accuracy and norm conservation are independent
  failure modes.
- Writing the validation after the implementation instead of alongside it —
  follow `superpowers:test-driven-development`.
