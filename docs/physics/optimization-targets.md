# Optimization targets (measured)

Roadmap Part 5 says "profile before cutting." `benchmarks/profile_hotpaths.py`
does that. This note records what it measured and the resulting plan.

## Measured hot path (TI)

`uv run python -m benchmarks.profile_hotpaths` on a representative N2 2-D problem
(electronic r_max=16 × nuclear R=22, order 6, 9-energy VE sweep):

```
21368 function calls in 14.313 seconds (cumulative)
  14.312  core/driven.py ve_cross_section
  14.083  scipy .../linsolve.py splu
  14.081  {scipy.sparse.linalg._dsolve._superlu.gstrf}   <-- 98.4% of total
  12.507  linalg/sparse_lu.py refactor  (per-energy, symbolic reuse)
   1.587  linalg/sparse_lu.py __init__  (first factorization)
```

**Finding: ~98% of the TI cost is the SuperLU numeric factorization** (`gstrf`),
reached through `qscat.linalg.SparseLU`. Everything else — grid assembly, the
per-energy solve, the cross-section projection — is noise by comparison. This
confirms, with numbers, that the sparse factorization is the single dominant
cost, exactly as the roadmap assumed.

## Ranked plan (unchanged by the measurement, now evidence-backed)

1. **Sparse factorization (`SparseLU`) — #1, by 60×.** MUMPS already beats
   SuperLU 72× in factor time / 9× in peak RSS on the production N2 matrices
   (`docs/physics/mumps-sparse-backend.md`), and `refactor` reuses the symbolic
   analysis across an energy sweep. Next: an **MKL PARDISO** backend (the eMoScat
   reference solver; the `cpu-mkl` Docker base is scaffolded for it) as a third
   dispatch option, and a **GPU sparse** backend (cuDSS/cuSOLVER) as a fourth —
   both slot into the existing `SparseLU` backend dispatch without touching call
   sites.
2. **The propagation inner loop** — profile with `--td`: thousands of
   cached-factor triangular solves + the extractor `record` projections. The
   ideal shape for the **first real Rust kernel** (high call count, small hot
   code, clean differential oracle).
3. **`mpmath` Coulomb functions** (`qscat.special.coulomb`) — the H2P DR Rydberg
   loop; vectorize/cache or move to a fast implementation.
4. **Hamiltonian assembly** — one-time; only worth it for the very largest decks.

The Rust kernel crate (`native/qscat-kernels`) is currently a stub; per the
lifecycle, the Python path stays as the differential oracle for every kernel.
