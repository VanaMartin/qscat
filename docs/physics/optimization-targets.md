# Optimization targets (measured)

Roadmap Part 5 says "profile before cutting." `benchmarks/profile_hotpaths.py`
does that. This note records what it measured and the resulting plan.

## Key result

Three measured findings set the agenda. **(1) The sparse LU is the whole
story:** ~98% of the TI cost is the SuperLU numeric factorization, and the
TD propagation is ~82% per-step triangular solves. **(2) MUMPS fixes
both:** 8.2× on the factor and 4.6× on the per-solve at N=20328, so MUMPS
is the correct default backend for TI and TD alike — and
`SparseLU(backend="auto")` already selects it where provisioned. On the
current direct-solver architecture that is near the practical optimum;
PARDISO and GPU cuDSS are incremental beyond it, not step changes.
**(3) There is NO pure-Python hot loop worth a first Rust kernel:** the
`c_product`/extractor loops the first-kernel spec targeted are ~0.1% of
runtime, so that spec's premise is invalidated. Separately: thread
oversubscription cost ~300× on a concurrent sweep — pin BLAS threads per
worker whenever processes are multiplied.

## Measured hot path (TI)

`uv run python -m benchmarks.profile_hotpaths` on a representative N2 2-D problem
(electronic r_max=16 × nuclear R=22, order 6, 9-energy VE sweep):

```text
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

## Measured hot path (TD) — and a correction

`profile_hotpaths.py --td` on the same problem (an 800-step order-3 Padé
propagation) tells the same story, with a twist:

```text
51.6 s total
  42.2 s (82%)  {method 'solve' of 'SuperLU' objects}   <-- per-step solve
   8.6 s (17%)  SparseLU factorization (Padé poles)
   ~0.05 s      c_product / extractor record (pure Python)
```

Two consequences:

1. **There is NO pure-Python or Rust-portable hot loop.** The `c_product` /
   extractor `record` loops the "first Rust kernel" spec targeted are ~0.1% of
   runtime; porting them to Rust would save nothing. That spec's premise is
   invalidated by this profile — the first-kernel opportunity does not exist for
   the current direct-solver architecture (it would only appear for a future
   iterative-solver or very-large-scale-assembly path).
2. **For TD, the per-step *solve* dominates the *factor* by ~5×.** This matters:
   MUMPS's headline win is *factorization* (72× — `docs/physics/
   mumps-sparse-backend.md`); whether its *solve/back-substitution* beats
   SuperLU's is a separate, unmeasured question that decides the optimal TD
   backend. `benchmarks/solve_throughput.py` measures exactly this.

## Measured: does MUMPS fix the solve-bound TD path? Yes.

`benchmarks/solve_throughput.py` on a CN/Padé shift of the N2 2-D Hamiltonian
(N=20328), factor once + 100 repeated solves, in the Docker `test` image:

| backend | factor | per-solve | 100 solves |
|---|---|---|---|
| SuperLU | 2.54 s | 15.8 ms | 1.58 s |
| **MUMPS** | **0.31 s (8.2×)** | **3.4 ms (4.6×)** | **0.34 s** |

MUMPS wins **both** — crucially the *solve* by ~4.6×, not just the factor. So the
solve-bound TD propagation is ~4.6× faster on MUMPS per step, and the whole
propagation (solve + factor) is ~5× faster. **MUMPS is the correct default for
both TI and TD**, and `SparseLU(backend="auto")` already selects it when present
(it is provisioned in the Docker `test`/`cpu` images; absent on a bare Mac, hence
these numbers are Docker-measured). Conclusion: on the current architecture we
are near the practical direct-solver optimum once MUMPS is used; PARDISO and GPU
cuDSS are *incremental* beyond this (parity/large-scale), not step changes — the
big win over the SuperLU fallback is already in hand.

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

## Thread oversubscription costs ~300× on concurrent sweeps (2026-08-17)

A dense figure-regeneration sweep that should take seconds took hours. The cause
was not the solver, not the backend, and not the VE path: it was running several
sweeps **concurrently** with unpinned multithreaded BLAS/MUMPS. Measured on the
same machine (32 cores), same image, same `backend: auto` (MUMPS), same configs:

| sweep | grid unknowns | energies | 3 concurrent, unpinned | 1 at a time, 8 threads | speedup |
|---|---|---|---|---|---|
| N₂ VE | 26 857 | 196 | 2899.8 s | **10.4 s** | **279×** |
| F₂ VE | 128 568 | 97 | >6840 s (killed) | **26.4 s** | **>259×** |
| NO VE | 78 804 | 117 | >6840 s (killed) | **20.8 s** | **>330×** |

Three processes on 32 cores cannot cost 300× by fair sharing. Each container's
BLAS/MUMPS sized its thread pool to the whole machine, so roughly 3 × 32 threads
contended for 32 cores; spin-waiting and memory-bandwidth saturation degrade a
sparse factorization catastrophically rather than proportionally. Load average
sat at 50–80 throughout.

**Operational rule: run sweeps one at a time, and pin the thread count.**

```bash
docker run --rm -e OMP_NUM_THREADS=8 -e OPENBLAS_NUM_THREADS=8 \
  -e MKL_NUM_THREADS=8 -e NUMEXPR_NUM_THREADS=8 ...
```

Correctness was unaffected: the N₂ sweep reproduces Houfek identically either way
(ratios 0.9996 / 1.0003 / 0.9994 at E = 0.10 from both the 2899 s and the 10.4 s
run). The contention cost wall-clock, not accuracy.

### Two hypotheses this disproved

Both looked plausible from the slow numbers, and both were wrong:

1. **"The VE path doesn't reuse the symbolic factorization."** It does.
   Instrumenting `SparseLU` shows VE and DA are identical — 4 energies gives
   `__init__` (full analysis) ×1 and `refactor` ×3 on *both* paths.
2. **"VE is intrinsically ~130× slower than DA on the same grid."** An artefact
   of the contention. Run properly, the per-energy costs are the same order:
   F₂ VE 0.27 s vs F₂ DA 0.53 s (identical 128 568 grid); NO VE 0.18 s vs NO DA
   0.37 s. VE is, if anything, cheaper per energy.

The lesson worth keeping: a 300× wall-clock anomaly invited structural
explanations, and two were constructed before the execution environment was
ruled out. Vary the run conditions before concluding anything about the code.

### Still open (small, and now much cheaper to test)

- Per-energy cost is **flat** on SuperLU — measured 3.16 s/energy at n = 1, 2, 4
  on the N₂ deck — because `scipy` re-runs `splu` with no reuse. That is
  documented behaviour, but it means the default laptop backend gets nothing
  from `SparseLU.refactor`, whose measured ~5× win is MUMPS-only.
- The thread count was not tuned: 8 was chosen to leave headroom on a shared
  machine, not measured as optimal.
