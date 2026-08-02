# MKL PARDISO backend for `qscat.linalg.SparseLU` — Design Spec

**Date:** 2026-08-02
**Status:** Approved design (derived from the hardening roadmap, Part 5 #1)
**Lifecycle:** optimize stage — a new solver backend behind an existing dispatch;
the SciPy SuperLU path stays as the differential oracle.

## Goal

Add `backend="pardiso"` to `qscat.linalg.SparseLU`: Intel MKL PARDISO as a third
direct sparse backend alongside SuperLU (`scipy`) and MUMPS (`mumps`), for the
complex-symmetric ECS matrices. This is the eMoScat reference solver (it did all
of N2/NO/F2 in <1 hr) and the measured #1 hot path — profiling shows **~98% of
TI time is the numeric factorization** (`docs/physics/optimization-targets.md`).

## Why

`SparseLU` already abstracts the backend; MUMPS beats SuperLU 72×/9× on the
production decks. PARDISO is the other top-tier complex-symmetric direct solver
and the historical reference; benchmarking it closes the loop on "are we at the
eMoScat bar?" and gives x86 users (the `cpu-mkl` Docker base) a fast path without
system MUMPS.

## Non-goals

- No API change beyond adding the backend name. `auto` resolution order and
  `set_default_backend`/`default_backend` semantics are unchanged.
- Not the GPU backend (separate spec).
- No change to `refactor` / symbolic-reuse semantics (PARDISO must honor them).

## Design

- New `libs/qscat/qscat/linalg/_pardiso_backend.py`, mirroring `_mumps_backend.py`:
  a factor object exposing `solve(b)` and honoring the complex-symmetric
  structure (PARDISO `mtype=6`, complex symmetric). Access via `pypardiso` if it
  supports complex, else a thin `ctypes`/MKL `pardiso` binding (decision point:
  evaluate `pypardiso` complex support first; the MUMPS backend used
  `python-mumps`, so prefer an existing wheel).
- Wire into `SparseLU.__init__` / `refactor` dispatch and `pardiso_available()`
  (mirror `mumps_available()`). `auto` stays `mumps → scipy`; PARDISO is opt-in
  by name (revisit `auto` order after benchmarking).
- Raise `qscat.exceptions.BackendError` when requested-but-unavailable (mirror
  the MUMPS path).
- `refactor(reuse_analysis=True)` must reuse PARDISO's symbolic phase (phase 11
  once, phase 22 per matrix) — the energy-sweep reuse `docs/physics/
  ti-energy-sweep-reuse.md` relies on.

## Validation

- **Differential oracle:** on the same matrices as `test_sparse_lu.py`, PARDISO
  `solve` matches SuperLU to `rtol=1e-10` (complex-symmetric RHS + the N2 2-D
  matrix). `@skipif` when MKL/PARDISO absent (mirrors the MUMPS Mac skips).
- **refactor equivalence:** a diagonal-shift sweep (`E·I − H`) via
  `refactor(reuse_analysis=True)` matches fresh-per-energy to round-off.
- **Benchmark:** extend `benchmarks/mumps_vs_superlu.py` to a three-way
  factor-time / solve-time / peak-RSS table on the real N2 decks (Docker
  `cpu-mkl` image). Report vs the eMoScat <1 hr bar.

## Deliverables

- `_pardiso_backend.py` + dispatch wiring + `pardiso_available()`.
- `pardiso` extra in `libs/qscat/pyproject.toml` (+ the `cpu-mkl` Docker base
  provisions MKL); tests `@skipif`-absent on a laptop, run in Docker.
- Three-way benchmark + a short `docs/physics/pardiso-sparse-backend.md`.

## Verification

`uv run pytest -k pardiso` green in the `cpu-mkl` Docker image; differential
+ refactor tests pass; benchmark table committed; `auto` resolution documented.
