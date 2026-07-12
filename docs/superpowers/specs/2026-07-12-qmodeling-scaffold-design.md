# qModeling Scaffold — Design Spec

**Date:** 2026-07-12
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending

## Context

`qModeling` is a new, ground-up quantum-mechanics research monorepo. It exists to
**resuscitate, clean up, and extend** prior quantum-scattering work (the C++/CUDA
`eMoScat` / QSCAT package) under a modern, Python-first, HPC-capable methodology.

The prior codebase (`~/src/eMoScat`) implements electron–molecule scattering with
FEM-DVR + Exterior Complex Scaling (ECS), Coulomb/Bessel special functions
(`coulcc.f`), Chebyshev methods, Crank–Nicolson time evolution, conjugate gradients,
Møller operators, and LCP/NRM modules. It is C++ (61 `.cpp`, 64 `.h`, 26 `.hpp`),
depends on a CUDA layer (`libXcuda`), OpenBLAS, and the Intel compiler, and carries a
Boost-based Python wrapper. Its README records known pain: a wrapper redesign, an
Input-library redesign, "Segfault on inconsistent LCP load", and "Evolution FAIL!".
C++ is judged increasingly obsolete for this work; the numerics are worth preserving,
the architecture is not.

This first phase does **not** port any physics. It delivers the **scaffold** — repo
structure, Claude Code skills, agents, and the operating manual (`CLAUDE.md`) — that
makes the intended development lifecycle the path of least resistance. Porting happens
in later, per-method phases.

## Development lifecycle (the organizing principle)

Every new capability follows five stages:

1. **Design** — imagine and specify the problem (physics + numerics + interfaces).
2. **Toy model in Python** — pure Python scaffolding and simple libs, correctness first.
3. **Validate** — prove it right (analytic benchmarks, convergence, conservation,
   differential tests vs. references).
4. **Optimize** — move computationally heavy paths into an optimal compiled language
   (default: Rust via PyO3/maturin), keeping the Python version as the oracle.
5. **Promote** — graduate validated, reusable code into the standard library (`qscat`).

Two invariants hold at every stage: **it runs locally on CPU**, and **it can be
containerized** for reproducible runtime.

## Tech stack (decisions)

- **Research/user language: Python.** Environment, workspace, and builds via **`uv`**
  (fast, fully reproducible lockfiles for HPC/Docker parity).
- **Compiled kernel language: Rust** (default), via **PyO3 + maturin**. Chosen per
  package; Rust is the default for its memory safety, Python interop, and clean static
  builds. The door stays open to alternatives when a package warrants it.
- **Intermediate optimization rung (optional): Numba/Cython** — validate a speedup
  inside Python before committing to a Rust kernel.
- **Numerics:** numpy / scipy; **mpmath** for high-precision reference oracles
  (Coulomb/Bessel); matplotlib for diagnostics.
- **Testing:** pytest + **hypothesis** (property-based) + numpy tolerance helpers +
  **pytest-benchmark**. Rust: `cargo test` + **criterion**.
- **Quality:** ruff + mypy (Python); clippy + rustfmt (Rust); pre-commit hooks.
- **Runtime:** multi-stage **Docker** (build Rust wheel → slim CPU runtime), AWS-ready
  from day one. AWS deployment is deferred to its own later phase.
- **GPU/CUDA:** deferred. `libXcuda` (recovered at
  `ssh://git-codecommit.eu-central-1.amazonaws.com/v1/repos/libXcuda`) is a read-only
  reference for eventual GPU kernels, not a current dependency.

## Repository structure

```
qModeling/
├── CLAUDE.md                 # operating manual: the lifecycle, layout, conventions
├── README.md
├── pyproject.toml            # uv workspace root
├── .claude/
│   ├── skills/               # custom domain skills (see below)
│   ├── agents/               # custom subagents (see below)
│   └── settings.json         # permissions + hooks
├── libs/
│   └── qscat/                # THE standard library (validated, reusable)
│       └── qscat/{special,dvr,ecs,evolution,linalg,units}/
├── native/
│   └── qscat-kernels/        # Rust crate(s) → Python wheels (PyO3/maturin)
├── projects/                 # per-problem research + toy models (lifecycle stages 1–2)
├── validation/               # analytic benchmarks, golden data, convergence studies
├── docs/
│   ├── superpowers/specs/    # design specs (this file lives here)
│   ├── physics/              # method derivations, unit conventions, references
│   └── adr/                  # architecture decision records
├── docker/                   # local-CPU test image + AWS runtime image
└── reference/                # read-only porting oracles
    ├── eMoScat/              # old C++ QSCAT
    └── libXcuda/             # recovered CUDA layer (for eventual GPU work)
```

- **`libs/qscat` is the standard library** — only validated methods are promoted here.
- **`projects/` is where toy models prove themselves** before promotion.
- **`reference/` is read-only** — an oracle for algorithms/math, never edited.

## Custom skills (`.claude/skills/`)

Domain- and methodology-specific; they complement (do not duplicate) Superpowers.

1. **`qm-method-lifecycle`** — flagship process skill enforcing design → toy → validate
   → optimize → promote for every new method.
2. **`numerical-validation`** — proving QM/numerical code correct where exact equality
   doesn't apply: analytic benchmarks (H atom, harmonic oscillator, Coulomb phase
   shifts), convergence/refinement studies, conservation checks (norm, unitarity of
   evolution), differential testing vs. `reference/` and vs. mpmath.
3. **`python-to-rust-kernel`** — profile → decide → scaffold a PyO3/maturin crate →
   mirror the Python API → benchmark speedup → keep Python as the differential-test
   oracle.
4. **`containerize-and-run`** — Docker patterns for local CPU test and AWS-ready
   runtime; reproducible builds via uv + maturin.
5. **`qscat-conventions`** *(reference skill)* — shared vocabulary: atomic units,
   FEM-DVR-ECS notation, tolerance conventions.

## Custom agents (`.claude/agents/`)

1. **`port-scout`** — read-only archaeologist over `reference/`. Given a method, extracts
   the math and algorithm (not the C++) for clean Python reimplementation.
2. **`physics-reviewer`** — reviews physical/numerical correctness: units, conservation,
   boundary conditions, ECS contour, convergence — what a generic reviewer misses.
3. **`rust-kernel-engineer`** — specialist for PyO3/Rust kernels with criterion
   benchmarks and differential tests against the Python oracle.

Superpowers already provides brainstorming, TDD, systematic-debugging, code-review,
writing-plans, etc.; only domain-specific gaps are added here.

## CLAUDE.md (operating manual)

Top-level manual covering: project purpose, the canonical five-stage lifecycle, the
repo map, the locked tech decisions, conventions (atomic units, test tolerances), and
when to reach for each skill/agent and the `reference/` oracles.

## Out of scope (this phase)

- Porting any physics/numerics from `eMoScat`.
- GPU/CUDA implementation.
- AWS deployment infrastructure (containers are AWS-ready, but no cloud infra yet).
- CI configuration (added when a remote host is chosen).

## Verification

The scaffold is "done" when:

- `uv sync` resolves the workspace, and `import qscat` works from a clean checkout.
- `pytest` runs green on a trivial placeholder test in `libs/qscat`.
- The Rust example kernel in `native/qscat-kernels` builds via maturin and is importable
  and differential-tested against a Python reference.
- `docker build` produces a runnable CPU image that executes the placeholder test.
- Each custom skill and agent is present, discoverable, and passes a smoke check
  (invocable / lints clean).
- `CLAUDE.md` accurately describes the above.
```
