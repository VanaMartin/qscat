# CLAUDE.md — qModeling operating manual

## Purpose

qModeling is a Python-first quantum-mechanics (QM) research monorepo, home to
**QSCAT**. It resuscitates and extends prior C++/CUDA work (`eMoScat`,
`libXcuda`) as a maintainable, CPU-first Python codebase with Rust for the hot
paths that actually need it. GPU/CUDA and AWS deployment are recovered as
reference material but deliberately deferred — not part of this phase.

## The lifecycle

Every method or numerical capability moves through five stages:

1. **Design** — work out the math/algorithm (see `docs/physics/`), decide the
   API shape.
2. **Toy model** — a small, readable Python implementation under `projects/`.
3. **Validate** — check it against analytic benchmarks, conservation laws,
   convergence studies, or a reference implementation (`validation/`).
4. **Optimize in Rust** — once validated and proven a hot path, move it to a
   PyO3/maturin kernel under `native/`, keeping the Python version as the
   differential oracle.
5. **Promote to `qscat`** — validated, reusable code graduates into
   `libs/qscat`, the standard library other projects depend on.

**Invariants that hold at every stage:** the code is CPU-runnable on a laptop
(no GPU required), and it is containerizable (builds and runs under
`docker/`). See the `qm-method-lifecycle` skill for the enforced workflow.

## Repo map

```
libs/       qscat — the standard library: validated, reusable QM code
            (units, linalg, dvr, ecs, special, evolution, core, model
            submodules)
            - qscat.linalg: dimension-general sparse linear algebra --
              `kron_sum` (Kronecker sum over arbitrary D), `SparseLU` (cached
              factorization with fill-in/memory diagnostics), and `c_product`
              (the bilinear, non-conjugated ECS inner product) -- see
              docs/physics/nd-tensor-hamiltonian.md. `SparseLU` dispatches
              between SuperLU (`backend="scipy"`, numpy/scipy-only, the fallback
              AND the differential oracle) and a complex-symmetric MUMPS `SYM=2`
              backend (`backend="mumps"`, the optional `qscat[mumps]` extra):
              `backend="auto"` (default) picks MUMPS if present else SuperLU,
              and `set_default_backend`/`default_backend()` override what
              `"auto"` resolves to process-wide. On the ECS-complex-symmetric
              N₂ matrices MUMPS beats SuperLU by 72× in factor time / 9× in peak
              RSS at the 143k production deck (3.6 s / 0.8 GB vs 260 s / 7.4 GB)
              -- see docs/physics/mumps-sparse-backend.md. `SparseLU.refactor(
              A_new)` reuses the symbolic analysis (same sparsity pattern, e.g. a
              diagonal shift `E_tot·I − H` across an energy sweep): MUMPS
              `factor(reuse_analysis=True)` skips the SCOTCH ordering; scipy
              re-runs `splu` (correct, no reuse); a pattern guard raises on a
              structure mismatch. On the N₂ working grid a reuse sweep is ~80% /
              ~5× faster than fresh-per-energy (the analysis dominates the cheap
              numeric factor there; the fraction shrinks for larger decks) — see
              docs/physics/ti-energy-sweep-reuse.md.
            - qscat.dvr: FEM-DVR-ECS radial grid (`FemDvrEcsGrid`), kinetic-
              energy assembly (`kinetic`), and diagonal-potential Hamiltonian
              + eigensolver helpers (`hamiltonian`, `eigen`) — see
              docs/physics/femdvr-ecs.md.
              Also `kinetic_sparse` and the N-dimensional tensor layer
              (`TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd`):
              H = sum_d I x .. T_d .. x I + diag(V) for any D, sparse (CSR),
              validated on analytic box/oscillator benchmarks at D = 1, 2, 3.
            - qscat.ecs: exterior-complex-scaling coordinate map (`ecs_map`),
              the single source of the `z(x) = x` / `R0 + (x-R0)e^{i theta}`
              transform used by qscat.dvr's complex tail; also
              `find_resonance_pole(eigs_a, eigs_b, window)`, the general
              two-spectrum resonance-pole matcher (promoted from the N2
              resonance project) — see docs/physics/n2-resonance.md.
            - qscat.evolution: `make_cn_stepper(H, dt)`, a general Crank-
              Nicolson time propagator for `d/dt psi = -i H psi` (any complex,
              possibly non-Hermitian `H`) — promoted from the N2
              time-dependent cross-section project — see
              docs/physics/n2-td-cross-section.md. Also
              `make_sparse_cn_stepper(H, dt)`, the sparse sibling (factors
              once with `SparseLU`, matches the dense stepper to round-off) —
              promoted from the N2 2-D time-dependent cross-section project —
              see docs/physics/n2-2d-td-cross-section.md. Also
              `make_pade_stepper(H, dt, order)` (+ `pade_roots`), the order-N
              diagonal-Padé generalization of the sparse CN stepper (order 1 ==
              CN; `O(dt^(2N+1))` per step). Order-1 CN under-converges badly
              over a long propagation (~100% accumulated error at dt=0.5-1.0 vs
              `expm`); order-3 Padé is what makes the TD cross section converge
              to the TI oracle to ~1-2% — see docs/physics/n2-2d-td-cross-section.md.
            - qscat.core: the model-INDEPENDENT electron–diatomic VE-scattering
              engine, promoted from the N2 sub-projects. `driven.ve_cross_section`
              (exact TI driven Lippmann-Schwinger, `SparseLU.refactor` sweep,
              σ=π|S−δ|²/2E) and `time_dependent.td_ve_cross_section` (order-3
              Padé propagation + Tannor-Weeks transform + elastic free-reference
              subtraction) — both take a `model`. Also
              `dissociation.da_cross_section` (exact TI **dissociative
              attachment**: the same driven Ψ₊ solve projected onto the nuclear
              dissociation channel with the rearrangement interaction
              `V_DR = V_int + v0 − V_int(r,R→∞)`, σ_DA=4π³|T|²/2E; plus
              `anion_electronic_states`, `v_dr_diag`) — the DA magnitude needs a
              per-molecule NUCLEAR grid (fast K_R~58 exit wave), built by
              `grids.segmented_grid` from eMoScat's per-molecule decks — see
              docs/physics/diatomic-ve-cross-sections.md. Also `lcp` — the
              **LCP (local-complex-potential) approximation** OF the DA (the
              approximation under test vs the exact `da_cross_section` oracle):
              `local_complex_potential` (model-independent R-dependent
              `V_d(R)/Γ(R)` via `qscat.ecs.find_resonance_pole`, seeded from the
              anion state) + `lcp_da_cross_section` (1-D TI resolvent, boundary-
              wavefunction-VALUE flux `ψ(X)=ψ_coeff[b]/√w_b`, on the fine
              per-molecule nuclear grid). F₂ σ_DA agrees with the exact-2D to
              ~11% away from threshold; documented departures near threshold
              (LCP misses the exact's spike) and in VE-elastic (LCP misses the
              non-resonant background) — see
              docs/physics/diatomic-ve-cross-sections.md. Plus `channels`, `grids`
              (parameterized FEM-DVR-ECS builders + `segmented_grid` for
              eMoScat's `(n_elem, endpoint)` deck format), `vibrational` (`v0`
              passed in), `wavepacket`, `correlation`, `plot`. **`qscat.core` never
              imports `qscat.model`/`projects` at runtime** (depends only on the
              `ResonanceModel` protocol; enforced by
              `test_core_no_model_import.py`) — see
              docs/physics/qscat-core-scattering.md.
            - qscat.model: everything tied to a specific model — the
              `ResonanceModel` protocol (the contract `qscat.core` depends on),
              `DiatomicResonanceModel` (the shared Morse+sigmoid+Gaussian form),
              and the `N2`/`NO`/`F2` parameter registry. `qscat.model.N2` is the
              single source of truth for the N2 model (the N2 projects consume it
              via thin shims). Adding a molecule = a registry entry + validation,
              never solver code — see docs/physics/qscat-core-scattering.md.
native/     Rust kernels (qscat-kernels crate) built with PyO3/maturin,
            mirroring validated Python APIs for hot paths
projects/   per-problem research and toy models — lifecycle stages 1-2
            - `n2_ti_cross_section`: time-independent (resolvent/driven-
              equation) N₂ vibrational-excitation cross-section solver
              (`nuclear_grid.py`/`vibrational.py`/`vres.py`/`cross_section.py`),
              built on `qscat.dvr`/`qscat.ecs` and the N₂ resonance pole
              finder — see docs/physics/n2-cross-section.md.
            - `n2_td_cross_section`: time-dependent (Crank-Nicolson
              propagation + energy transform) route to the same N₂
              vibrational-excitation cross section (`propagator.py` — thin
              re-export of `qscat.evolution.make_cn_stepper`;
              `td_cross_section.py` — doorway wavepacket propagation under
              `H_res`, correlation function, energy transform), validated
              against the TI solver as an exact differential oracle — see
              docs/physics/n2-td-cross-section.md.
            - `n2_2d_cross_section`: the exact 2-D (electronic r × nuclear R)
              driven Lippmann-Schwinger solver for the same N₂
              vibrational-excitation cross section — no local-complex-
              potential reduction (`electronic_grid.py`/`channels.py`/
              `hamiltonian2d.py`/`cross_section_2d.py`/`convergence.py`/
              `nuclear_density.py`), validated standalone (free-particle and
              first-Born limits, S-matrix reciprocity/unitarity) and then
              gated against Houfek's data as an independent implementation of
              the same model/method; once gated, it is the ORACLE the
              1-D LCP solver is compared against — see
              docs/physics/n2-2d-cross-section.md. `ve_cross_section_2d` sweeps
              energies analyze-once/refactor-per-energy (via
              `SparseLU.refactor`) — same σ, cheaper sweep. The generic,
              experiment-agnostic `cross_section_plot.plot_cross_sections`
              (no physics, reference passed as an argument — reusable for
              F₂/NO) renders the dense σ_{0→v'}(E) curves; the committed N₂
              curve vs Houfek is docs/physics/figures/n2-2d-ti-cross-section.png
              — see docs/physics/ti-energy-sweep-reuse.md.
            - `n2_2d_td_cross_section`: time-dependent (sparse Crank-Nicolson
              propagation + Tannor-Weeks energy transform) route to the SAME
              exact 2-D N₂ vibrational-excitation cross section as
              `n2_2d_cross_section` — an incident Gaussian wavepacket
              `g(r) chi_0(R)` propagated under `H_2D`
              (`wavepacket.py`/`td_propagation.py`/`correlation.py`/
              `td_cross_section.py`/`convergence.py`/`observation.py`),
              validated against the exact 2-D solver as an exact differential
              oracle (σ_TD/σ_TI = 0.93 at E=0.10, 1.10 at E=0.15). The
              ELASTIC (v'=v_init) channel subtracts a free-particle (V_int=0)
              reference S_free(E) instead of a literal 1 — the transform's
              outgoing normalization makes S_free≈2π²≠1, so |S−1|² left a ~500×
              spurious elastic background; `td_ve_cross_section_2d(...,
              subtract_free_reference=True)` (default) runs the reference
              propagation. The propagator is the order-3 diagonal Padé
              (`make_pade_stepper`, dt=1.0, eMoScat's setting) — order-1 CN
              under-converged, capping TD-vs-TI at ~10-15%; with order-3 the
              TD cross section matches the exact TI (and Houfek) to ~1-2%
              median across 0.04-0.18 Ha for elastic + first excitations,
              boomerang oscillations resolved point-by-point — see
              docs/physics/n2-2d-td-cross-section.md and the
              td-elastic-wavepacket-normalization note.
validation/ analytic benchmarks, golden datasets, convergence studies
            - `validation/n2/`: N₂ electron-scattering harness; its C5 group
              anchors this solver's σ_{0→v'}(E) against Karel Houfek's
              independent `CSVE.V00.J00` data (documented cross-model
              tolerance, not exact agreement); its D1 group cross-checks the
              TD solver against both the TI solver (rtol<=0.10) and the same
              Houfek data at the GATED anchors (`td_check.py`); its E1 group
              (`exact2d.py`) reruns the same 6 anchors through the exact 2-D
              solver, GATED at `GATED_RTOL=1e-3` against Houfek (a tight
              differential-oracle bound, not the LCP's cross-model factor-3
              band) with the two DOCUMENTED-LIMITED LCP anchors reported as
              NOTE rows showing how the exact model closes that gap; its F1
              group (`td_exact2d.py`) reports the time-dependent 2-D solver's
              σ_TD-vs-σ_TI agreement at the two validated anchors as
              recorded, cited NOTE rows rather than a live in-harness
              propagation — a full run costs ~210–250s (measured), far over
              the harness's per-group budget, so the genuine PASS/FAIL gate
              lives in `n2_2d_td_cross_section/test_td_cross_section.py`'s
              `@slow` tests instead. Run the harness with
              `uv run python -m validation.n2.experiment`. `ti_curve.py` is the
              N₂-specific driver behind the dense exact-2D σ(E) figure: it reads
              Houfek (`loader`) and calls the generic `plot_cross_sections`
              (validation may import projects; projects must not import
              validation), and `test_ti_curve.py` gates the dense curve against
              Houfek at the anchors — see docs/physics/ti-energy-sweep-reuse.md.
            - `validation/diatomic/`: the NO and F₂ exact-2D VE cross sections —
              the model port, computed entirely through `qscat.core` +
              `qscat.model` (the first consumers of sub-project A beyond N₂).
              `config.py` holds per-molecule grid/energy config (incl. the
              eMoScat per-molecule NUCLEAR deck + `MoleculeConfig.da_grid()` used
              for DA); `curves.py` (`compute_ti_curve`, `main`) computes the
              exact-2D TI VE σ(E) oracle + figures (`docs/physics/figures/{no,f2}-
              2d-ti-cross-section.png`); `da_curves.py` (`compute_da_curve`,
              `main`) computes the exact-2D TI σ_DA(E) oracle on the eMoScat
              per-molecule nuclear grids + figures (`docs/physics/figures/{f2,no}-
              2d-ti-da-cross-section.png`); `lcp_da_curves.py` (`compute_lcp_da_curve`,
              `main`) computes the LCP-approximation σ_DA(E) and overlays it on the
              exact-2D oracle (`docs/physics/figures/{f2,no}-2d-da-lcp-vs-exact.png`)
              — the sub-project-B deliverable (LCP ~11% of exact away from
              threshold; documented near-threshold departure) + `MoleculeConfig`'s
              `lcp_elec_grids()`/`lcp_nuclear_grid()`. No independent golden data
              exists for NO/F₂ (only N₂ has Houfek's), so the exact solver IS the oracle;
              NO/F₂ have lower/sharper resonances than N₂, so their TD-vs-TI
              check is a documented follow-on — see
              docs/physics/diatomic-ve-cross-sections.md.

`projects/` and `validation/` (and their sub-project directories) are real
Python packages (`__init__.py` present at every level); all intra-repo
imports are package-absolute (e.g. `from projects.n2_resonance.potential
import v0`, `from validation.n2 import loader`) — no bare intra-dir imports,
`sys.path.insert`, or `importlib.util.spec_from_file_location` workarounds.
Root `pyproject.toml` sets `pythonpath = ["."]` under
`[tool.pytest.ini_options]` so pytest resolves these packages from the repo
root without any path hacking. A module meant to be run directly (not just
imported), such as `validation/n2/experiment.py`, is invoked as
`python -m validation.n2.experiment`, not by file path.
benchmarks/ standalone measurement scripts (a real package, run via
            `python -m benchmarks.<name>`; imports `projects`, never imported by
            `projects`/`qscat`) -- `mumps_vs_superlu` measures the MUMPS vs
            SuperLU factor/solve/RSS/fill on the real N₂ 2-D matrices;
            `sweep_reuse` measures the energy-sweep symbolic-reuse speedup
            (`SparseLU.refactor` vs fresh-per-energy); both run in the Docker
            `test` image (need system MUMPS + `qscat[mumps]`).
reference/  read-only oracles: eMoScat (C++/CUDA snapshot), libXcuda
            (CUDA submodule) — for porting reference only, never imported
docs/       specs/plans (docs/superpowers), physics notes (docs/physics),
            and ADRs (docs/adr)
docker/     layered CPU images: base (architecture/vendor) + app (build/
            test/runtime)
.claude/    skills and agents (below) plus repo permissions settings
```

## Tech decisions

- **Python >=3.12**, managed with `uv`; the exact interpreter is pinned in
  `.python-version` (`3.12`) for reproducibility and Docker parity.
- **uv workspace**: `pyproject.toml` declares `members = ["libs/*", "native/*"]`
  with `package = false` at the root. The canonical setup command is
  `uv sync --all-packages` — a plain `uv sync` prunes the `qscat` and
  `qscat_kernels` workspace members.
- **Rust + PyO3/maturin** for compiled kernels (`native/qscat-kernels`),
  built and installed into the active environment via
  `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml`.
- **Testing**: pytest, hypothesis for property-based tests, pytest-benchmark
  for kernel benchmarks; every Rust kernel gets a Python differential test
  against its oracle.
- **Quality tooling**: ruff (lint), mypy --strict (types), cargo clippy (Rust
  lint). Pre-commit (`.pre-commit-config.yaml`) currently enforces only
  **ruff**, **ruff-format**, and **cargo fmt --check**; mypy and clippy are
  run manually / in CI, not (yet) pre-commit hooks.
- **Units**: atomic units throughout (`libs/qscat/qscat/units.py`); no ad hoc
  unit conversions scattered through method code.
- **Docker is layered, not a single Dockerfile**:
  - `docker/base.Dockerfile` builds `qmodeling-base` — the architecture /
    BLAS-FFT vendor layer: OpenBLAS + LAPACK(E) + FFTW3 dev libs, pkg-config,
    a Rust toolchain, and uv/Python 3.12. Code targets the standard CBLAS /
    LAPACKE / FFTW3 ABIs, so this layer is swappable. Default vendor choice
    (OpenBLAS/FFTW3) keeps ARM/Graviton open; an MKL x86-64 variant is a
    planned future alternative, not implemented yet. It also provisions system
    **MUMPS** (`libmumps-seq-dev` + `libscotch-dev`, the sequential build) for
    `qscat.linalg`'s optional MUMPS backend, and **synthesizes the pkg-config
    `.pc` files** Debian omits (the conda-forge `{d,z,c,s}mumps_seq` names
    python-mumps looks for) so the extra builds against the system library.
  - `docker/Dockerfile` is `FROM ${BASE_IMAGE}` and layers `build` → `test` →
    `runtime` on top, using `uv sync --all-packages` for setup. The `test`
    stage adds `--extra mumps` so the MUMPS backend is exercised; `runtime`
    deliberately omits it, keeping python-mumps out of the production image.
    MUMPS tests run in the container and `@skipif`-absent on a MUMPS-less box
    (a bare Mac), so the Mac suite stays green while the same tests run + pass
    in Docker — see docs/physics/mumps-sparse-backend.md.
  - `docker/build.sh [test|runtime]` builds the base image then the
    requested app target. Verified working: `test` prints `5 passed`;
    `runtime` prints `qscat 0.0.0 ready`.
- **GPU/CUDA and AWS deployment are deferred.** `reference/libXcuda` is kept
  as future-GPU-kernel reference only.

## Skills & agents

| Name | Kind | Use when |
|---|---|---|
| `qm-method-lifecycle` | skill | Adding or porting any QM method/capability — enforces design → toy → validate → optimize-in-Rust → promote-to-qscat. |
| `numerical-validation` | skill | Validating quantum/numerical code where exact equality doesn't apply: analytic benchmarks, convergence studies, conservation checks, differential testing. |
| `python-to-rust-kernel` | skill | A validated Python method has a proven hot path — scaffolding a PyO3/maturin crate under `native/`, mirroring the API, benchmarking, keeping Python as the oracle. |
| `containerize-and-run` | skill | Packaging a capability to run in Docker locally or prep it for AWS — reproducible multi-stage `uv` + `maturin` builds, CPU-only. |
| `qscat-conventions` | skill | Unsure how the project names or measures things — atomic units, FEM-DVR-ECS notation, tolerance defaults, standard-library layout. |
| `port-scout` | agent | Before porting anything from `reference/eMoScat` or `reference/libXcuda` — read-only archaeologist that extracts the math/algorithm, not the C++. |
| `physics-reviewer` | agent | Before promoting a method into `qscat` — reviews for physical/numerical correctness (units, conservation, boundary conditions, ECS handling, convergence), not style. |
| `rust-kernel-engineer` | agent | During the optimize-in-Rust stage — builds PyO3/Rust kernels in `native/` mirroring a validated Python API, with benchmarks and differential tests. |

For general engineering workflow — brainstorming, TDD, systematic debugging,
writing/executing plans, requesting code review — use the **Superpowers**
skills (`superpowers:brainstorming`, `superpowers:test-driven-development`,
`superpowers:systematic-debugging`, `superpowers:writing-plans`,
`superpowers:executing-plans`, `superpowers:requesting-code-review`, etc.).
These are general-purpose and complement the qModeling-specific skills above.

## Reference oracles

`reference/` holds read-only prior codebases: `reference/eMoScat` (C++/CUDA
FEM-DVR-ECS electron-molecule scattering code, snapshot copy) and
`reference/libXcuda` (recovered CUDA layer, git submodule). Never edit,
build, or import these as dependencies — they exist only to be read. Always
use the `port-scout` agent before porting a method out of them; it extracts
the underlying math/algorithm so it can be reimplemented cleanly in Python.

## Common commands

```bash
# Setup — installs qscat and builds the Rust qscat_kernels in one step.
# Plain `uv sync` is NOT sufficient: it prunes the qscat/qscat_kernels members.
uv sync --all-packages

# Run the full test suite
uv run pytest

# Right after a maturin build in the same step, skip the (already-satisfied)
# sync to avoid a redundant resolve:
uv run --no-sync pytest

# Rebuild the Rust kernel extension in place during kernel iteration
uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml

# Lint / type-check
uv run ruff check .
# mypy is type-clean over the qscat library (type stubs for the Rust
# qscat_kernels extension are pending, so repo-wide strict mypy isn't
# claimed to pass yet):
uv run mypy libs/qscat

# Build and test the CPU Docker images (base, then test target)
docker/build.sh test
```
