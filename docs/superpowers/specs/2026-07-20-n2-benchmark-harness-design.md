# N₂ LCP Benchmark Harness — Design Spec

**Date:** 2026-07-20
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending

## Context

qModeling needs a concrete verification target for the electron–N₂ ²Π_g shape-resonance
problem (the Local Complex Potential / LCP model) before any scattering solver is ported.
This harness establishes that target now: it encodes the N₂ experiment definition and the
reference values, runs in the CPU Docker image, and reports which checks are already
computable versus pending a ported solver. As methods land (via `qm-method-lifecycle`),
pending checks flip to real pass/fail — the harness *is* the acceptance test those methods
must satisfy.

Two independent sources of truth anchor it:
1. **Closed-form model + literature** — the neutral N₂ Morse potential and the
   electron–molecule interaction are closed forms (extracted from `reference/eMoScat`); the
   ²Π_g resonance sits at ≈2.3–2.4 eV, Γ≈0.4–0.5 eV.
2. **Golden data** — `CSVE.V00.J00`, time-independent vibrational-excitation cross sections
   from prior work of **Karel Houfek** (not from this repo), used as regression anchors.

## Provenance & data

- **`CSVE.V00.J00`** (currently at repo root; to be moved into the benchmark):
  - Source: Karel Houfek, **time-independent** calculation. External to `eMoScat`.
  - 400 rows (energy points), 32 whitespace-separated columns, Fortran `E` notation.
  - Col 1 = collision energy in **Hartree** (5×10⁻⁴ … 0.2, step 5×10⁻⁴).
  - Col 2 = elastic / vibrationally-elastic (v=0→0); cols 3–32 = v=0→1 … v=0→30.
  - Cross sections in **atomic units (bohr²)**. Initial state v=0, J=0.
  - Higher-v channels are exactly 0 below their energetic threshold.
- **Model parameters** (from `reference/eMoScat/input/experimental/N2-model.json`, atomic
  units): `reduced_mass=12766.36`, `impulsemomentum=2` (l=2), `D_0=0.75102`,
  `alpha_0=1.1535`, `R_0=2.01943`, `lambda_inf=6.21066`, `lambda_1=1.05708`,
  `R_lambda=-27.9833`, `lambda_c=5.38022`, `R_c=2.405`, `alpha_c=0.4`.
- **Closed-form model** (extracted by port-scout; full detail in
  `.superpowers/sdd/n2-lcp-model-extraction.md`, to be promoted to `docs/physics/`):
  - `V0(R) = D_0*(exp(-2*alpha_0*(R-R_0)) - 2*exp(-alpha_0*(R-R_0)))` — Morse, min −D_0 at R_0.
  - `lambda0(R) = (lambda_c - lambda_inf)*(1 + exp(lambda_1*(R_c - R_lambda)))`;
    `lambda(R) = lambda_inf + lambda0(R)/(1 + exp(lambda_1*(R - R_lambda)))` — λ(R_c)=λ_c.
  - `V_int(r,R) = -lambda(R)*exp(-alpha_c*r**2)`; `V_eff_el(r,R)=l*(l+1)/(2 r**2)+V_int`.
  - E_res(R), Γ(R) are **not** closed form — they are the ECS complex-eigenvalue pole of the
    fixed-R electronic Hamiltonian (a small eigensolver, deferred).

## Structure

```
validation/n2/
├── model.py              # closed-form V0(R), lambda(R), V_int(r,R), V_eff_el(r,R) — pure numpy
├── config.json           # N2 parameters (copied from eMoScat N2-model.json, provenance noted)
├── data/
│   ├── CSVE.V00.J00      # Houfek golden data (moved from repo root)
│   └── MANIFEST.md       # source, units, columns, initial state
├── reference.py          # literature resonance values + selected golden-data anchors + tolerances
├── experiment.py         # harness entrypoint: PASS/PENDING/FAIL table, nonzero exit on FAIL
├── test_n2.py            # pytest wrapper (runs under `uv run pytest` and docker test stage)
└── README.md             # physics, references, provenance, how to run
```

`model.py` stays local to the benchmark for now (YAGNI); it may be promoted into
`qscat` later once a method depends on it.

## Checks

`experiment.py` groups checks and prints one row each (PASS / PENDING / FAIL):

**Group A — closed-form model (green now):**
- A1 `V0(R_0) == -D_0` (rtol 1e-12).
- A2 argmin of V0 over a fine R grid ≈ R_0 (within grid spacing).
- A3 `V0(R→∞) → 0` (small at R=20 a₀).
- A4 `lambda(R_c) == lambda_c` (rtol 1e-12).
- A5 `V_int(r,R)` is a negative, decaying-in-r well; `V_eff_el` includes the l=2 centrifugal term.

**Group B — resonance position (PENDING, needs small ECS eigensolver):**
- B1 `E_res(R_0) ≈ 2.3–2.4 eV`, `Γ(R_0) ≈ 0.4–0.5 eV`. Reported PENDING until the ECS
  eigensolver is ported. (Literature + port-scout prototype ≈2.44 eV / 0.46 eV.)

**Group C — Houfek golden data:**
- C1–C4 **data integrity (green now):** file parses to 400×32; energy strictly increasing
  over 5e-4…0.2 Ha; all cross sections ≥ 0; per-channel thresholds are ordered (channel
  v=0→(j+1) opens at higher energy than v=0→j; each is 0 below its threshold).
- C5 **value anchors (PENDING time-independent solver):** a small fixed set of
  `(energy, channel)` coordinates (e.g. elastic and v=0→1,2,3 near E=0.2 Ha, plus one
  mid-range point and one near-threshold point). Reference values are **looked up from
  `CSVE.V00.J00`** (never hardcoded). When the time-independent solver lands, its output at
  those coordinates is compared within `rtol` (default 5%, tunable). Reported PENDING now.

**Group D — time-dependent model (PENDING, later):** same anchor comparison once the
time-dependent LCP propagation is ported; only if/when TD golden data or TI/TD agreement is
available.

`experiment.py` exit code: **0** if no FAIL (PENDING does not fail); **non-zero** if any
integrity/closed-form check FAILs (e.g. data file missing/corrupt, model formula regressed).

## Running

- Locally: `uv run python validation/n2/experiment.py`; `uv run pytest validation/n2`.
- CPU Docker (image already built): `docker run --rm qmodeling:runtime python validation/n2/experiment.py`.
- Convenience wrapper: `docker/run-n2.sh` (builds/uses the runtime image and runs the experiment).
- `test_n2.py` runs the green checks under the normal suite and the docker `test` stage.

## Out of scope (this phase)

- The ECS resonance eigensolver (Group B stays PENDING).
- The time-independent and time-dependent LCP solvers (Groups C5/D stay PENDING).
- Promoting `model.py` into `qscat`.
- Any change to the existing scaffold beyond adding `validation/n2/` and `docker/run-n2.sh`.

## Verification

- `uv run pytest validation/n2` passes (green groups A + C1–C4).
- `uv run python validation/n2/experiment.py` prints the table and exits 0, with Group A + C
  integrity as PASS and Groups B, C5, D as PENDING.
- `docker run --rm qmodeling:runtime python validation/n2/experiment.py` produces the same
  table and exit 0 inside the CPU image.
- `CSVE.V00.J00` no longer at repo root; lives under `validation/n2/data/` with a manifest.
