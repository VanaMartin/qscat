# Validation harnesses

Every method in qscat is checked against an analytic benchmark, a conservation
law, a convergence study, or an independent reference. This note is the index
of those checks: what each harness in `validation/` gates, against which
oracle, at what tolerance, and the command that runs it.

Atomic units throughout.

## The harnesses

| Harness | Oracle | Gate | Run it |
|---|---|---|---|
| `validation/n2` | Houfek's independent `CSVE.V00.J00` data | `GATED_RTOL = 1e-3` for the exact 2-D solver (`validation/n2/exact2d.py:82`); the LCP (Group C5) is held only to the cross-model `ANCHOR_FACTOR = 3.0` band (`validation/n2/reference.py:41`); the TD-CN route (Group D) gates on a conjunction — `abs(ratio_td_ti - 1.0) <= 0.10` AND the same `ANCHOR_FACTOR = 3.0` band against Houfek (`validation/n2/td_check.py:93-95`) | `uv run python -m validation.n2.experiment` |
| `validation/h2plus` | the published `omega_i^j` level table (`reference_levels.LEVELS`) | overlap verdicts against a Born-Oppenheimer basis (`qscat.core.assignment.pair_by_overlap`) | `uv run python -m validation.h2plus.bo_overlap` |
| `validation/diatomic` | none published — the exact 2-D solver **is** the oracle | well-posedness (finite, non-negative, right shape) plus a per-molecule scanned resonance floor on the fine `da_grid` (F2: sigma_v'=1 >= 0.15 bohr^2 over 0.010-0.044 Ha; NO: >= 10.0 bohr^2 over 0.017-0.045 Ha — `validation/diatomic/test_diatomic.py`) | `uv run pytest validation/diatomic -q` |
| `validation/tuning` | the eMoScat per-molecule decks | `_2D_SPOT_CHECK_RTOL = 0.02`, the once- vs twice-refined 2-D sigma_DA agreement (`validation/tuning/test_emoscat_decks.py:257`); `_2D_CONVERGENCE_RTOL = 0.15`, the base resonance-aware grid vs once-refined agreement (`validation/tuning/test_resonance_aware.py:107`) | `uv run python -m validation.tuning.calibrate` |

Each constant above is transcribed from the source file and line cited beside
it; re-read those lines before relying on a value here, since a gate can be
retuned without this table noticing.

## What each one actually establishes

### `validation/n2`

This is the only harness in the repo with an external, independently
computed oracle: Karel Houfek's `CSVE.V00.J00` cross-section table, produced
by a different implementation of the same model and method. Its Group E1
(`validation/n2/exact2d.py`) holds the exact 2-D driven-Lippmann-Schwinger
solver to `GATED_RTOL = 1e-3` against Houfek at the four anchors clear of
their own vibrational threshold — a differential-oracle bound, tight because
exact-2D and Houfek are the *same* physics computed two independent ways.
Group C5 (the older 1-D local-complex-potential route,
`validation/n2/cross_section.py`) is held only to the much looser cross-model
`ANCHOR_FACTOR = 3.0` band, because the LCP is a *known*, energy-dependent
approximation to the exact model, not a candidate for exact agreement — a
large, well-characterized LCP discrepancy there is a finding, not a bug.
Group D (the time-dependent Crank-Nicolson route, `validation/n2/td_check.py`)
is held to more than that: its gate is a conjunction of `abs(ratio_td_ti -
1.0) <= 0.10` — TD-vs-TI agreement, the field comment in `td_check.py` calls
this "the real gate" — AND the same `ANCHOR_FACTOR = 3.0` band against
Houfek (`td_check.py:56-57, 93-95`). The 10% TI bound is the one that
actually constrains the TD solver day to day; the Houfek band is inherited
from C5 mainly because the TD route shares the LCP's `V_d(R)`/`Gamma(R)`
machinery. What this harness does **not** establish: that N2 the molecule is
correctly modeled. The model potential is a given testbed; agreement with
Houfek certifies that this solver correctly solves *that* model, not that the
model is physically realistic.

### `validation/h2plus`

This harness answers a question energy proximity alone cannot: whether an
angle-stable pole found by `qscat.core.exact_resonance_states` is a genuine
resonance at all. ECS two-angle stability is necessary but not sufficient —
on the H2+ campaign four of the angle-stable poles (window 0) scored overlaps
of 6e-4 to 7e-3 against the Born-Oppenheimer product basis, where genuine
states score 0.87-0.99 (`libs/qscat/qscat/core/assignment.py:25`,
`libs/qscat/qscat/core/resonance.py:86`). `bo_overlap.py` builds that basis
(`qscat.core.bo.bo_basis`) and pairs every pole to it by c-product overlap
(`qscat.core.assignment.pair_by_overlap`), reporting one of seven verdicts
per pole (including `box-limited`, caught only by the separate `real_weight`
check — see `docs/physics/h2plus-resonance-states.md`). What it establishes
is a verdict on *individual poles*, not a cross-section curve: it is a
position anchor, not a `sigma_DR(E)` anchor, because the resonances are
narrower than the published sweep's own energy sampling.

### `validation/diatomic`

**No independent golden data exists for NO or F2** — only N2 has Houfek's
table. For these two molecules the exact 2-D solver *is* the oracle, so
agreement of its output with itself (a converged grid vs a finer one, or a
physically-expected resonance peak appearing where the model's own
Born-Oppenheimer curve predicts one) is self-consistency, not external
validation. `test_diatomic.py` gates two things instead: well-posedness (the
computed sigma is real, finite, non-negative, and the right shape on a small
grid) and a physical sanity floor — each molecule's known low-lying
resonance must drive a v'=1 excitation cross section above a measured floor
somewhere in a scanned energy window on the converged fine `da_grid` (F2:
>= 0.15 bohr^2 over 0.010-0.044 Ha at step 0.004; NO: >= 10.0 bohr^2 over
0.017-0.045 Ha at step 0.002 — NO's structure is sharp enough that the
window is scanned rather than sampled at a few points, since a single 0.002
Ha step swings the value by a factor of 5.7). This is real, useful protection
against a broken port — but it is not the same claim `validation/n2` can
make, and the note that motivated this harness (`config.py`) says so plainly
in its own docstring.

### `validation/tuning`

This harness calibrates and gates `qscat.tuning`, the automatic FEM-DVR-ECS
discretisation tuner, against the eMoScat reference decks
(`reference/eMoScat/input/{N2,NO,F2}/grids.txt`, transcribed into
`validation/diatomic/config.py`). `calibrate.py` sweeps the mesh's
de-Broglie phase constant `C` and picks the smallest value at which
`propose_grid`'s F2 nuclear grid reproduces-or-beats the eMoScat F2
dissociative-attachment deck — F2 is the deciding case because it is the
only molecule with a genuinely open DA channel in its tested range.
`test_emoscat_decks.py`'s `@slow` 2-D spot-check found a genuine gap the
1-D convergence probes could not see: the reproduce-and-beat grid passes
both 1-D probes but gives an unconverged sigma_DA (one nuclear h-refinement
moves it ~5x). `test_resonance_aware.py` closes that gap with a
resonance-aware nuclear mesh, gated at `_2D_CONVERGENCE_RTOL = 0.15` against
one further refinement. What this harness does not establish for N2/NO: the
1-D channel-representation floor it also checks is a deliberately
conservative bound that not even N2's or NO's own committed eMoScat decks
clear at `rtol=1e-3` — reported as a finding, not silently loosened.

## Where the `@slow` boundary falls

Groups within `validation/n2/experiment.py` share a per-group time budget of
roughly 60 seconds. A full time-dependent 2-D propagation at
`TD_WORKING_GRID` costs `~210-250s` (`validation/n2/experiment.py:189`) —
even the *shortest* configuration on the sub-project's own T-scan (T=600)
costs `~85s`, over budget, and is also the least-converged point that scan
measured (sigma_TD/sigma_TI = 0.760 there vs 0.931 at the converged T=1500).
So Group F1 (`validation/n2/td_exact2d.py`) does not run a TD propagation at
all: it reports the already-validated `sigma_TD` as a cited, literal
constant, recomputing only the cheap `sigma_TI` side live every harness run.
The genuine PASS/FAIL gate on the TD-vs-TI agreement this NOTE row reports
lives outside the harness, in
`projects/n2_2d_td_cross_section/test_td_cross_section.py`'s
`@pytest.mark.slow` tests, run explicitly with
`uv run pytest projects/n2_2d_td_cross_section -m slow` rather than as part
of a default harness pass. The same pattern — a cheap, always-live check in
the harness paired with an expensive, opt-in `@slow` pytest gate elsewhere —
recurs in `validation/tuning`, whose 2-D spot-checks and convergence tests
(`test_emoscat_decks.py`, `test_resonance_aware.py`) are both
`@pytest.mark.slow` and excluded from a default `not slow` pytest run.

## Reading a NOTE row

`validation/n2/experiment.py` prints four status values: `PASS`, `FAIL`,
`PENDING` (no result computed yet), and `NOTE`. A `NOTE` row is not a
skipped check — it is "a result exists but is a documented, non-gating
observation" (the module's own status-value comment), and it is never
counted toward `FAIL` or the process exit code. Two examples show the
distinction is doing real work, not hiding a failure: Group C5's elastic and
near-threshold anchors are NOTE because the LCP model is *known* to miss the
elastic background and to get the threshold law wrong there — reporting the
ratio anyway documents the size of that known gap without failing the
harness on a limitation the model was never expected to clear. Group F1's
rows are NOTE for a different reason — not a known model limitation, but a
measured execution-cost decision (the `~210-250s` figure above) to report a
result computed elsewhere rather than pay for it again in-harness. Both are
legitimate uses of NOTE: a documented, reported fact carried in the output,
never a silent omission.
