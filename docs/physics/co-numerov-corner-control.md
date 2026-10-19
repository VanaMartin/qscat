# Cross-anchor energy corners and endpoint control

## Decision and finite scope

The [two-point increment control](../../projects/ukrmol_co/NUMEROV_INCREMENT_CONTROL.md)
passes its declared gates, but its coarser-step 800-bohr/native64 K discrepancies
remain 1.59e-6 / 1.46e-6. Before broader use, test whether the same algebra works
at the energy extremes of every saved sector, and whether the endpoint growth
at the two large-K points is consistent with fixed-step truncation.

This is a finite Python validation, with no change to either propagator,
boundary decoder, tail, native reference or physical inputs. Rydberg coefficients
and bohr lengths follow the existing real-axis problem; `surface_matrix` receives
total Hartree energy. Unsymmetrized K and separate absolute K/S/phase/invariant
budgets of 1e-6 are retained. Historical false/unset verdicts and original exits
remain independent records. No native or inner-region executable is invoked.

## Analytic prerequisite

Use the same binary64 scalar near-K=-185 and rotated two-open/one-closed
boundaries as the previous analytic control. At energy 0.0735 Ryd and a=18 bohr,
evaluate the increment formulation at the Cartesian grid
`h=[0.0125,0.00625]`, `b=[200,400,800]`, order 8: exactly twelve solves.
Independent continuous boundary matching and the scalar closed-form discrete
Numerov oracle are evaluated at 80 and 110 digits. Require:

- Oracle precision drift <1e-20 at every configuration.
- Scalar candidate/discrete K error <=1e-7 at every configuration.
- Finest-step continuous K/S/phase error <=5e-7 at every endpoint.
- All conservation checks <=1e-6.

The coarser-step continuous errors are **reported without a pass gate**:
the hypothesis predicts length-dependent truncation there, not that all those
rows are already accurate. No measured result is used to choose the endpoint
or grid. Failure/unset of this prerequisite prevents submission of CO controls.
Budget: four reserved physical cores, 4 GiB, expected 20 s / 120-s hard cap.

## CO contract

After analytic qualification, freeze its hash with the numerical code and input
hashes before submission. Twelve corner identities are indices 0 and 98
(requested 0.1 and 5 eV) for all compressed/equilibrium/stretched B1/B2 sectors.
Their six configurations reuse the original passing probe spacings:

```
(h,b,order) = (.05,400,8), (.025,400,8), (.0125,400,8),
              (.0125,200,8), (.0125,800,8), (.0125,400,6)
```

Additionally retain the two stretched index-18 identities. Reuse their three
existing order-8/400-bohr step rows at h=.025/.0125/.00625, with their existing
native64 references. Add only three solves per identity: h=.00625 at b=200/800,
order8, and at b=400, order6. This is a separate, declared 2-D spacing/endpoint
control, not an automatic spacing ladder. Exactly **78 new CO solves** are
allowed (72 corners + 6 endpoint/tail controls); six rows are reused unchanged.
No intermediate energy or resonance curve is inferred from these corners.

Native32 references for the corners retain exact saved 99-point energy identity,
complete binary/header/EOF validation and their native16→32 controls. Stretched
index18 uses the exact single-energy native64 reference qualified previously.
Native refinement is conditional reference evidence, not an absolute error bound.

Each identity requires complete six-configuration coverage. Checks are final
step, endpoint, order and finest/native differences **separately** for K/S/phase;
unitarity, reciprocity; and the prior S-step trend criterion (eightfold reduction
or 1e-8 floor, with both raw changes reported). Also require **every finest-step
endpoint row** (200/400/800) to meet the native K/S/phase budget. This explicitly
closes the prior all-row gap at the final spacing; coarse step rows retain their
diagnostics and are not silently relabelled qualified. No Richardson bound.

Submit the two endpoint identities first as fixed three-solve bundles. Only if
both pass, submit the twelve corner bundles (six solves each), at most eight in
flight on physical cores 4–11, reserving cores 0–3 for smaller work. Stop new
submissions at the first false/unset identity; already executing bundles finish.
No retry, changed gate, extra resolution or full-grid successor is automatic.

## Cost forecast

The passing two-point run measured 1334.834 CPU-seconds for 368160 intervals,
approximately 0.00363 CPU-seconds/interval. This contract has **2364640 new
intervals**, predicting about 8570 CPU-seconds / 18–25 min on eight cores,
allowing imbalance and lower parallelism in its initial two-bundle stage.
Expected wall 1500 s, upper 2400 s, hard cap 3600 s; memory cap 8 GiB. The
forecast is conditional on the prior measured host/image and may fail.
Preserve before/after source/input/reference and all 242 engine artifact hashes,
original runner/container exits, 0.2-s cgroup samples and a dated `/home` archive.
No paid provisioning or rebuilt engine is needed.

## Analytic prerequisite result

The frozen source `66ebd1a` completes twelve analytic solves with **32/32 checks
passing**, exits0/0, in15.0354 s /14.451588 sampled CPU-seconds /
46,170,112-byte recorded kernel high-water. At 800 bohr/h=.00625, scalar
continuous K error is1.26839e-7; the discrete oracle attributes1.25185e-7 to
truncation and1.65453e-9 to floating-point drift. Both 80/110-digit oracle
comparisons pass. The [independent verification](../../projects/ukrmol_co/verification-numerov-endpoint-analytic-20261019.json)
reconstructs every error from saved matrices and 110-digit analytic oracles.
All36 archived sources and242 installed engine artifacts verify before/after.
The [CO contract](../../projects/ukrmol_co/numerov-corner-contract.json) pins this
prerequisite,24 raw files,12 full-grid native16/32 references, two native64
references and the previous increment report. The 55-test focused suite passes.
