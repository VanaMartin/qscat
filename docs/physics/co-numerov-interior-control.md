# Interior midpoint and native-K peak control

## Question and selection

The [corner/endpoint study](../../projects/ukrmol_co/NUMEROV_CORNER_CONTROL.md)
qualifies fourteen previously visited identities, not the intervening energies.
Read all594 saved native32 matrices before any new candidate computation. For
each of the six sectors, select the largest absolute K entry on interior indices
1–97 (lowest index breaks a tie), plus the fixed index48 midpoint. Retain an
already-qualified peak as reuse rather than choosing a different peak.

The read-only inventory selects compressed index87 (requested4.45 eV, maximum
entry207.915, largest eigenvalue≈313), equilibrium index53 (2.75 eV, entry26.710,
eigenvalue≈−51.446) and stretched index18 (1 eV, already qualified). Thus the
finite new scope is **four peak identities plus six midpoint identities**.
Four peaks are outside the original132-identity independent coverage; the six
midpoints were original probes but have not been increment-qualified.

Selection depends only on saved native K, with exact energy bits and input hashes,
not favourable candidate errors. It is a stress sample, not a global error bound
or an argument that the largest K necessarily maximizes every observable error.

## Prospective spacing rationale

For the asymptotically free open channels, the leading accumulated Numerov phase
defect is `L E^(5/2) h^4 /480` in native Rydberg coefficients and bohr lengths.
For a common infinitesimal phase shift, `dK/d(delta)=I+K²`. Use
`maxabs(I+K²) E^(5/2)` only as a spacing-cost diagnostic, not a CO error bound:
coupling, finite-radius terms and native-reference uncertainty remain.

The compressed peak's diagnostic is roughly100 times that of the already tested
stretched peak. Halving h reduces the leading defect by16. A prospectively fixed
compressed sequence `.00625/.003125/.0015625` is therefore justified before a new
run rather than discovering failure through an automatic spacing ladder. At the
equilibrium peak use `.025/.0125/.00625`; at midpoints use the previously passing
corner sequence `.05/.025/.0125`. All use the existing six-coordinate control:
three spacings at400/order8, finest200/800/order8 and finest400/order6.

## Analytic prerequisite before CO use

Test scalar and rotated two-open/one-closed free boundaries at the two new peak
energies, with target tangent+320 (compressed) and−55 (equilibrium). These bracket
the saved dominant eigenvalues without copying a CO surface map or potential.
Each uses its six declared coordinates: **24 analytic solves**. Compare continuous
boundary matching and scalar discrete Numerov transfer at80 and110 digits.
Require oracle drift<1e-20, candidate/scalar-discrete K drift<=1e-7 on every row,
finest-endpoint continuous K/S/phase error<=1e-6, and invariants<=1e-6. Coarser
continuous errors are diagnostics, not claimed passes. Failure/unset stops before
CO submission. Four physical cores0–3/4 GiB, forecast120 s /cap600 s.

The propagator, tail, unsymmetrized matching, physical inputs and native decoder
remain unchanged. All K/S/phase/invariant budgets remain separately1e-6; retain
the original rejections and exits independently.

## Conditional finite CO contract and resource forecast

Only after analytic qualification, pin its hash and all inputs/references with a
separate frozen contract before launch. At most **60 new solves**, all six
coordinates at the ten selected identities. Check all native16→32 reference
controls before propagation. Peak stage precedes midpoint stage. Use eight row
workers on physical cores4–11 (reserving0–3), with8 GiB and single-threaded BLAS;
schedule the longest fixed rows first so compressed endpoint800 does not become
an avoidable serial tail. Every complete identity retains all21 prior case gates,
including every finest200/400/800 endpoint versus native32.

Stop new submissions at the first false/unset row invariant or completed-identity
gate; fixed in-flight rows finish and are preserved. No new native invocation,
inner solve, retry, extra rung, changed tolerance or automatic full-grid successor.

There are **4190160** new radial intervals. The prior measured8441.309248 CPU s /
2364640 intervals predicts about14960 CPU s. Heavy compressed rows and sequential
scientific stages limit parallelism: forecast2400 s /upper3600 s /cap5400 s.
Preserve0.2-s cgroup samples, original exits, dated `/home` artifacts and before/
after all source/input/reference and242 installed engine artifact hashes.

If complete, historical solved coverage can increase from132 to136 identities,
leaving458 uncomputed, while increment qualification can increase14→24. Those
counts must be derived from actual completed rows and decisions, not assumed at
submission. The original archived132/462 counts and false verdict stay intact.

## Recorded prerequisite rejection — 19 October, 09:41 UTC

The analytic source freeze **`b7f2278`** completes all24 solves with **64/68 checks
passing**, aggregate **false / original runner-container exits1/1**. The four
vetoes are confined to the compressed(+320 tangent/high-momentum) controls:

- Scalar h=.0015625/800: candidate/discrete K drift **1.20892527e-7**, above1e-7,
  despite continuous error5.76310e-8 passing its separate1e-6 budget.
- Rotated open/open/closed h=.0015625/200 and400: continuous K errors
  **1.12543592e-6 /1.11113104e-6**, above1e-6. The400/order6 row duplicates the
  same numerical veto; it remains a separate declared check.

The rotated400-bohr error worsens from4.92619e-7 at h=.003125 to1.11113e-6 at
h=.0015625. All S/phase/invariant gates and80/110-digit oracle checks pass.
Those passes cannot waive either the scalar precision or rotated K vetoes.
The scalar precision failure identifies floating-point drift there; the coupled
continuous failure is not by itself a separately measured roundoff bound.

Measured **72.952831 s /72.491243 sampled CPU-seconds /47,329,280-byte kernel
high-water**,363 samples at0.2 s, four physical cores0–3/4 GiB. All46 frozen
sources and242 installed engine artifacts verify before/after. The10-file
collection preserves original source, outputs, exits and resources.

The stopping rule is applied: **zero CO successor submissions, zero new CO
independent solves and zero native/inner solves**. No CO contract was frozen;
the conditional60-solve proposal is withheld. Increment qualification remains
14 identities; historical coverage132/594 and462 uncomputed remain unchanged.
No tolerance relaxation, new rung or retry follows this contract. Further useful
work requires a separately costed fixed-grid precision/algebra attribution or
independent propagation control in this analytic regime before returning to CO.
