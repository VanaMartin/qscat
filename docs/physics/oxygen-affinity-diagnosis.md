# Oxygen affinity: separate implementation, basis and correlation errors

## Completed orbital-convention control

The saved-state 5Z replay gives **1.426635554 eV** after semicanonicalization,
versus **1.426606459 eV** originally. The change is **+0.029095 meV**, below the
predeclared 0.1-meV meaningful-effect threshold. Original/native triples logs,
energies, source/image/input pins, rotation/density witnesses and original exits
are independently reconstructed by `projects.ukrmol_co.verify_oxygen_affinity`.
The reference density and CCSD energy are invariant; spin-Fock off-diagonal
blocks fall from up to 0.164 Ha to below 1.94e-13 Ha. This rules out this
orbital-convention effect as the explanation of the affinity deficit, not all
possible coupled-cluster implementation errors. Original endpoints stay intact.

No SCF or CCSD iterations were performed. Four physical Sadaharu cores, 16 GiB:
owner 76.515963 s, worker CPU 138.717374 s, cumulative kernel high-water
12,688,596,992 bytes (11.817177 GiB). All original exits zero; no OOM/timeout.
An immediate terminal observation occurred before the first 600-s periodic check.
See [the control design](oxygen-affinity-orbital-convention.md).

## What is and is not incorrectly set

The saved isolated atoms have consistent electron counts, 2s2/2p4 and 2s2/2p5
occupations, spin, one frozen 1s pair and E(O)−E(O−) subtraction. No R-matrix
radius enters these calculations. There is no demonstrated configuration bug.

The unsupported assumption would be the **accuracy expectation**: a frozen-core,
nonrelativistic finite-basis CCSD(T) calculation is not a complete prediction of
an isotope/lowest-J experimental threshold. The saved 5Z affinity decomposes as
HF −0.538835931 eV, CCSD 1.263421920 eV, and triples +0.163184539 eV. Correlation
supplies the binding; diffuse extent alone need not resolve the residual error.

[de Oliveira et al.](../../reference/literature/deoliveira-1999-pra60-1034.md),
pp. 1036–1037, Table I and the oxygen discussion, report slow valence-basis
convergence (about 16 meV still added beyond AV6Z) and a roughly 12-meV
post-CCSD(T) contribution, mostly quadruples/higher rather than iterative triples.
Their valence correlation limit exceeds our saved 5Z correlation contribution
by 28.467610 meV. This makes the observed deficit plausible without a state bug;
it is **not** our measured correction or a validated error bound. Core and
relativistic/fine-structure contributions have different signs and remain separate.
No published correction is added to our endpoints.

## One bounded next control: radial diffuseness at fixed QZ

**Completed attempt, stopped at the anion; no affinity result.** On3November2026
at14:25:19UTC the frozen source`2b21b43` finished with original worker exits0/1
and science/container/owner/launcher exits2. Neutral O converged, passed internal
stability and CCSD, and gave −74.99518604197657Ha. O− oscillated through100 ROHF
cycles and failed convergence; neither its stability nor CCSD/triples was run.
Its last SCF iterate is not an endpoint and must not be subtracted from neutral O.
Owner15.860873s,CPU44.236479s,kernel high-water5,896,859,648bytes,without OOM or
timeout. Integral-only parent-basis controls passed;no overlap rank deletion or
convergence-gate relaxation occurred. The original failure and successful neutral
are preserved in [the partial-evidence bundle](../../projects/ukrmol_co/qualification-evidence/oxygen-affinity-diffuse-stop/README.md).

The diffuse contribution is therefore **unmeasured**, not small or zero. The
positive occupied orbital energies and occupied/virtual order changes in the
failed log are consistent with a difficult unbound HF reference, but do not prove
that no stationary ROHF solution exists or invalidate the original converged
QZ/5Z states. No retry, optimizer switch or silently accepted unconverged reference.
The following is the exact scope attempted, not continuing launch authorization.

Recalculate only O/O− with the **unchanged original atomic calculation function**,
same occupations, frozen core, convergence/stability gates and noncanonical
ROHF-reference UCCSD(T) convention, adding the second diffuse shell from BSE
d-aug-cc-pVQZ version 0. Keep QZ angular maximum g; do not run 6Z. Original QZ
endpoints are reused, not repeated. The input ledger records the BSE API URL,
response digest and five uncontracted exponents (s/p/d/f/g). The pinned-image
integral-only control proves all parent shells and their AO overlaps unchanged,
80→105 spherical functions, overlap minimum >1e-9; no SCF in that control.

Observable: matched extended−original-QZ affinity, with per-atom HF/CCSD/triples
changes reported separately. A shift exceeding 1 meV is a practically resolved
diffuse contribution; a small shift rules out a large diffuse omission **at QZ**,
not all radial/angular basis error or its effect at 5Z. Neither outcome establishes
absolute accuracy or authorizes a basis ladder. No experimental target tuning.

Caps: two new SCF/CCSD/triples endpoints, sequential, four physical cores0–3,
16-GiB container/12-GB PySCF, 24-GiB available RAM/10-GiB disk floors; 900 s
per atom, 1800 s science, 1950 s container, 2100 s owner service. Conditional
1–5-minute forecast based on smaller-than-5Z AO count, not a speedup promise.
Implementation cap one hour. Durable 600-s observations and immediate completion;
preserve all original exits and first-failure stop, no retries/rescue or successor.
No molecular/scattering/fit/holdout work. A later small-basis FCI correlation
control requires its own release; it is not silently chained to this run.

## Current decision and remaining deliverable

**Stop this diffuse-SCF route. Absolute oxygen affinity remains unresolved.**
The defensible conclusion is a finite-basis/correlation accuracy limitation,
with the orbital-convention contribution ruled out at0.1meV and a separate
diffuse-SCF convergence failure,not a demonstrated occupation bug.

The next bounded deliverable is at most60minutes of saved-source/method selection,
zero new quantum calls: choose a state- and core-consistent correlation control
(for example frozen-core FCI versus CCSD(T) in AVDZ),with analytic spin checks,
a residual/convergence gate and a measured-cost pilot contract. Such a control
measures higher correlation in its finite basis,not a CBS affinity;selecting it
does not qualify a corrected5Z endpoint. Any complete threshold prediction must
also resolve basis limit,core,scalar relativistic,fine-structure and isotope
conventions. Do not add the benchmark's corrections to our data or report that
experimental agreement has been calculated. No active quantum job remains.
