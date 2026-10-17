# Native exterior error attribution passes — 17 October 2026

**All ten frozen checks pass; original command/container exits0/0.** Native
subrange refinement accounts for **99.72% of the measured5-eV S mismatch**
against the fixed finer independent reference. The original pilot remains
false/exits1/1. See [verification](verification-native-attribution-20261017.json)
and [public evidence](qualification-evidence/native-attribution/README.md).

## Completed result — 14:01:09 UTC

| Comparison at requested5eV | K-entry difference | S-entry difference | Phase / rad moduloπ |
|---|---:|---:|---:|
| Original8-subrange native vs finer independent | 7.23362e-7 | 1.44664e-6 | 2.58085e-6 |
| Native16→32 subranges | 1.43254e-11 | 2.85943e-11 | 3.84278e-11 |
| Independent step0.025→0.0125 at400bohr | 3.26498e-8 | 4.71630e-8 | 1.40888e-7 |
| Independent endpoint100→400 at step0.0125 | 1.73821e-9 | 2.54094e-9 | 7.51674e-9 |
| Native32 vs finer independent | **2.68388e-9** | **4.08632e-9** | **1.15175e-8** |

The step S change improves by16.07-fold over the original last step change.
Native baseline reproduces all99 original K matrices **exactly**; baseline
S difference is zero and determinant-phase roundoff is2.475e-17. The native
16/32 runs each retain99 complete energy records and differ only by their
single declared input line. Across the entire99-point grid,16→32 differences
are5.58703e-9 in K,3.13854e-10 in S and3.47600e-10rad in phase. This native
refinement evidence does not establish independent accuracy at every energy.

The99.72% reduction is measured against the **same h=0.0125,400-bohr** reference
on both sides. The earlier h=0.025 pilot's1.49310e-6 S discrepancy and false
verdict are not replaced. This is a controlled attribution to native radial
resolution for the measured case,not proof that native32 is exact or that no
other uncertainty remains. Native asymptotic defaults were unchanged.

The independent100/400-bohr matching raw condition estimates are5.64612e8 /
1.21054e4; equilibrated1.89222e8 /2.93565e3. Independent unitarity/reciprocity
errors stay≤1.121e-12/1.327e-12. Full matrices are retained unsymmetrized on
the Python side. A separate right-Cayley/eigenphase recomputation verifies the
finest comparison within1e-14. The two-line journal exactly matches both cases.

Measured Sadaharu cost is **1339.551s /0.56024GiB sampled kernel high-water**
(601,554,944bytes),6689 samples at0.2s,sampled CPU1334.964s. Native8/16/32
cost198.495/347.997/659.113s; independent100/400 cost23.628/109.689s. The
total is below the1667s estimate/3600s cap; sampling is not a post-exit peak.

Ten public payloads retain all297 native energy records,two independent cases,
decks/logs,eleven frozen source/design inputs,original native inputs and the
prior rejection/public-replay prerequisites.
[Fresh public replay](verification-native-attribution-public-20261017.json)
verifies all ten payloads,eleven sources,nine native output files,four original
native inputs and297 energy records. Both independent cases are recomputed;
all ten gates and the attribution pass,exit0. Maximum host/public K-entry /
S-component drift is9.48917e-10/1.52518e-9. Four public attribution tests pass.
This reexecutes independent propagation and reconstructs saved-native gates;
it does not rerun the native executable on the laptop.

The [14:04:19 UTC host check](progress-native-attribution-result-20261017.json)
verifies17 new hashes/297 prior hashes/242 installed artifacts. No CO job runs
on Sadaharu;unrelated `comfyui` remains and free persistent disk is40.76GiB.

**Next decision:** prospectively design the full saved-native reconstruction
using adequately resolved native outer references. The passing scope is one
sector at5eV,with native-only refinement over99 energies. Freeze controls for
the other anchors/sectors,especially the stretched5-eV point near its first
excited threshold,plus unchanged K/S/phase budgets and costs before execution.
Original8-subrange outputs remain evidence rather than assumed1e-6 references.
No new inner solve is justified. CO poles and production remain unreleased.

## Decision

The [original pilot](NATIVE_EXTERIOR_PILOT.md) rejects S/phase reproduction at5eV.
[Fresh public replay](verification-native-exterior-pilot-public-20261017.json)
reproduces every case/check,false and exit1 from nine checksum-verified payloads.
Maximum host/public K-entry/S-component drift is1.61643e-9/3.23265e-9. The false
gate remains; the published method controls remain separately passing.

The [source-derived attribution design](../../docs/physics/co-native-exterior-attribution.md),
[contract](native-exterior-attribution-contract.json),runner and tests are frozen
at **`666c056` before execution**. **21 focused tests pass**, including four
new attribution tests for literal deck fidelity,full-record comparisons and the
separate phase gate. No engine rebuild or new inner-region calculation occurs.

## Fixed experiment and gate

- Same pinned address-audit engine and equilibrium/nine/B1 saved298-channel /
  2048-pole inputs.
- Exactly three native outer-only runs: literal original deck (automatic8
  subranges),then one-line16/32-subrange controls. All retain the original99
  energies,100-bohr endpoint,NLEG10,Gailitis defaults and multipoles.
- Two independent Numerov cases at5eV:step0.0125bohr,endpoints100/400bohr,
  tail order8. The validated solver is unchanged.
- Literal native baseline must reproduce all99 saved K/S/phase records within
  1e-10 before the two refined native runs execute.
- Final native refinement,independent step/endpoint and native32-to-independent
  comparison retain separate1e-6 K/S/phase gates as defined in the contract.
- Attribute most of the measured discrepancy to native subrange resolution only
  if all gates pass and the S mismatch falls by at least80%. This is a
  conditional attribution,not an absolute error bound or an exclusive cause.

Native triangular K output cannot measure unsymmetrized native reciprocity.
Independent K is left unsymmetrized. Full594-matrix reconstruction,CO poles,
model/channel completeness and production remain gated even if this control
passes. Execution/interface failures stop before successors; false/unset permits
no additional subrange/order/radius ladder or relaxed budget.

## Cost

The historical B1 native outer stage costs218.9348s for8 subranges/99 energies.
Linear subrange scaling predicts1532.54s for all three native runs. The two
independent cases add37,120 intervals,about133.94s by the completed pilot's
measured rate. **Total estimate1667s; upper2700s; hard cap3600s; memory4GiB**,
physical cores0–3 reserved and one BLAS thread. Outputs,inputs,original exits
and sampled resources are retained in the persistent campaign directory
`prepared/native-exterior-attribution-source-20261017/`.

## Launch observation — 13:39:45 UTC

The [host capture](progress-native-attribution-launch-20261017.json) verifies
**14 new immutable remote hashes /283 prior hashes /242 installed artifacts**,
the read-only native data,prior reports and engine mounts,physical cores0–3 and
4-GiB cap. The container starts **13:38:50 UTC**. The literal8-subrange baseline
is running; no baseline or attribution verdict is yet available. No independent
case has started. The unrelated `comfyui` workload remains active; persistent
free disk is40.78GiB. The completion hook will return original exits and results.
