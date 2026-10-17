# Native exterior error attribution — running 17 October 2026

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
