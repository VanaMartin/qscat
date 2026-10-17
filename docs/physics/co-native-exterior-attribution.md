# Attributing the saved-native CO exterior discrepancy

## Decision-changing observation

The [original native pilot](../../projects/ukrmol_co/NATIVE_EXTERIOR_PILOT.md)
rejects: at requested5eV, finest native S-entry and phase errors are
1.4931025067e-6 and2.7217428483e-6rad moduloπ, above the frozen1e-6 budgets.
K-entry, unitarity, reciprocity and independent step/endpoint/order changes pass.
The original false verdict and exits1/1 are retained. This separately designed
study tests **native radial-propagator resolution** against a modest additional
independent spatial refinement; neither side is presumed exact.

## Source and model controls

The pinned `outer/source/compak/rsolve_compak.f` has SHA256
`93c462b3a76a35cf068d195c3469ae57f1ef12aa19af118a43dc33cb57b59684`.
`ASYM1_compak`, lines1321–1324,1328–1333,1344–1353,1369–1372,1395–1405
and1439–1442 establish:

- `NLEG` defaults to its maximum10;
- `NRANGE=0` chooses8 subranges for this saved energy range and radius;
- explicit positive `NRANGE` retains the supplied outer endpoint when it is
  greater than the18-bohr inner radius;
- `RPROP1` builds the global propagator and `CURLYR` evaluates it at each energy;
- `IASY=20` and `EPS=5e-5` are fixed source defaults, absent from the `BPROP`
  namelist. This study does not claim that setting another input named EPS
  changes those defaults.

The saved native log explicitly reports8 subranges and Gailitis matching at
100bohr. Its zero-K count does not supersede the complete99-record binary K
dataset. Native controls use the **same original99-point0.1–5 requested-eV
grid**, the same298-channel/2048-pole exports, the same100-bohr endpoint and
the checksum-pinned address-audit engine. Matching conventions are unchanged.

The baseline uses the literal original deck bytes. Two new decks insert exactly
one `nrange=16` or `nrange=32` line after the `&bprop` header. All other bytes
remain identical. No CI, target optimization, boundary export, eigenphase fit,
engine rebuild or complex-pole search is involved. Each native invocation has
a fresh directory, original exit, stdout, full K and input hashes.

## Independent control and success criteria

Retain the validated Python Numerov/tail method unchanged. Add only two5-eV
cases: step0.0125bohr, endpoints100/400bohr, tail order8. Compare with the
already saved0.1/0.05/0.025 families; do not choose a different step after
seeing these outputs. Both new cases keep all298 channels and closed columns.

The prospective [contract](../../projects/ukrmol_co/native-exterior-attribution-contract.json)
defines all inputs, the source and engine pins, budgets and stopping rules.
Before new numerical work, the original pilot's fresh public replay must
reproduce all21 case records, all ten check verdicts, false and exit1.

1. Reproduce all99 saved baseline K/S matrices and phase sums at1e-10 absolute
   budgets with the same engine/deck. Check complete records and energy sequence,
   not only process exit. A baseline failure stops before refined native runs.
2. Require finite complete99-record outputs from native16/32, correct printed
   subrange counts, unchanged inputs and original exits0. Preserve warnings and
   stop on an execution or interface failure.
3. At5eV, require native16→32 K/S/phase changes≤1e-6. Report all99-point
   differences separately; this targeted study does not qualify their agreement
   with independent propagation.
4. Require the new independent step change≤1e-6 in S and phase and its S change
   to improve by≥8 relative to the prior step change, subject to1e-8 floor.
   Require the new100→400 endpoint change≤1e-6 in S and phase.
5. Require native32 versus the finest independent400-bohr result to satisfy
   **each** original K/S/phase1e-6 budget. Require real-axis unitarity and
   reciprocity≤1e-6 for both new independent matrices and the native32 matrix.

Native K records are triangular; reconstructing their symmetric form cannot
measure unsymmetrized native reciprocity. Independent Python K is left
unsymmetrized, as in the original pilot.

Report the default→refined native correction, the change on the independent
side, the residual vector and fractional reduction of the original S mismatch.
If gates pass and the baseline mismatch is reduced by≥80%, label the result
“native subrange refinement accounts for most of the measured5-eV discrepancy”.
Otherwise retain an unresolved attribution. This classification is conditional
on this fixed model, saved spectrum and asymptotic implementation; a successive
difference is not an absolute error bound or proof of an exclusive cause.

## Cost and stopping scope

The historical native equilibrium-nine/B1 `rsolve` stage measured
218.93477273s for8 subranges and99 energies. With fixedNLEG/channel count,
rough linear scaling in the subrange count gives1532.54s for the baseline,
16 and32 combined. The two independent cases have37,120 radial intervals;
the completed pilot's734.60s/203,580 intervals imply133.94s. The combined
estimate is **1667s**, upper estimate2700s, hard cap3600s,4GiB, one BLAS
thread and physical cores0–3 reserved on Sadaharu. This is a cost forecast,
not a new timing benchmark; allocation, serial diagonalization and conditioning
may alter the realized cost. All artifacts and sampled resources stay under
the persistent campaign root.

This authorizes exactly three native outer-only invocations and two independent
cases. Any false/unset gate stops without retries, additional subranges, changed
tail order, larger radii, relaxed budgets or automatic successors. Even a pass
qualifies only this5-eV attribution control; the original pilot remains false,
full six-case/594-matrix reconstruction remains unreleased, and native poles,
electronic/channel convergence and production remain unqualified.
