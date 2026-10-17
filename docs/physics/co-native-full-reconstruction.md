# Six-case CO native exterior reconstruction

## Decision and fixed model

The [completed attribution](../../projects/ukrmol_co/NATIVE_EXTERIOR_ATTRIBUTION.md)
and fresh public replay pass: native subrange refinement removes 99.72% of the
equilibrium-B1 5-eV S mismatch against the same finer independent reference.
This motivates testing all six nine-A2 saved anchor/sector combinations with
resolved native outer references. The original pilot false/exits 1/1 remain.

The [contract](../../projects/ukrmol_co/native-full-reconstruction-contract.json)
pins 24 inputs and the attribution/replay prerequisites. Each case retains the
original forty-component orbital average, 48 retained components, 298 channels,
2048 uncorrected poles, neutral 18-bohr zero-Bloch boundary and both coupled
multipoles. All closed columns remain. B1/B2 K output units come from their
decks (921/922). The validated [independent equations and solver](co-native-exterior-reconstruction.md)
are unchanged. Native NLEG10, 100-bohr matching and Gailitis defaults stay fixed.

## Prospective stages

1. For compressed/equilibrium/stretched × B1/B2, run the literal native deck
   and exactly two one-line NRANGE16/32 variants. The literal baseline must
   reproduce all 99 saved K/S/phase records within 1e-10 before its refinements.
   Validate exits, printed subranges, complete records/EOF and input fidelity.
   Require native16→32 K/S/phase differences ≤1e-6 at every energy.
2. After all six native gates pass, probe indices 0/48/98 (requested
   0.1/2.5/5 eV) in each case. Seven independent one-factor configurations are:
   steps 0.05/0.025/0.0125 bohr at endpoint400/order8; endpoints200/400/800 at
   step0.0125/order8; orders4/6/8 at step0.0125/endpoint400. The common comparison
   point is 0.0125/400/8. Retain both successive changes in every family.
   Require final step, endpoint400→800 and order6→8 changes to satisfy **each**
   K/S/phase 1e-6 budget. Require S spatial improvement ≥8, subject to 1e-8
   floor, and common finest versus native32 K/S/phase ≤1e-6.
3. Only after **all 18 probes pass**, reconstruct the remaining 96 energies
   per case at steps0.025/0.0125, endpoint400, order8. Require separate final-step
   and finest-versus-native32 K/S/phase ≤1e-6 for every energy. Reuse probes;
   success requires **594 distinct energy identities**, **1278 independent
   solves** and **18 native invocations**.

Unitarity and unsymmetrized independent reciprocity ≤1e-6 are separate gates
on **every** independent matrix, including coarse controls. Native triangular
storage cannot check unsymmetrized native reciprocity. Keep matching conditions,
residuals, full unsymmetrized K, S components and determinant phases.

The stretched first excited threshold is 5.025239863945225 requested-eV, only
0.02523986 eV above 5 eV. Its Rydberg gap is approximately0.00185513 and slowest
closed decay length23.22 bohr. This motivates the800-bohr and tail-order controls.
Compressed/equilibrium first thresholds are7.495618966/6.385015673 eV. No saved
energy crosses a threshold. Small successive differences are conditional
resolution evidence, not absolute error bounds or uniform tail bounds.

## Cost, resources and stopping

Measured attribution native stages sum1205.605 s per case: six cases forecast
7233.63 CPU-seconds. Probes total3,450,960 radial intervals; the paired-step
remainder totals26,403,840. At133.3177 s/37,120 intervals the independent forecast
is107,225 CPU-seconds. Total114,459 CPU-seconds gives3.97 h ideal eight-worker
throughput; allow **4.5 h expected /8 h upper /10 h hard cap** for gates,
six-way native occupancy and contention. These are conditional forecasts,
not an eight-worker benchmark. Native per-stage timeout1800 s.

Use Sadaharu physical cores4–11, eight single-threaded processes, 16-GiB
aggregate container memory; reserve cores0–3 for small work. Record each solve
to its own persistent per-energy journal and flush parent reports as bundles
return. Stop new submissions on the first false/unset gate. Already executing
fixed bundles finish and retain evidence, but cannot authorize a later stage.
Execution/interface failures leave the stage unset; numerical rejection is false.
Record the subprocess's actual exit separately from timeout classification.
No retries, larger grid ladders, relaxed gates or uncontracted successors occur.

Verify the pinned 242-artifact installed engine, frozen source, read-only input
mounts, host topology/workloads and memory/disk before launch. Keep all partial
and final evidence under a new persistent `/home` directory. Reverify inputs and
engine afterward. This is saved-data outer work; no new inner solve is justified.

A passing verdict qualifies these finite-grid real-axis comparisons with the
declared probe/paired-step coverage. Full channel/electronic convergence,
omitted-spectrum error, global continuity, widths and complex poles remain
unqualified. Five-A2 unset readings, original twelve-case aggregate false/exits1/1
and the original native pilot rejection remain. Production remains gated.
