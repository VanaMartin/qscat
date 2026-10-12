# CAS12 audited-engine baseline control

The completed address audit exported4096B1 roots after passing all structural
and address guards. It did not reproduce either prior fault. Assertions,
compiler choices and debugger/runtime instrumentation can change layout or
timing even when computational statements remain unchanged. The result supports
a reproduction control; it does not establish a fault cause or engine correction.

The [separate declaration](cas12-audited-engine-baseline-contract.json) tests
the exact242installed artifacts of that completed build against the preserved
successful2048-root scientific baseline. It retains `-O1 -g -fbacktrace
-fcheck=bounds`, all address assertions and Class1/Class5–7 markers. The engine
is mounted read-only, without a rebuild or patch. Numerical-library bytes are
pinned independently of the image tag. All12serial/MPI engine checks and three
debugger tests remain evidence for this unchanged build.

## Runtime and model

Use the original scientific `sparse_pilot` driver with direct eight-rank MPI
and the original read-only SLEPc telemetry interposer. Unlike the completed
GDB diagnostic, this runtime records solver coverage, selection and absolute
action residuals. This is an explicitly declared runtime change; the purpose
is to qualify the engine/runtime intended for scientific comparison.

Both B1 and B2 use2048roots, the same frozen target/integral/CONGEN inputs,
CAS12-specific DENPROP map and native/outer decks as the successful baseline.
The electronic model remains equilibrium2.1323bohr/cc-pVDZ/CAS12[6,3,3,0],
two frozen A1 core orbitals,ten active electrons,40-component averaging and
40explicit target states. Radius18bohr,maxl4,deletion1e-6,propagation100,
crite1e-12,maxiter1000 and the uncorrected pole sum are retained.

The control uses eight physical cores4–11,one library thread per rank,80-GiB
cgroup,8-GiB shared backing and16-GiB native `memp`. Native-stage/whole-native/
owner limits are24/50/72hours, with a one-hour independent reader. Physical
cores0–3 remain available for the small readers and other small jobs. An
observation-only hourly owner is bounded to72hours and starts only after the
calculation owner's execution record is readable.

## Ordered gates and verdicts

1. Reverify the completed engine/source/check evidence, installed artifacts,
   successful scientific baseline and the fresh source/input freeze. Check
   physical-core topology, other workloads,≥96GiB available RAM and≥50GiB
   free disk. Verify numerical-library bytes and telemetry symbol loading.
2. Independently read the baseline against itself under the new reproduction
   budgets, proving that the retained raw data and same-window fits are usable.
3. Run one fresh two-sector2048-root replay. Preserve native exits, raw logs,
   timings and memory peaks. Native errors/STOP diagnostics, pointer guards,
   bounds errors, timeout or export/residual/Pi violations end the owner.
4. Independently re-read both runs and verify all five native input hashes per
   sector, exact decks, root coverage, energies, channel/threshold/multipole
   digests, full99-point phase grids and modulo-π Pi agreement. Compare the
   same1.8–3.3requested-eV fit window for background orders1–4.

Same-root reproduction must be much tighter than a scientific resolution
comparison: **0.1meV position and full width**, **1e-6rad full-grid phase**.
These budgets exceed printed-phase rounding while staying well below the
scientific0.05eV/5%/0.05rad budgets. Shared energies and absolute residuals
retain1e-7Hartree; CI/boundary energies retain1e-10Hartree. Native selection,
dimension and solver type are checked. Empty/MAXFIT/failed/bound fits leave
the reproduction verdict unset; false or unset verdicts stop the owner.

`engine_baseline.py` labels same-root agreement as
`matched_engine_baseline_verdict`. Root-pair stability, omitted-spectrum,
engine-correction and electronic-model verdicts remain unset. This comparison
uses the same target and discretization, so it qualifies engine reproduction,
not electronic accuracy or model balance.

The owner schedules only this replay and its independent reader. `STOP` in the
fresh source directory stops the active step. After a pass, separately declare
a4096-root scientific run on the same qualified build/runtime. Any continued
fault investigation, source correction or omission test remains a separate
experiment; no success automatically releases the sweep.
