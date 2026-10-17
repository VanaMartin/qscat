# Six-case saved-native exterior reconstruction

**Status: prospective.** The passing
[subrange attribution and public replay](NATIVE_EXTERIOR_ATTRIBUTION.md)
justify this bounded all-anchor/sector reconstruction. The original native
pilot remains false/exits1/1. No scientific result from the new study exists yet.

The [design](../../docs/physics/co-native-full-reconstruction.md) and
[machine-readable contract](native-full-reconstruction-contract.json) define:

- Six saved nine-A2 anchor/sector cases, all99 energies and all298 channels.
- Literal native baseline and one-line16/32-subrange controls,18 invocations.
- Eighteen three-rung independent probes in step, endpoint and tail order.
- Only after all native/probe gates pass, paired-step completion of594 matrices.
- Separate K/S/determinant-phase1e-6 vetoes and1e-6 invariant gates.
- 1278 independent solves if complete;4.5 h forecast,10 h cap,16 GiB,
  eight single-threaded workers on Sadaharu physical cores4–11.

The stretched5-eV probe is only0.02524 eV below its first excited threshold.
The endpoint800-bohr and tail-order controls therefore apply to every probe;
equilibrium agreement does not qualify the stretched case.

Pre-launch validation passes **28 focused tests**, including phase-only vetoes,
every independent refinement knob, coarse-matrix invariants, unique coverage,
stop-on-false/unset scheduling, stage authorization and actual native timeout /
zero-exit-incomplete-output preservation. All24 saved-input hashes verify locally;
all six decoded models have99 six-open-channel K records and the expected
298-channel/2048-pole boundary. Lint, formatting and diff checks pass.

Execution uses a frozen source archive with the unchanged validated Numerov
method. `native_full_reconstruction.py` writes original native exits, per-energy
JSONL journals and incremental reports. `native_reconstruction_checks.py`
evaluates fixed configurations and coverage. False/unset stops new submissions;
already executing fixed bundles finish and are retained. No retry or uncontracted
successor is authorized. Work uses saved exports only, with no new inner solve.

```bash
python -m projects.ukrmol_co.native_full_reconstruction \
  --directory /native --prerequisites /prerequisites \
  --engine /opt/ukrmolp/bin.double/rsolve \
  --contract projects/ukrmol_co/native-full-reconstruction-contract.json \
  --output /outputs/result
```

The mounted input/prerequisite/engine paths must match contract hashes. This
finite-grid comparison cannot establish absolute omitted-spectrum errors,
electronic/channel convergence, global resonance identity or complex poles.
Five-A2 interfaces, width qualification and production remain gated.
