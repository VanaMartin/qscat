# Six-case saved-native exterior reconstruction

**Status: stopped on two K-matrix rejections, 17 October 2026, 17:48:36 UTC.**
The original command/container exits are **1/1** and the reconstruction verdict
is **false**. All six native controls and all 18 independent probes pass. The
paired-step remainder stops after 112 passes and two rejections: stretched B1/B2
at requested **1.0 eV**. Coverage is **132/594 distinct energies**, with **354
independent solves**. The remaining 462 energies were not computed.

The source/design/contract/tests freeze is **`fd7c871`**, committed and pushed
before execution. The [passing predecessor](NATIVE_EXTERIOR_ATTRIBUTION.md) and
original pilot false/exits 1/1 remain. See the
[saved-decision verification](verification-native-full-20261017.json) and
[published evidence](qualification-evidence/native-full-reconstruction/README.md).

## Completed controls and strict rejection

All 18 native invocations exit 0. The six literal baselines reproduce **594 saved
K matrices exactly**. Native 16→32 maxima across all six complete grids are
**1.13753e-7 K / 3.94956e-10 S / 3.47599e-10 rad phase modulo π**, below the
separate 1e-6 budgets. All 18 independent step/endpoint/order probes pass,
including the stretched 5-eV points only 0.02524 eV below the first excited
threshold. That probe coverage did not predict the subsequent K-only rejection.

| Stretched 1-eV comparison | B1 | B2 |
|---|---:|---:|
| Final step change in K | 1.05931e-5 | 1.22856e-5 |
| Finest independent versus native K | 8.16850e-6 | 6.10584e-6 |
| Final step change in S | 6.29827e-10 | 6.34026e-10 |
| Finest independent versus native S | 3.20011e-9 | 2.92038e-9 |
| Finest independent versus native phase / rad modulo π | 3.48745e-9 | 3.50252e-9 |

Both `step_k` and `native_k` reject in each sector. All other checks pass,
including invariant errors below 7.94e-12 at the rejected points. Separate K,
S and phase gates were frozen; an S/phase pass cannot waive the K rejection.
The runner stopped new submissions and retained already executing bundles.

## Saved-matrix diagnosis

The independent [read-only replay](native_reconstruction_replay.py), frozen at
**`64226b9`**, uses a right-solved Cayley transform and eigenvalue-arctangent
phase sum, separate from the producer's left solve and determinant ratio. It
reconstructs all **1782 native records**, all 354 journalled solves, every
per-energy decision and the false aggregate from the saved matrices. Fresh
public bytes reproduce those decisions. Three analytic
tests check the large-tangent-channel limit, noncommuting matrix identity and
phase wrapping.

[Fresh public verification](verification-native-full-public-20261017.json)
checks **13 payloads / 22 producer, replay and test-support files / 24 inputs**.
Reexecuting the two rejected bundles (four solves, unchanged numerical settings)
on Sadaharu reproduces K and stored S components **exactly**, with identical
false decisions. Verification exits 0/0; original scientific exits 1/1 remain.
The replay took **190.720 s**, sampled CPU 374.781 s, recorded kernel memory
high-water 141,156,352 bytes, two single-threaded physical cores and a 4-GiB cap.
It does not reexecute native Fortran or all 354 independent solves.

All **31 focused tests** pass with two checksum-pinned test-only dependencies
supplied in a separate archive. The original producer source archive remains
unchanged; initial public test-package failures are retained in `public-tests.json`.

At 1 eV, maximum native K entries are about **141.309** and the largest-magnitude
K eigenvalue is **−185.358**. Its eigenphase lies **0.00539490 rad** from a tangent
singularity; adjacent saved energies have maximum K entries about 2.20/2.40.
Relative Frobenius K discrepancies are only **5.77e-8 / 4.31e-8**. The projection
onto that large native eigenchannel accounts for **0.99999865 / 0.99999713** of
the discrepancy's Frobenius norm.

The exact identity

\[
S(A)-S(B)=2i(I-iA)^{-1}(A-B)(I-iB)^{-1}
\]

reconstructs the S difference to **5.69e-15 / 2.63e-15** maximum-entry residual.
The scalar local Cayley derivative in the large channel is about **5.82e-5**.
This explains the disparity between absolute K and S changes in these saved
matrices. It does not locate a complex pole, establish an absolute propagation
error, identify the source of the remaining numerical discrepancy, or convert
the strict result into a pass. Any successor needs its own prospective scope,
observable/error model, cost and stopping rule.

## Measured resources and preservation

The original run used physical cores 4–11, eight single-threaded workers and a
16-GiB container limit. Wall time was **6877.185 s (1.91033 h)** and sampled CPU
time **45030.741 s**. The 0.2-s resource record has 34325 samples. Recorded kernel
memory high-water is **3,378,905,088 bytes (3.14685 GiB)**; this is the partial
stopped study's measured cost, not a full-grid forecast.

Post-run checks verify 242 installed engine artifacts, 24 inputs and 18 frozen
source files. A separate collection verifies and preserves 335 output/metadata
files, including all journals, logs, original exits and resource samples.
The [18:14 UTC remote observation](progress-native-full-result-20261017.json)
reverifies 348 collected producer/replay files, 447 prior immutable remote hashes
and all 242 engine artifacts. No scientific container remains running.

## Frozen launch provenance and scope

The [launch observation](progress-native-full-launch-20261017.json) verifies
22 new /314 prior remote hashes /242 installed artifacts, read-only mounts,
physical cpuset4–11 and16-GiB limit. All six literal8-subrange native baselines
were executing; no stage or independent solve had completed at that cutoff.
Start15:53:59 UTC. Persistent directory:
`prepared/native-full-reconstruction-source-20261017/`; durable owner11237.
Source archive SHA256:
`42b7f5b666b8312a3cf335f387c5208742026996bf61d61bef0b7e64377f0ac2`.
Original exits, per-solve journals and0.2-s sampled resources are retained there.

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
