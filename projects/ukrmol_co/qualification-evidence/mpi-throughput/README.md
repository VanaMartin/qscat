# Matched CAS(10,10) MPI scaling and throughput

All seven tight equilibrium CAS(10,10)/cc-pVDZ scattering replicas pass the
matched numerical gates. Each uses the same final checkpoint from the
[published matched baseline](../README.md), forty equal-weight orbital-ensemble
components, forty computed target roots and forty retained target states,
l=4, radius 18 bohr, deletion 1e-6 and 99 energies. Retained target states
are distinct from their partial-wave scattering channels. External energy
labels span 0.1–5.0 eV under the pinned 18.44-ppm conversion convention;
phases are radians and cross sections bohr². All target/import, grid,
cross-section-sum and Pi checks pass;
both contracted scattering dimensions are **8350**.

## Single-job MPI scaling

The one-/two-/four-rank jobs run sequentially on the same four physical cores,
CPUs **12–15**, with 16-GiB container limits and one thread per numerical library.

| MPI ranks | Runner wall (seconds) | Aggregate CPU (seconds) | Kernel memory peak (GiB) | Full-job speedup | Scattering SCATCI (seconds) | RSOLVE (seconds) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 3105.91 | 3105.86 | 2.849 | 1.000× | 2109.96 | 260.87 |
| 2 | 1830.62 | 3435.09 | 2.743 | 1.697× | 1065.44 | 151.72 |
| 4 | 1239.61 | 4268.46 | 2.822 | 2.506× | 560.58 | 109.48 |

Whole-job parallel efficiency is **84.83% / 62.64%** at two/four ranks.
The scattering eigensolve alone speeds up **1.980× / 3.764×**; serial stages
and less scalable work reduce the full-pipeline gain. The batch wall is
**6177.72 seconds** including container startup and final analysis. Native
SCATCI logs independently report the requested process counts.

## Four concurrent one-rank jobs

Four replicas pinned individually to CPUs **12 / 13 / 14 / 15** use 8-GiB
limits each and finish together in **3444.36 seconds / 57.41 minutes**.
Runner walls are 3442.00/3412.68/3380.81/3443.62 seconds; their per-job kernel
peaks are 2.835/2.837/2.835/2.835 GiB.

Against four sequential repetitions of the measured four-rank runner time,
`4 × 1239.61 / 3444.36 = 1.43958`: concurrency delivers **43.96% higher
throughput**, or **30.54% less elapsed time for four jobs**. This denominator
is the measured concurrent batch; the sequential four-job reference is an
extrapolation from the measured single job. Aggregate CPU is **13675.56 seconds**,
versus `4 × 4268.46 = 17073.86 seconds`, **19.90% less**.

Summed per-job kernel peaks are **12178423808 bytes / 11.342 GiB**, and summed
run-disk peaks are **5368893440 bytes / 5.000 GiB**. These are conservative
aggregate envelopes, not synchronized instantaneous peaks: profiles record
per-run elapsed times without a shared absolute sample clock.

These are single-shot measurements of one compact electronic model. Other
calculations on disjoint core groups share the host; exact workload samples
are retained before both batches. The result supports using four independent
workers for multiple qualified jobs of this size and MPI when single-job
latency matters. Larger active spaces, bases or channel counts need their own
memory/scaling measurements.

## Numerical equivalence and evidence

The replicas start from the baseline's final checkpoint, while the baseline
itself started from an earlier continuum-run checkpoint. That difference is
explicitly checked rather than assumed irrelevant:

- All forty QC roots match the baseline within **6.119e-9 Hartree** and ground
  dipoles within **2.986e-8 a.u.**.
- Minimum active-subspace overlap is **0.9999999999998141**; core subspaces
  also satisfy the same numerical restart gate.
- All forty imported roots match each replica's QC within **5.044e-10 Hartree**,
  and DENPROP dipoles within **4.295e-11 a.u.**.
- B1/B2 phases differ from the baseline by at most **2e-7 rad modulo pi**.
- Native fitted position/full width change by at most **0.137 / 0.041 micro-eV**.
- Total-cross-section changes are at most **1.3e-5 bohr²**; complete grids,
  nonnegativity, final-state sums and Pi agreement pass.

The controller's original exit is **zero**. Both frozen batches, all seven runs,
raw profiles/native outputs, seed/final checkpoints, baseline payloads, source
hashes and acknowledged CPU-slot lease/resumption are retained. The archive's
verifier reanalyzes disposable copies of all seven replicas and the baseline,
preserving original records; it checks **98997 raw profile samples**, reconstructs
native comparisons and recomputes timing/throughput ratios. Orbital singular
values are retained measurements with the executable PySCF comparator.

All **1221 payloads** verify by size/digest, with byte-identical repackaging and
a byte-for-byte public fetch. Numerical equivalence here does not qualify the
electronic model or the fitted candidates as certified resonance poles.
Large integral and native intermediate files remain on Sadaharu.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/mpi-throughput
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/mpi-throughput/mpi-throughput.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-mpi-throughput.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

Input recipes are [MPI scaling](../../calibration-sa10-tight-mpi-scaling.json)
and [concurrent throughput](../../calibration-sa10-tight-throughput.json).
The [continuation](../../CONTINUATION.md) records the resumed continuum queue.
