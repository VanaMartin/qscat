# Resource-corrected stretched K attribution

## Why a new fixed execution is justified

The [original fixed attribution](co-native-k-attribution.md) produced no K
matrix: both native32 single-energy slices reached its180-s cap. Original
verdictunset, native exits−9/−9 and runner/container exits1/1 are
[published and independently verified](../../projects/ukrmol_co/qualification-evidence/native-k-attribution-timeout/README.md).
The unresolved decision remains the same stretched1-eV K discrepancy.

The native cost estimate incorrectly amortized a setup over energy count.
Checksum-pinned source establishes that `ASYM1 → RPROP1 → SETMTR/AMPLTD` runs
before the energy loop and K-file initialization. Each radial subrange
diagonalizes a2980-dimensional Hamiltonian, even for one incident energy.
`RPROP1` uses fixed MRANGE directly; the energy range only selects its count
when MRANGE=0. Thus slicing the input does not avoid its32/64 setup. Partial
buffered logs do not locate the live routine when killed; no stack was sampled.

## Cost model and changed execution resources

The prior full-grid stretched native32 calls cost866.028/920.267 s. Conservatively
retain those totals as one-energy setup forecasts, with a linear subrange
forecast for64. For two sectors this gives:

- native32:1786.295 CPU-seconds;
- native64:3572.591 CPU-seconds, conditional on linear subrange cost;
- ten independent solves:about1510 CPU-seconds from398720 radial intervals;
- total:about6869 CPU-seconds; the two native bundles run concurrently, with
  longest-bundle forecast2760.8 s, followed by at most600 s independent work.

Use **3400 s expected /5400 s upper /9000 s hard batch cap** and a **3600-s
per-native-stage cap**. These are conservative forecasts, not measured one-energy
or64-subrange timings. Four single-threaded workers on physical cores0–3 and
the4-GiB container limit remain. The native64 surface arrays scale linearly
with subranges; two native32 processes previously occupied0.84241 GiB in total,
so doubling their arrays is consistent with this cap, without proving sufficiency.

The [resource-successor contract](../../projects/ukrmol_co/native-k-attribution-resource-contract.json)
is the original scientific contract plus only the changed status/resource
allowance and checksum-pinned timeout/public/source-trace prerequisites. Compare
all other fields exactly before launch. Freeze and push this contract before
execution in a new persistent directory, using the unchanged producer source.
Verify engine, source, inputs, prior outputs, resource/workload state and the
published timeout prerequisites before launch and afterward.

## Gates and stopping

There are exactly four possible native invocations and ten independent controls.
The native32 slice must reproduce saved32 at1e-10 with identical energy bits;
native32→64 and each independent step/endpoint/order/reference comparison retain
their separate1e-6 K/S/phase gates. An unset/false native gate stops submissions;
the fixed in-flight bundles finish. No added rungs, further retries, relaxed
numerical budgets, inner solve or automatic full-grid continuation follows.

The original180-s attempt stays unset/exits1/1. The original six-case rejection
stays false/exits1/1. A passing resource-successor point study would qualify
only its two identities and controls, not an absolute error or a CO pole.

## Recorded outcome

The [completed experiment](../../projects/ukrmol_co/NATIVE_K_ATTRIBUTION.md)
finishes in1669.647 s /1.973831 GiB with **false / original exits1/1**. Native32
slices reproduce the saved matrices exactly,32→64 passes,and all ten independent
controls finish with passing invariants. Both sectors reject the separate step,
endpoint and finest/native K gates. A smaller step worsens the discrepancy;
tail-order sensitivity is tiny. This rejects a simple monotone step cure and
requires separately validated precision/algebra or independent-propagation
attribution before further native work. No successor follows automatically.
