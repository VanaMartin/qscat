# Independent coverage of the completed fifty-component frame

## Completed — 21 October,14:43 UTC

**The new fixed orbital frame passes all48 probes and all74 saved native
roots/vectors/header/EOF checks**, original command/container/owner exits0/0/0.
Maximum native/reference energy error **8.24230e-13 Hartree**, native-vector
Gram error **1.26565e-14**. Maximum independent/fresh-QC ensemble-root difference
**1.42109e-13 Hartree**. Physical/penalized residual maxima **7.20000e-10 /
9.94442e-10 Hartree**, penalty-action maximum6.88784e-10; all remain below the
unchanged1e-9 gates. These are finite residual/coverage checks, not an absolute
electronic-model error bound.

Checkpoint/Molden physical-overlap minimum singular values are
**0.9999999999999499 core /0.99999999999998 active**, above1−1e-8. The mandatory
Hubbard-dimer exact/perturbed controls pass. Local focused tests:11 pass, one
PySCF skip; its analytic kernel is exercised successfully in the remote run.
The unchanged reference helper retains exact SHA256 `095e81c…` below.

Measured runtime **308.657s /5.144min**, sampled CPU979.610s /16.327CPUmin,
peak **290,611,200 bytes /.27065GiB**,1540 samples; cores0–3 /four threads /4GiB
cap. All10 frozen sources and10 original input pins verify before/after. Source
freeze **`55e89a6`** precedes execution. No orbital optimization, native engine
invocation or scattering solve is repeated.

[Public evidence](qualification-evidence/fifty-frame/README.md) fresh-verifies
all5 new payloads /108 original files and reauthenticates4 prerequisite payloads.
The [saved reconstruction](verification-fifty-frame-public-20261021.json)
recomputes all74 native vector Grams and headers/EOF, energies versus the saved
independent CI references,48 declared probe decisions and the retained original
seed rejection. CI residual actions and checkpoint/Molden overlaps are producer
observations, not repeated physical solves in the saved-byte reading.
[Dated result](progress-fifty-frame-result-20261021.json) preserves the new pass
beside the original six-root seed rejection/exits1/1/1.

The [explicit selected40-property calculation](FIFTY_SELECTED_PROPERTIES.md)
subsequently starts14:52UTC under a separately frozen contract `96b8ed9`. It
addresses the global-low50 property's missing singlet-A2 roots4–5, reusing these
qualified vectors. Physical cross-frame state identity and matched scattering
remain open.

## Prospective scope retained

This successor qualifies the **new fixed checkpoint**, not preservation of the
old seed. The original six-root seed comparison remains rejected with original
verifier/container/owner exits1/1/1. Its source, bytes and1e-7-Hartree gate stay
unchanged in the [published target evidence](qualification-evidence/fifty-import/README.md).

## Finite question and independent reference

The [contract](fifty-frame-coverage-contract-20261021.json) pins the completed
checkpoint, Molden, binary target/property exports, QC/config/results/resources
and original failure records. No orbital optimizer, native engine or scattering
Hamiltonian is invoked. Previous coverage from the old checkpoint is not reused
as a qualification of the new one.

The CI/residual helper is byte-identical to the independently validated earlier
two-frame96-probe implementation, SHA256
`095e81ca0b57ded9a5eb6c6f082860864fdd641bd72eda247cd696ea6e2f9cfd`, now tracked as
`coverage_reference.py`. The analytic Hubbard-dimer benchmark and perturbed-vector
negative control run before the new probes. The new reference differs from the
original orbital optimizer's CI by independently constructing fixed-frame CASCI
solves, root counts/trial spaces and physical/penalized residual actions; native
eigenvectors/eigenvalues come from UKRmol's separate CSF representation.

All **48 fixed-orbital probes** execute: both spins /four irreps, roots10/13 in
A1 and5/8 otherwise, spaces40/80/160, each with the already qualified600-cycle
fixed-orbital cap. The ordinary200-cycle orbital-optimization contract is not
changed. CI energy/lindep tolerances remain1e-12/1e-22. Physical and penalized
residuals and penalty-action norms must each be≤1e-9 Hartree, spin error≤1e-6,
CI-vector Gram error≤1e-9, and energies across spaces/root counts and versus the
new50-component QC record≤1e-7 Hartree.

The new checkpoint is compared to its own imported Molden through the physical
AO cross metric: both frozen-core and active-subspace singular values must be
≥1−1e-8. After coverage passes, read every vector/header and complete EOF of all
**74 saved native roots**, with energy gate1e-7 Hartree and Gram gate1e-9. These
are the original binary checks factored into a shared helper without changing
their conditions. The original verifier retains its seed comparison and still
rejects that comparison; the successor has a separately named verdict.

## A pass does not yet give matching scattering channels

The original global-low50 DENPROP export contains per-spin-sector counts:

| |A1|B1|B2|A2|
|---|---:|---:|---:|---:|
| Singlet |7|6|6|**3**|
| Triplet |6|8|8|6|

Thus it lacks singlet-A2 roots4–5 required by the existing five-A2 comparison.
All74 eigenvectors include those roots, but missing transition properties cannot
be manufactured by relabeling IDs. A separately frozen explicitly selected
property export can reuse the native target vectors and integrals without
reoptimizing or rerunning an inner Hamiltonian. Its cost and physical subspace/
root correspondence must be established before matched scattering. This successor
reports the known property inventory and makes no matched-channel claim.

## Cost and stops

One sequential four-thread worker on reserved physical cores0–3 /4GiB. Prior
two-frame96 probes measured1280.58s /.658GiB. Forecast **900s expected /1800s
upper /3600s hard cap**; no attempt is made to repeat the10436-second QC/import.
Require8GiB available RAM and30GiB disk, persistent outputs, read-only completed
inputs, before/after ten-input and frozen-source integrity checks, all logs,
resource samples and original exits. Freeze/push source, tests and contract before
execution. A failing gate stops with no repair ladder or tolerance change.

A passing new-frame result releases matched-property/channel-map and cost design,
not a production sweep. The larger electronic/basis/channel, extraction,
continuity, neutral-domain and restartable-pilot gates remain independent.
