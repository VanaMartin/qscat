# Independent native exterior pilot rejects — 17 October 2026

**The frozen21-case pilot rejects: eight of ten checks pass, but native S and
phase reproduction exceed1e-6 at5eV.** All21 cases complete without execution
failure. Original scientific command/container exits **1/1** remain. There is no
retry, higher-order/endpoint successor or relaxed budget. See
[verification](verification-native-exterior-pilot-20261017.json) and
[public evidence](qualification-evidence/native-exterior-pilot/README.md).

## Completed result — 13:14:24 UTC

| Requested energy / eV | Finest K-entry error | Finest S-entry error | Phase error / rad moduloπ |
|---:|---:|---:|---:|
| 0.1 | 2.95023e-9 | 5.90030e-9 | 1.03832e-8 |
| 2.5 | 2.44019e-8 | 1.14479e-8 | 3.29286e-8 |
| 5.0 | 7.46593e-7 | **1.49310e-6** | **2.72174e-6** |

The K-entry budget passes at all three energies. It cannot replace the separate
S/phase budgets: the Cayley transform changes entrywise sensitivities, and the
phase sums responses across six open channels. Finest unitarity/reciprocity
errors are≤1.844e-11/2.172e-11. Those algebraic checks do not settle agreement
with the native result.

A [separate saved-matrix algebra audit](verification-native-exterior-pilot-algebra-20261017.json)
recomputes every case's metrics using the right Cayley form and a separate native
eigenvalue/arctangent phase path. Maximum reported-vs-recomputed metric drift is
7.78e-16; all ten gate decisions reproduce. This checks the report's algebra,
not the propagation error.

All three independent S refinement families pass their1e-6 budget:

| Family | First change | Final change |
|---|---:|---:|
| Step0.1→0.05→0.025bohr | 1.2130121887e-5 | 7.5777337175e-7 |
| Endpoint100→200→400bohr | 1.3462577919e-8 | 2.6819140924e-8 |
| Tail order4→6→8 | 6.0206521230e-9 | 1.4953838483e-12 |

Spatial improvement is16.01-fold. At5eV, the finest phase mismatch remains
2.602e-6/2.642e-6/2.722e-6 across endpoints100/200/400, and2.722e-6 across
orders4/6/8. The saved controls therefore isolate a persistent native-comparison
discrepancy; they do not assign it to the native code, the independent method,
finite-radius asymptotics or any single normalization convention. A successive
change is not an absolute error bound. Further resolution is not automatic.

The method's17 focused tests and independent public toy replay remain passing.
[Fresh public native-pilot replay](verification-native-exterior-pilot-public-20261017.json)
verifies nine payloads,eleven sources,four native inputs,all21 case records and
all ten check verdicts. It reproduces **false and exit1**, with maximum host/public
K-entry/S-component drift1.61643e-9/3.23265e-9. This verifies the rejection;
full native K and CO poles remain unreleased.

Measured Sadaharu cost is **735.161s /0.08564GiB sampled kernel high-water**
(91,959,296bytes),3671 samples at0.2s; sampled CPU usage734.512s. Sampling is
not a post-exit peak. This is below the975.5s estimate and2400s cap. The
21-line JSONL journal exactly matches all case records in the final report.
Source archive SHA256:
`65117f1aeab98eed68a8bb3ced3371e288c8f4e93524bc44730ac6b5286897d8`.

**Next decision:** use the measured5-eV S/phase discrepancy to justify a separate
native-outer/method error-attribution contract. Source defaultsIASY20/EPS5e-5
are provenance, not proof that either side is accurate to1e-6. Freeze any
independent-method or native outer-only control, its cost and unchanged
comparison budgets before execution. Existing exports suffice; a new inner
solve is not justified. The full594-matrix reconstruction remains unreleased.
The separately [designed subrange-attribution control](../../docs/physics/co-native-exterior-attribution.md)
is frozen at **`666c056`**, with unchanged comparison budgets and a1667s cost
estimate/3600s cap. Its three native outer-only invocations and two independent
cases test native radial resolution; they do not replace this original result.

## Decision and fixed scope

The [combined toy resolution gate](COMBINED_EXTERIOR_RESOLUTION.md) passes;
the next unresolved question is whether independent exterior propagation of the
saved surface/threshold/multipole inputs reproduces **all open-open K entries**,
rather than only a phase sum. The [native sign/flux source trace and numerical
method](qualification-evidence/native-exterior-method/README.md) now pass
independent analytic controls and fresh public replay.

The [pilot contract](native-exterior-pilot-contract.json) is frozen at
**`00c5839` before native comparison**. The method is frozen at **`0a9b0ce`**.
It uses **equilibrium/nine/B1**, all298 channels, all2048 uncorrected exported
poles, and requested energies **0.1/2.5/5.0eV** (native indices0/48/98).
The original forty-component orbital averaging and selected48-component target
retention remain inputs; no new electronic or inner-region calculation occurs.

Seven configurations per energy vary **step0.1/0.05/0.025bohr**,
**endpoint100/200/400bohr** and **tail order4/6/8** independently around one
shared finest point, for **21 cases**. Full native K/S entry errors and phase
moduloπ, real-axis unitarity and reciprocity use **1e-6** budgets. Each final
independent S refinement must be≤1e-6; spatial refinement must improve by at
least eightfold, subject to a1e-8 roundoff floor. Matrices are not symmetrized.

The source-derived convention fixes native sin/cos signs,1/sqrt(k) flux,
total-Hartree surface denominators and already-Rydberg multipoles. All closed
channels are retained. The finite inverse-power tail is not assumed converged:
nearby distinct closed thresholds may amplify its coefficients. The metadata
inventory has33 distinct thresholds, exact degeneracies, and minimum distinct
gap0.000988054Ryd. The independent endpoint/order families test this explicitly.
Native Gailitis defaultsIASY20/EPS5e-5 are preserved as provenance, not a bound
on native error.

## Resources and stopping

The prospectively measured synthetic298-channel/100-interval kernel costs
0.479163s. The pilot has203,580 radial intervals, giving **975.5s estimated**,
**2000s upper estimate**, **2400s hard cap** and **4-GiB memory limit** on
Sadaharu, physical cores0–3 reserved and one BLAS thread. Inputs are mounted
read-only from the existing saved-matrix collection. Each completed or failed
case is flushed to a persistent JSONL journal; original exits and sampled
resources are retained.

One fixed matrix is authorized. False/unset stops without a retry, larger
endpoint, higher order or relaxed budget. Even a pass qualifies only this
three-energy/one-sector pilot, not the six-case594-matrix reconstruction, native
complex poles, general nonzero-Bloch paths, spectral completeness or production.

## Launch observation — 13:05:57 UTC

The [dated host capture](progress-native-exterior-launch-20261017.json) verifies
**33 new immutable remote hashes /238 prior hashes /242 installed artifacts**,
the read-only native mount, physical cores0–3 and4-GiB cap. The container starts
at **13:02:09 UTC**. Six cases have completed with no execution failure, all
at0.1eV. Their K/S errors are≤4.416e-9/8.832e-9, but the full pilot verdict is
unset until all three energies and independent refinement families complete.
At endpoint100/order8 the raw/equilibrated matching conditions reach
3.982e9/7.645e8. Arbitrary closed-basis conditioning is recorded separately
from the entrywise observable comparison.

The unrelated `comfyui` container remains active; free persistent disk is
40.78GiB. Persistent source/output directory:
`prepared/native-exterior-pilot-source-20261017/`. No native engine solve,
observer service or automatic successor is launched.
