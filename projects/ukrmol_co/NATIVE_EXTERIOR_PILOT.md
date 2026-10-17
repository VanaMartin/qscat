# Independent native exterior pilot — launched 17 October 2026

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
The original forty-component orbital averaging and selected48-channel target
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
