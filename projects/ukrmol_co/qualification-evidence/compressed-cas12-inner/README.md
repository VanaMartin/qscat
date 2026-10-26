# Compressed CAS12 B1/2048 native export

**Own finite export independently qualifies:**2048 selected roots,250 boundary
channels,86352 contracted states and1410 stored continuum components. Recorded
native maximum residual1.34620e-10Ha; CI and boundary energies agree within the
1e-10Ha gate (observed difference0). All99 own native K references decode and
their phase sums agree with native EIGENP within4.97418e-8rad moduloπ, below the
unchanged1e-6rad same-root budget. No S matrix or outgoing pole was constructed
by this reconstruction; the114 protected route holdouts remain untouched.

[Result/resources](../../progress-compressed-cas12-inner-result-20261026.json) ·
[Public reconstruction](../../verification-compressed-cas12-inner-public-20261026.json) ·
[Exact native contract](../../compressed-cas12-inner-contract-20261026.json)

Source/tests/exact owner froze/pushed at `6ae5b79`, before26October04:42:17UTC
launch. Finish08:56:07UTC. Native stages: SCATCI15086.519829s,
RSOLVE124.850008s, EIGENP.015859s, RESON.015592s; all native exits0. Original
scientific/command/container/owner/service/resource exits0, no timeout/OOM/retry.
Profiled container15217.721929s (4.227h), owner15230.767316s;
181161.455029CPU seconds, kernel peak24148611072bytes (22.490147GiB), allocated
output disk2759098368bytes. Twelve physical cores0–11 /64GiB /8GiB shared backing.
Systemd7.8G /17.102CPU seconds are supervisor/hash-cache accounting, not scientific
container resources. Previous8-rank/equilibrium cost is not a controlled speedup
comparison with this12-rank/compressed calculation.

Eight public payloads authenticate100 result/source files and50 prerequisites.
All44 source entries also match frozen Git blobs. The replay independently reads
Fortran CIDATA/channel/boundary/K records, assembles all six arrays, checks native
per-root recorded residuals, metadata/closed gaps, CI-boundary energies and
K-to-EIGENP phase arithmetic. All copied inputs/deck and242 installed engine
artifacts remain unchanged according to the pinned original owner. This is not
an independent Hamiltonian action, CIRMAT phase audit, exact-Lambda proof,
omitted-spectrum bound or electronic-model/curve qualification.

The first local reader exited1 because it assumed two center records as in
target CIDATA. Scattering CIDATA contains C/O plus an uncharged continuum center
at the origin. The corrected reader requires all three, charges6/8/0 and actual
C–O distance1.9bohr. [Correction record](reader-correction.json); initial reader
bytes and failed local replay remain preserved. Corrected local/public replay
exits0, with **no native rerun or scientific-file change**.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-cas12-inner
mkdir cas12-inner-replay
cp projects/ukrmol_co/qualification-evidence/compressed-cas12-inner/artifacts.json cas12-inner-replay/
uv run --no-sync python projects/ukrmol_co/qualification-evidence/compressed-cas12-inner/reconstruct-saved.py \
  cas12-inner-replay/artifacts.json cas12-inner-replay/verification.json
```

The reader streams large payloads into a fresh pointer-only directory using
exclusive creation. No QC, target diagonalization, DENPROP or preparation was
repeated. The original native pipeline includes one RSOLVE propagation to100bohr;
the adapter's zero new-propagations counter means zero new complex/Python maps.
Native RESON fits are diagnostics, **not the coupled outgoing-pole target**.
B1-only Pi pairing/pole/model verdicts remain unset; no B2 or sweep is automatic.
