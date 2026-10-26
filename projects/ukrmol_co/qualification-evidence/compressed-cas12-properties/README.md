# Compressed CAS12 matched first-five properties

**Qualified missing input, not electronic-model accuracy:** one DENPROP export
from the saved64 native target vectors produces40 energy-ordered state IDs,
exactly five in every singlet/triplet C2v sector, and725 unique finite multipole
records. Native target energies agree within4.76760e-10Ha; comparison with saved
qualified QC energies gives4.42394e-9Ha. Ground dipole0.177380573649au agrees
within2.07486e-9au. All original/copy orbital/CI/integral bytes remain unchanged.

[Compact result](../../progress-compressed-cas12-properties-result-20261026.json) ·
[Independent reconstruction](../../verification-compressed-cas12-properties-public-20261026.json) ·
[Property contract](../../compressed-cas12-matched-properties-contract-20261025.json) ·
[Matched-model rationale](../../../../docs/physics/co-compressed-cas12-model-control.md)

Source/tests/exact owner froze/pushed at `bfa3a88` before launch25October23:40UTC.
Finish26October04:06:02UTC; original native/scientific/command/container/owner/
service exits all0. No timeout/OOM or scientific retry. One DENPROP call only;
zero QC, target diagonalization/reimport, integrals, CONGEN, scattering or S calls.

Measured native15946.447569s, profiled container15948.199065s, owner15948.793576s;
container CPU15979.841532s and kernel/sample peak2043465728bytes =1.903126GiB.
Allocation8 physical cores4–11 /one library thread /64GiB /8h native cap. This is
roughly one-core execution, not an eight-core speedup. Systemd's23.9M peak /
7.126CPU seconds describe only its supervisor, not the scientific process.

Fresh stdlib-only public-byte reconstruction authenticates nine payloads,
63 archived files,23 source entries,30 original inputs and12 unchanged copies.
All23 source digests also match the frozen Git blobs. It reads native Fortran
CIDATA records independently, reconstructs state selection/map, all moment
references/ranks, energies/dipole, one changed deck assignment, original exits
and sampled/kernel resource arithmetic. It does not recompute transition moments
or establish electronic accuracy. The own contraction-to-property-ID map is a
permutation of1…40, not CAS11's sparse five-from48 ID subset.

The first public reader exited1: the original-input archive contained a
`fort.16 -> moints` symlink without its target. A local hard-link-only parser
control also exited1. Both failures are preserved; the same30 pinned original
inputs were recollected as regular byte payloads, with unchanged digests.
[Correction record](reader-correction.json). Original native outputs/exits and
first archive/collection remain published unchanged. Corrected fresh public
reconstruction exits0; **DENPROP was not repeated**.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-cas12-properties
uv run --no-sync python projects/ukrmol_co/qualification-evidence/compressed-cas12-properties/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/compressed-cas12-properties/artifacts.json \
  verification-compressed-cas12-properties-fresh.json
```

**Next:** one separately frozen same-frame full-continuum integral/two-CONGEN
preparation. No inner Hamiltonian, pole, geometry, fit or sweep is released by
this property result. Retain nine as working reference; CAS11/five remains only
a local matched comparator. All114 holdouts/equilibrium exclusion stay protected.
