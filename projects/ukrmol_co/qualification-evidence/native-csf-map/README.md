# Signed target-CSF correspondence: saved-only evidence

[Outcome](../../NATIVE_CSF_MAP.md): all8 sectors /158,994 CSFs /1,670,879 determinant
terms per layout pass exact signed address correspondence. Only16 singlet-A1 CSFs
change sign. All74 native target vectors round-trip exactly in memory, preserving
energies/root IDs/Gram matrices. **No converted CIDATA or native computation.**
The preferred native phase-aware route retains original CIDATA unchanged;
consumer execution, integral equivalence and scattering-channel maps remain unset.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-csf-map
```

Eight new payloads contain original/reconstruction results, engine provenance,
11 source entries across producer91a4e7f/reconstructione2e126c freezes, exact
input-member pointers and prerequisite artifact indices. `source-sha256.json`
authenticates every archived Python/design/contract entry. No prior CSF/vector
capsule is republished. `input-members.json` names17 exact members of three
authenticated prior archives:

* target-only: `co-fifty-import-20261021/results.tar.gz`, SHA256
  `dc055c3964aa52474715c70ea405185c7f75b8bd5cd1907bdc8e0577010dd389`;
* first full CSF: `co-immutable-continuum-20261021/initial-results.tar.gz`, SHA256
  `271a547ee1be5e343016e965e91599d57fe91f8e0353a2db1b4238ac6f5d5f7d`;
* remaining full CSFs: `co-immutable-continuum-20261021/completion-results.tar.gz`,
  SHA256 `2ad26a64eacf27cb0e746f77f79bc5e6f91cd2633dd93941dbd11398583a3918`.

Fetch the prerequisite indices in separate fresh directories with `qscat-run
fetch`, then extract only the members named in `input-members.json` to their
keys under a fresh `inputs/` directory, verifying hashes. Copy the fresh engine
provenance there. Obtain the UKRmol-in archive from the pinned URL in
`manifest.json`, verify its SHA256 and save it as `inputs/ukrmol-in-3.3.0.tar.gz`.
Extract the six contract-listed native source anchors beneath
`inputs/ukrmol-in-3.3.0/source/`; authenticate their hashes. Extract the new
`source.tar.gz` under a fresh `source/`, checking `source-sha256.json`.

Run from this fresh directory (not a checkout that shadows the archived modules):

```sh
PYTHONPATH=source/producer python -m projects.ukrmol_co.native_csf_verification \
  --root inputs \
  --contract source/producer/projects/ukrmol_co/native-csf-map-contract-20261021.json \
  --output replayed-original.json
PYTHONPATH=source/reconstruction python -m projects.ukrmol_co.native_csf_reconstruction \
  --root inputs \
  --contract source/producer/projects/ukrmol_co/native-csf-map-contract-20261021.json \
  --original result-original.json --output replayed-reconstruction.json
```

Use Python≥3.12 with NumPy/SciPy for the producer; reconstruction is stdlib-only.
Each output must be fresh. Compare each result's `sectors` to the original
published result; elapsed time/RSS need not agree across machines. These commands
read saved files only and cannot invoke native QC/Hamiltonian/property programs.

[Fresh public verification](../../verification-native-csf-map-public-20261021.json)
has already reproduced every sector metric exactly, exits0/0, after authenticating
all8 payloads /3 prerequisite archives /17 inputs /11 source entries /6 native
anchors. Original producer16.496s /0.60803GiB; independent full-expansion
reconstruction38.100s. The physical two-A1 subspace rejection and original
seed-preservation rejection are unchanged. This is not integral-accuracy,
electronic-model or scattering qualification.
