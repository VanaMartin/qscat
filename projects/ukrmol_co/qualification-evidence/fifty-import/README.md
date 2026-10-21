# Fifty-component native target import evidence

The [result](../../FIFTY_COMPONENT_IMPORT.md) retains a **target/analyzer pass
and seed-frame preservation rejection**. Original target command/container
exits0/0; analyzer0/0; independent verifier1/1 and owner1. Native energies match
the newly optimized frame, while six averaged QC roots move more than1e-7
Hartree from the independently covered seed. No tolerance changes or retries.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-import
```

All **four public payloads** fresh-download-verify. Authenticate against
`artifacts.json` before extraction. `results.tar.gz` contains153 original files,
all pinned by `collection-sha256.json`:

- `run/`: the entire completed target run, including native binaries, checkpoints,
  templates, config, QC/engine logs, result and original resource records.
- `owner/`: unchanged failed owner, execution record, container inspections,
  launch, stage/resource logs and12-file source archive/manifest.
- `reference/`: original QC checkpoint/Molden, target/config/result/resources.
- `coverage/result.json`: original two-frame96-probe coverage,48 for this seed.

`result-summary.json` distinguishes current-frame native agreement from rejected
seed preservation. [Fresh public-byte verification](../../verification-fifty-import-public-20261021.json)
recomputes the50 averaged native-table comparisons and six seed-frame energy
changes. The overall import verdict remains **unset**; the original verifier
never reaches all74 binary-vector/EOF or subspace checks. No numerical rerun is
part of this saved-file reconstruction.

For the original verifier's expected failure, extract the collection into a
fresh root; put `run/` under `runs/co-eq-ccdz-sa11-fifty-slepc-target-mpi8-20261021/`
and `reference/` under `runs/co-eq-ccdz-sa11-50-qc-hostprep/`. Place the coverage
file under `diagnostics/equilibrium-fifty-component-ci-residual-coverage-20261007/`.
Overlay `owner/source.tar.gz` into `source/` and verify its12 hashes. In the pinned
Linux x86-64 image from the contract, mount the root at `/work`, source read-only
at `/opt/qmodeling`, set `PYTHONPATH=/opt/qmodeling` and use the original command:

```sh
python3 -m projects.ukrmol_co.fifty_import_verification --root /work --contract /opt/qmodeling/projects/ukrmol_co/fifty-import-contract-20261021.json --output /work/fresh-verification.json
```

Expected outcome: exit1 on the seed-frame comparison; no successful verification
JSON. This only reads saved data. Source freeze `7827867`, original image,
complete commands and original exits are in `owner/execution.json`. Reproducing
the three-hour native/QC production itself requires a new run directory and a
separate authorized execution; it is not needed to reproduce this rejection.
