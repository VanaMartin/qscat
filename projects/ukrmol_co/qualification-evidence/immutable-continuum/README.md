# Immutable-frame continuum integrals and target CSFs

[Result](../../IMMUTABLE_CONTINUUM.md): one full-continuum integral call plus
eight target CONGEN calls generate154 MOs (13 target /141 continuum) and eight
complete target-CSF files without another QC, target diagonalization, property
export or scattering Hamiltonian. The original wrapper stops on an incorrect
CONGEN trailer check after two successful native calls; its **1/1/1** exits remain.
A separately frozen seven-sector completion passes **0/0/0** using those saved
integrals and first CSF. No completed native call is repeated.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/immutable-continuum
```

All **6 payloads** fresh-verify. `initial-results.tar.gz` preserves25 files /
44,141,420 uncompressed bytes, including the failed owner/logs,13-file source
freeze `bc1a0c3`, full integral data and completed first CSF. Its SHA256 is
`271a547ee1be5e343016e965e91599d57fe91f8e0353a2db1b4238ac6f5d5f7d`.
`completion-results.tar.gz` preserves32 new files /157,941,080 bytes, including
the15-file corrected source freeze `14c4c61`, seven remaining CSFs/logs/decks,
final/stage reports and resources. SHA256:
`2ad26a64eacf27cb0e746f77f79bc5e6f91cd2633dd93941dbd11398583a3918`.

`collection-sha256.json` pins all57 generated/owner files and lists **eight
unchanged completion copies** that are authenticated but not republished. It
identifies each exact initial-archive member needed to reconstruct them. The
original input capsule contains all17 original input files /381,881 bytes;
the successor's11 additional pins resolve to the initial collection/source.
No prior CI-wavefunction or property computation is needed for saved reading.

Extract both result collections into separate fresh `initial/` and `complete/`
directories. Authenticate all files using the collection manifest. Restore the
eight unchanged copied members into `complete/` using
`unchanged_completion_copies_not_republished`; these point only to exact
authenticated initial members, never to another calculation. Extract each
directory's `source.tar.gz` into its own `source/`, and verify13/15 source hashes.
The original failing source is not overwritten by the corrected successor.

From `complete/`, saved native reading requires no native or PySCF calculation:

```sh
PYTHONPATH=source python - <<'PY'
import json
from pathlib import Path
from projects.ukrmol_co.immutable_continuum import (
    congen_completion, csf_records, native_inventory,
)
w = Path('outputs/preparation')
r = json.loads((w / 'result.json').read_text())
counts = native_inventory((w / 'integrals.native.out').read_text())
for sector, saved in r['target_csf_records'].items():
    spin, irrep = sector.split('.')
    unit = (710 if spin == 'singlet' else 730) + ('A1','B1','B2','A2').index(irrep)
    dimension = {'singlet.A1':15548, 'triplet.A1':24366,
                 'singlet.B1':15196, 'triplet.B1':24540,
                 'singlet.B2':15196, 'triplet.B2':24540,
                 'singlet.A2':15044, 'triplet.A2':24564}[sector]
    congen_completion((w / f'congen.{sector}.out').read_text(), unit, dimension)
    assert csf_records(w / f'fort.{unit}', counts['total'], sector, dimension) == saved
assert r['completed_layout_generation_verdict'] is True
assert r['original_preparation_exits'] == [1,1,1]
assert r['retained_subspace_dominance_verdict'] is False
print('Eight complete target-CSF files; original wrapper and physical rejections retained')
PY
```

[Fresh public verification](../../verification-immutable-continuum-public-20261021.json)
also uses a separate SciPy Fortran-record reader to reconstruct48 record digests,
headers and complete EOF, native inventories/markers, all nine native exits,
exact card/Molden bytes,17 parent inputs and11 successor pins. It performs no
native replay. Integral contents are authenticated opaque bytes, not an
independent integral-accuracy test. All eight new CSFs match recorded host-observed
baseline file hashes; this reading does not independently regenerate/download
the baseline CSF definitions. Native vector signed-basis/address conversion and
Hamiltonian agreement remain **unset**, not supplied by those equal hashes.

Original/completion costs11.530s /3.519s, peaks2.46218GiB /.28282GiB, disjoint in
time on four reserved physical cores. These do not forecast a target reimport or
scattering solve. The two-A1 physical subspace rejection remains unchanged;
no automatic new target solve, property export or Hamiltonian is released.
