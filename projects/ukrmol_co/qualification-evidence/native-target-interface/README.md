# Native target consumer and saved integral comparison

[Outcome](../../NATIVE_TARGET_INTERFACE.md): original native CIRMAT reads all74
target roots exactly under the proved CSF sign map; source/full target h and
overlaps are exact, max ERI difference4.44089e-15Hartree, conservative common-basis
operator bound6.26708e-13Hartree. **No QC, target eigensolve, property export or
scattering Hamiltonian.** Original CIDATA remains unchanged.

Three failed attempts preserve compiler relocation, missing static-link
dependency and a quad-header/double-library ABI mismatch, each original
command/container/owner1/1/1. The completed native target read in the third attempt
is saved and never repeated. A separate matching-double-header completion passes
0/0/0, reading only the two integral files. The future Dockerfile copies matching
headers too; no image rebuild/test is claimed by this experiment.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-target-interface
```

The seven payloads include `results.tar.gz` (SHA256
`10c3925dfbe70391ae164d5ae8f73e49c927f2d18217ce617251d78fedef023d`),
all-file collection hashes, input/prerequisite pointers, summary/manifest and a
standalone saved-only reconstruction script. The archive preserves **115 files**:
four owners/source freezes, original exits/logs/resources, old/new compiled I/O
harnesses, completed target/integral diagnostic arrays, empty failed integral
output, recovered double-module capsule/provenance and original inherited quad
module metadata. The source freezes are eb092a3 /87860b8 /c1c3026 /1a785ee,
with16/17/18/22 source entries. Five linked native libraries are hash-pinned.

`collection-sha256.json` lists **six unchanged copied files** not republished:
source/full integral read-copies in two attempts and the final copied completed
target.bin/target.out. They resolve to exact prior input-capsule members or the
third attempt's member. The original14 input pins and16 successor pins resolve
through the four source archives, raw failed-owner files, header capsule and
`input-members.json`; no current remote calculation directory is needed.

For a fresh **stdlib-only, no-native** reconstruction:

```sh
python projects/ukrmol_co/qualification-evidence/native-target-interface/reconstruct-saved.py \
  --directory projects/ukrmol_co/qualification-evidence/native-target-interface
```

The script authenticates all7 payloads, downloads and authenticates3 prerequisite
archives as named in `prerequisite-artifacts.json`, checks all115 files /73 source
entries /30 input pins /46 double modules, independently decodes native records,
flips IEEE sign bits for every coefficient, and reconstructs the scalar integral
bound/symmetries. It writes a fresh `verification-saved.json`; an existing output
is not overwritten. It does not execute the archived native binaries or import
the producer's numerical reader. The prior CSF/vector/integral capsules are
referenced, not republished.

[Fresh verification](../../verification-native-target-interface-public-20261021.json)
already passes this full reconstruction. Native target-channel scattering phases,
root/property maps and the scientific old40/new50 contrast remain **unset**.
Both A1 physical retained subspaces still reject. The bound compares saved
numerical target Hamiltonians, not independent integral accuracy or a scattering
error. No scientific job remains.
