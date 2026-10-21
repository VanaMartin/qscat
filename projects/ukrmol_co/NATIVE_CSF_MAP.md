# Exact same-frame target-CSF correspondence

## Outcome — 21 October 2026,18:20 UTC

**The signed target-only → full-layout CSF correspondence passes in all eight
sectors.** The 74 saved target vectors can be represented exactly in the new
layout by a diagonal ±1 map: **16 singlet-A1 CSFs flip sign; all other CSFs are
unchanged**. No CSF/root permutation or rotation is needed. This removes the
signed-basis obstacle without another QC/CASCI, native target diagonalization or
property calculation.

The preferred route is **retain the original CIDATA unchanged**, not generate a
padded/reheadered replacement. The pinned native consumer already implements
the source/destination phase product. Actual native consumer execution and the
eventual scattering target-channel map are **not yet qualified** by this
saved-only result. Integral equivalence is also unset. The physical old40/new50
two-A1 subspace rejection is a different question and remains false.

## What is proved

All **158,994 CSFs /1,670,879 determinant terms per layout** are read to complete
EOF. Each reference and every packed determinant substitution maps exactly by
preserving target orbital (irrep, ordinal, alpha/beta). Only unused continuum
addresses are inserted, so the occupied-address map is increasing and introduces
no extra determinant sorting sign. All term coefficients match byte-exactly
after the per-CSF sign product. Canonical first-determinant phases reconstruct;
orbital metadata, spin populations, counts, offsets and unit CSF norms pass.

| Sector | CSFs | Determinant terms | Negative CSF signs |
|---|---:|---:|---:|
|Singlet A1|15,548|165,784|16|
|Triplet A1|24,366|196,626|0|
|Singlet B1|15,196|198,044|0|
|Triplet B1|24,540|226,343|0|
|Singlet B2|15,196|198,044|0|
|Triplet B2|24,540|233,486|0|
|Singlet A2|15,044|205,288|0|
|Triplet A2|24,564|247,264|0|

For all **74 CIDATA roots**, native integer phase metadata agrees with the source
CSFs; root IDs/order, dimensions, physical centres, energy/diagonal records,
finite coefficients and complete EOF pass. Applying the signs twice restores
the coefficient bytes exactly. Energies and root ordering are untouched; the
maximum vector Gram error is **1.17684e-14**, and mapped Gram change is **zero**.
No overlap normalization or root reidentification is used.

The source target CIDATA remains SHA256
`f81547d9776f6b5f23fd8921d4ba0202ddbb912d069d313e23b629a01527ec59`.
**No converted native file is written.** The map is evaluated in memory and its
sign/vector/energy checksums are recorded. It certifies representations in the
**same new orbital frame**, not physical-state agreement with the old frame.

## Native source interpretation

The downloaded UKRmol-in3.3.0 source archive matches the running image's recorded
source provenance, SHA256
`1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2`.
Six source-file hashes and qualified routine anchors are retained in the
[design/contract](../../docs/physics/co-native-target-address-map.md) and
[machine contract](native-csf-map-contract-20261021.json).

`wrnfto` stores packed substitutions, coefficients and per-CSF phases; `dophz0`
computes the canonical first-determinant phase. `SolutionHandler` writes
unmodified eigenvectors and **integer** phase metadata separately. The native
`Target_RMat_CI.read_ci_mat → CIRMAT → READCIP/CIVIO → CPAKPZ` path multiplies
source/destination CSF phases on target import. `READCIP` reads the CIDATA file's
own centre count, then that many geometry records; it does **not** require the
two-centre target file to acquire the continuum integral file's third centre.

These are source-verified consumer semantics, **not an execution certificate**
for an unrun scattering driver. The full scattering CONGEN `dophz` target-channel
phases, root/property maps and integral target subblock remain to be checked.

## Independent reconstruction and measured cost

Producer freeze **`91a4e7f`** checks saved inputs in **16.496s**, peak
**0.60803GiB** (Darwin process RSS), original command exit0. Independent freeze
**`e2e126c`** takes **38.100s**, exit0. It uses a separate `struct` record reader,
scalar native-address → (irrep,ordinal,spin) decoder, canonical determinant keys,
sorting parities and coefficient dictionaries. It does not import the producer
reader or derive its signs from producer phase metadata. Every complete canonical
CSF expansion gives the same exact sign map and term counts.

Both are one-worker local saved-file checks, not remote native calculations;
CPU affinity and a hard memory bound were not imposed. The prospective180s/2GiB
estimate and300s shell cap are retained, not relabelled measured limits.

[Fresh public reconstruction](verification-native-csf-map-public-20261021.json)
authenticates **8 new payloads /3 prior result archives /17 native input files**,
11 frozen Python-source entries and six native-source anchors. Original and
independent routines run from their distinct archived source trees on those
fresh input bytes, exits0/0,16.599s /37.457s, reproducing every sector metric
exactly. Prior target/full CSF capsules are referenced, **not republished**.
An initial fetch helper used the wrong multi-suffix object naming and returned
404 before any calculation; that partial download remains separate from the
successful fresh reconstruction. No scientific/native failure is repaired here.

Analytic/regression controls cover fermionic slot parity, target address
embedding, negative CSF orientation, a changed nonleading determinant, invalid
phase/occupation/count/EOF, typed CIDATA phases, all74 synthetic root records and
malformed root/energy/centre/trailing data. See the two `test_native_csf_*` files.

## Next production decision

**19:55UTC advance:** [the installed native target consumer and saved integral
target blocks](NATIVE_TARGET_INTERFACE.md) now pass. All74 vectors/relative
energies agree byte-exactly; h/overlaps exact, common-basis Hamiltonian difference
bound6.26708e-13Ha. The completed target read is reused while matched double GBTO
headers repair the separate integral I/O ABI. Actual scattering target-channel
phase/root/property maps and the physical A1 rejection remain separate.

Do **not** repeat QC, target diagonalization or the37-minute property export
merely to change orbital addressing. Retain the original phase-aware target
CIDATA and saved40-state properties. A bounded native **read/import-only** check
and saved target-integral subblock comparison can now qualify the remaining
technical interface without assembling a scattering Hamiltonian.

Before scattering, declare combined orbital/channel prescription sensitivity
versus a separately justified matched-manifold selection, then qualify the
actual target/property/channel phase maps. No matched physical manifold is
created by this sign result; both A1 physical rejections remain. No dense inner
Hamiltonian, root/channel ladder, production sweep or fitting is authorized.
No scientific job is running.

[Public reproduction](qualification-evidence/native-csf-map/README.md) ·
[Dated result/resources](progress-native-csf-map-result-20261021.json)
