# Saved scattering input map — 21 October 2026,21:55 UTC

**The saved B1/B2 phase-algebra and forty-state property-ID audit passes.**
All **317,988 phase-defining first determinants** reproduce the stored target
phases exactly. The saved CONGEN definitions are reusable candidates without
another CONGEN generation: their inventory, prototype grouping and target-sector
dimensions match the qualified full layout. This is not a new physical
cross-frame eigenstate identity result or a full N+1 determinant audit.

## Production-relevant findings

- **Do not copy the old `idtarg`.** Exactly **25/40 entries change** for the new
  lowest-five-per-sector property file. The replacement is recorded verbatim in
  [the result](progress-scattering-input-map-result-20261021.json), not installed
  into production cards. All40 property IDs map to the original CIDATA roots1–5;
  maximum absolute-energy difference **4.91639e-10 Hartree**, gate1e-7, rtol0.
- The unchanged new property file has **676** finite unique moment records;
  the old file has **679**. Both retain complete state/rank coverage and valid
  moment indices. Neither is padded, relabelled or recomputed. This checks
  selection and indexing, not independent accuracy of every transition moment.
- **Scattering phases are not target-only phases.** Continuum placement and
  spin coupling create sector-dependent signs. B1/B2 destination phases also
  differ. Use the native source/destination phase product with original CIDATA,
  not the previous target-only sign vector or a converted CIDATA.
- CONGEN stores **two prototype continua per sector**; SCATCI expands to
  **1,410 channel coefficients** in each symmetry. Do not confuse prototype
  counts with the physical continuum inventory or channel convergence.

| Destination | Expanded continuum counts, interleaved spin sectors |
|---|---|
|B1|34,34,51,51,22,22,34,34|
|B2|34,34,22,22,51,51,34,34|

Each saved CONGEN file has344,124 CSFs, including26,136 L2 CSFs. All six native
records reach EOF, all CSF count/index endpoints and monotonicity checks pass.
The phase audit reconstructs one first determinant per target CSF, not all
6.2 million determinant terms or any scattering Hamiltonian.

## Preserved execution and cost

| Freeze | Saved-only outcome | Original command/container/owner | Wall |
|---|---|---|---:|
|`1899892`|B1 witness complete, then integer-card parser absorbs digit in following `nob0`; B2 unrun|1/1/1|3.52888s|
|`07a64e8`|Parser-only correction; all phase/card/property gates pass|0/0/0|6.83320s|

The original failure is not overwritten or retroactively passed. The successor
inherits all19 original input pins and adds four predecessor pins. All23 verify
before/after; all seven frozen source entries remain unchanged. Native calls,
QC, integral generation, CONGEN generation, target diagonalizations, property
exports and scattering solves are **zero**. Successor Python wall6.45214s,
Linux peak RSS737,088KiB (**.70294GiB**); this is process RSS, not sampled cgroup
memory. Four physical cores0–3 reserved, one process,4GiB limit. Two-job wall
10.36208s; no aggregate CPU-time measurement is claimed.

Persistent Sadaharu evidence:
`prepared/scattering-input-map{,-parser}-source-20261021/`; publication mirror
`prepared/publication-scattering-input-map-20261021/`. No scientific job remains.

## Evidence and its boundary

[Fresh public compact reconstruction](verification-scattering-input-map-public-20261021.json)
authenticates **6 payloads /50 archived files**, plus a fresh prerequisite archive
containing original CIDATA. Its independent permutation-cycle parity algorithm
reproduces all317,988 phases, their raw destination phase records, sector hashes
and the25 changed IDs; absolute property energies are checked against freshly
authenticated CIDATA. It imports no producer numerical reader and runs no native
program. [Fetch/replay instructions](qualification-evidence/scattering-input-map/README.md).

The1.4GB original scattering CSF files remain persistent and their complete
SHA256s were checked during the remote audit. The compact public capsule contains
derived first-term witnesses and exact header/phase slices, **not those complete
CSF files**. Compact replay does not independently reauthenticate their full
hashes or validate all N+1 determinants. Do not describe it as a full raw-CSF
public reconstruction.

## Remaining finite interface check

The earlier target-only native read does **not** qualify the different scattering
destination phases. The smallest remaining native step is two read-only
74-root CIRMAT consumers, one per B1/B2 destination, with their own frozen
contract and matching double headers. No new CONGEN, QC, target diagonalization
or property export is justified by the saved map.

Downstream source inspection locates the consumer in pinned UKRmol-out3.3.0.1,
`source/libouter/swinterf.f:SWTARG`, file SHA256
`303f6fa5867af6c28e91f4489c54609d4932ec6215ddadc687750c63af50b952`:
it matches state IDs through `IDTARG`, sorts energies with `sort_outer`, remaps
moment indices through that sorting permutation, and loads both symmetric
moment entries with `chargesign`. `ORDERC` maps energy-ordered channels back to
CI order; optional `construct_channel_phases` maps target phases using IDTARG.
This trace is **not native execution qualification**, especially for energy
ties/order and downstream moment/channel phases. The source came from retained
build image `6f080f2f678b`, with the already pinned outer archive provenance.

Both physical A1 dominance screens remain rejected. A future lowest-five
old40/new50 comparison is **combined orbital/channel-prescription sensitivity**,
not a clean orbital-average effect. A matched-manifold alternative needs separate
physical projection/selection evidence. Only after consumer semantics and the
scientific contrast are fixed should a sparse2048-root comparison be costed.
No dense Hamiltonian, root/channel ladder, production sweep or potential fitting
is released by this saved-only result.

[Algorithm and contract](../../docs/physics/co-scattering-input-map.md)
