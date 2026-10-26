# Geometry-free CAS12 preparation cards

Copied byte-for-byte from the retained completed CAS12 dimension preflight:
`diagnostics/cas12-scattering-dimension-preflight-nbmx/co-eq-ccdz-sa12-slepc-target/output/CO/geom1/inputs/`
under the campaign root recorded in `CONTINUATION.md`. SHA256 pins are in
`../compressed-cas12-preparation-contract-20261026.json`.

These are **algebraic recipes, not equilibrium physical data**. The integral
card reads geometry/orbitals exclusively from `./co.molden`, supplied by the
hashed compressed CAS12 run. Its radius18bohr, lmax4, deletion1e-6 and retained
target/core counts8/3/3/0 match the control. CONGEN defines the same15-electron
doublet groups/core-active shell algebra and native contraction order. Only
full `nob` addressing may update after the actual continuum inventory is measured.
The retained SCATCI cards are not executed by this preparation; only `notgt`
may resolve to actual radial counts, to count contracted dimensions.

Do not copy equilibrium integrals, vectors, properties, target IDs or eigenpairs.
These cards alone neither qualify target-CSF phases/CIRMAT contraction nor a
scattering interface. All physical inputs come from the qualified compressed
checkpoint/Molden/64 vectors/new40 properties. One integral and two CONGEN calls
are the entire native scope; the later inner solve needs its own exact freeze.
