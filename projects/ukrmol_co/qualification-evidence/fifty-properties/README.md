# Explicit forty-state properties on the qualified fifty-component frame

[Result](../../FIFTY_SELECTED_PROPERTIES.md): one DENPROP invocation passes
native command/container/owner **0/0/0**,40 targets /five per sector /676 moments,
37.203min /.43713GiB. No orbital optimization, target diagonalization or scattering
inner solve is repeated. Original seed-preservation rejection/exits1/1/1 stay
separate; physical cross-frame channel matching remains unset.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-properties
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-frame
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-import
```

All **8 new payloads** fresh-download-verify;4 target and5 frame prerequisite
payloads reauthenticate. The two prerequisite pointers name their immutable
bytes. `results.tar.gz` contains **22 generated/owner files**, all pinned by
`collection-sha256.json`, including raw40-state `fort.24`, native output/deck,
intermediates, final and stage-checkpoint reports, source/archive manifest,
commands/exits and resources. All generated files are preserved. Ten unchanged
copied native inputs are authenticated but not duplicated; their byte counts,
digests and exact prerequisite target-archive paths are included in the collection
manifest. The original copies remain untouched on Sadaharu.

Source freeze `96b8ed9` produces the property result. Overlay the separate
`verification-support.tar.gz` onto the producer source to add the independent
reader/tests frozen at `5cbefa0`; verify its two hashes against
`verification-support-sha256.json`. This leaves the12 frozen producer sources
unchanged. The independent parser does not call the producer property reader.

After authenticating the new files and13 original input pins, extract the new
collection into a fresh directory and its `source.tar.gz` into `source/`.
Extract the frame prerequisite collection separately into `frame/`. Then:

```sh
PYTHONPATH=source python -m projects.ukrmol_co.fifty_property_reading --properties outputs/properties/fort.24 --frame-report frame/outputs/coverage/result.json --output fresh-property-reading.json
```

This performs **no native or CI solve**. It reconstructs the40 IDs, exact
per-sector/root map, finite/ordered energies versus independent fixed-frame CI,
676 unique moment records, operator/state/symmetry-index consistency, full
state/rank coverage and ground-dipole convention. Maximum energy discrepancy is
4.91610e-10 Hartree and dipole discrepancy−1.12928e-11a.u., below unchanged
1e-7/1e-5 budgets. The
[public verification record](../../verification-fifty-properties-public-20261021.json)
also authenticates all fetched bytes and original/copied inputs, reconstructs
the exact one-line selection change and checks producer/reader map equality.

A separate same-frame subset observation finds38 shared state labels and620
shared moment keys against the original global-low50 export; all620 values are
exact. It is not a new acceptance gate or physical cross-frame correspondence.
Inventory checks do not independently validate every transition-moment value.
The old `pilot.json` is the intermediate stage checkpoint with a pending verdict;
`result.json` and the owner record contain the final pass. Preserve both.

To reproduce native DENPROP itself, all ten reused input units are available in
the prerequisite target archive. Such a new calculation requires a fresh scope
and output directory; it is not part of this saved reading. No inner solve or
production sweep is released by the property pass.
