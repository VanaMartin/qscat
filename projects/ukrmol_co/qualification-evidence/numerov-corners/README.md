# Finite cross-anchor / endpoint evidence

The [corner and endpoint study](../../NUMEROV_CORNER_CONTROL.md) passes both
stages: **294/294 case checks**,36/36 reference-control checks, original exits0/0.
It covers14 identities with78 new solves and six unchanged reused rows; no native
or inner-region executable runs. Earlier original rejections remain preserved.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-corners
```

[artifacts.json](artifacts.json) pins public bytes; [manifest.json](manifest.json)
records producer freeze`7ea624347ab36e6d54f962005cea9f20370757a5` and separate
verification support. Extract:

| Payload | Destination |
|---|---|
| `results.tar.gz` | `collected/` |
| `source.tar.gz` | `source/` |
| `inputs.tar.gz` | working root, creating `raw/`, `saved/`, `prior/`, `increment/` |
| `analytic-results.tar.gz` | `analytic/` |

First verify all41 frozen source hashes against `source-manifest.json`.
Then extract `verification-support.tar.gz` over `source/` and verify its two
hashes from `manifest.json`. This adds the separate saved-decision verifier
and updates only its test module with two supplementary tests. **No numerical
producer module is replaced**; the original frozen archive stays intact.

With NumPy/SciPy available, run from the working root:

```bash
PYTHONPATH="$PWD/source" python -m projects.ukrmol_co.numerov_corner_replay \
  "$PWD" --output fresh-saved-replay.json
```

The verifier uses separately expressed K/S/phase comparisons, not the producer
evaluator or propagator. It checks every84 saved rows, disk/report/stream equality,
exact energy bits, complete configuration coverage, references and all decisions.
`analytic-results.tar.gz` preserves the12-solve high-precision endpoint prerequisite,
32 passing checks/exits0/0, original source and sampled resources. The separately
checksum-pinned mpmath wheel enables its analytic tests if needed; extract it to
a separate dependency directory on PYTHONPATH.

[Public-byte verification](../../verification-numerov-corner-public-20261019.json)
passes saved decisions and57 tests from fresh downloads. The prospective numerical
replay selection is each corner's finest400-bohr/order8 row and both stretched1-eV
finest800-bohr/order8 rows: **14 reexecutions**, not all78 producer solves.
Selection is independent of scientific verdict, with a pinned worker and separate
450-s forecast/1200-s cap, eight physical cores/8 GiB.

All200/400/800-bohr finest-step/native comparisons pass, worst K2.15388e-7.
The older h=.0125/800-bohr discrepancies remain preserved; their finer h=.00625
checks now pass. This qualifies the14 identities only. No new unique energy is
covered: original132/594 coverage and462 uncomputed identities remain. No full
reconstruction, absolute native-error bound, complex-pole or production claim.

## Completed public numerical subset

All **14 selected rows** reproduce every stored K and S component exactly,
replay exits **0/0**. The21-file replay archive includes the owner, row outputs,
report and sampled resources; separate owner/report/execution payloads expose
the same bytes. Measured441.291700 s /2229.002709 sampled CPU-seconds /
396,722,176-byte recorded kernel high-water,2200 samples at0.2 s, eight physical
cores/8 GiB. All **15 final payloads** are freshly downloaded and verified after
supplementation. Saved-decision verification covers all84 producer/reused rows;
numerical reexecution covers only the predeclared14-row subset.
