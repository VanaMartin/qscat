# Rejected interior analytic prerequisite: public evidence

The [interior stress prerequisite](../../NUMEROV_INTERIOR_CONTROL.md) completes
24 analytic rows with **64/68 checks passing**, aggregate **false / original
runner-container exits1/1**. Its stopping rule withholds the conditional CO run.
No new CO/native/inner solve or unique-energy coverage is produced.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-interior
```

[artifacts.json](artifacts.json) pins public payloads;
[manifest.json](manifest.json) preserves producer freeze`b7f2278` and separate
verification support. Extract:

| Payload | Destination |
|---|---|
| `source.tar.gz` | `source/` |
| `results.tar.gz` | `collected/` |
| `selection-references.tar.gz` | working root, creating `saved/` |

Verify all46 frozen source hashes against `source-manifest.json`. Extract
`verification-support.tar.gz` over `source/`; verify its two hashes from
`manifest.json`. It **adds** two saved-decision verifier modules, without replacing
any frozen module. Extract the pinned mpmath1.3.0 wheel into `dependencies/`.

With NumPy/SciPy installed, run from the working root:

```bash
PYTHONPATH="$PWD/dependencies:$PWD/source" python -m projects.ukrmol_co.numerov_interior_replay \
  "$PWD" --output fresh-saved-replay.json
```

The verifier checks native reference/history hashes, deterministic identity
selection, all24 analytic coordinates, exact energy bits, continuous and scalar
discrete80/110-digit oracles, and all68 saved decisions. It shares model/oracle
definitions with the producer and uses separately expressed scattering/error
calculations; it does not reexecute propagation. Saved-decision verification
exit0 reproduces scientific **false/exits1/1**. All62 focused tests pass from fresh
public bytes. Full numerical replay uses the separately pinned worker/contract,
all24 identical analytic rows, four physical cores0–3/4 GiB/600-s cap.

Four rejected checks are scalar finest800 floating-point K drift1.20893e-7>1e-7,
and rotated finest200/400/order6 K errors1.11–1.13e-6>1e-6. S/phase/invariant
passes and the earlier14-identity corner pass do not waive the new prerequisite.
Historical coverage132/594,462 uncomputed identities and14 increment-qualified
identities remain. No absolute native-error or complex-pole claim follows.

## Completed full numerical replay

All **24 analytic rows** reproduce every stored K component exactly, with all68
check decisions and the false scientific verdict identical. Replay exits0/0
preserve original scientific exits1/1. All **16 final payloads** fresh-verify.
Extract `public-replay.tar.gz` separately; its8-file manifest pins the complete
reexecuted report, comparison report, owner, exits and sampled resources. The
standalone `public-owner.py`, `public-execution.json` and `public-replay.json`
are byte-identical to their archived copies.

Measured72.218529 s /71.695078 sampled CPU-seconds /44,679,168-byte kernel
high-water,359 samples at0.2 s, four physical cores0–3/4 GiB/600-s cap. No new
CO/native/inner solve occurs. [Verification](../../verification-numerov-interior-public-20261019.json)
records fresh payload hashes,62 tests, saved decisions and full numerical replay.
