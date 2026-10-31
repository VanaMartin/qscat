# Stretched DZ full-continuum/B1 preparation qualifies

The source/tests/exact owner froze and pushed at `b776f4c` before the31October
17:02:14UTC launch. Completion17:03:45UTC: exactly one four-rank integral call
11.080651s and one B1 CONGEN76.726857s. All original native/scientific/profile/
command/container/owner/service exits0, no timeout/OOM/retry or successor.

Measured target8/3/3/0, no PCOs, continuum51/34/34/22, total59/37/37/22;
990990 raw CSFs,84942 L2 plus2256 target-contracted continuum coefficients =
**87198 contracted states**, eight target roots per sector. Profile90.402706s,
kernel peak8898768896bytes (8.288GiB), CPU123.056646s, allocated peak3381870592bytes.
The read-only30-second monitor records normal progress at17:02:52/17:03:22 and
finished success17:03:52; it makes no scientific or scheduling calls.

Seven fresh public payloads/51 raw archived files/452 resource samples and30
original Git-matched sources pass the [independent stdlib streaming reader](../../verify_stretched_preparation.py).
Full2.67GB CSF and integral bytes are streamed/authenticated, not just private
hash attestations. Native inventories, complete CONGEN group/range/count/unit
markers,64-state contraction arithmetic and the native304-byte CSF header
reconstruct independently. Own64-state property/native-energy/dipole checks pass.
The executed contract digest remains `d00a24a308e3b8da289e631ceba937f7086a61d37e8984ac79cb6da151cf6f4f`.

[Public independent success receipt](../stretched-dz-preparation-receipt/README.md).

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-dz-preparation
python3 -S -m projects.ukrmol_co.verify_stretched_preparation projects/ukrmol_co/qualification-evidence/stretched-dz-preparation/evidence.tar.gz projects/ukrmol_co/qualification-evidence/stretched-dz-preparation/collection.json
```

No CIRMAT-phase audit, native inner/export or electronic-accuracy claim follows.
The [separate conditional TZ release](../../stretched-tz-preparation-reuse-release-20261031.json)
may reuse these geometry-free algebraic CSFs only after exact resolved deck and
all eight target-CSF hashes agree; it still requires own TZ integrals, CIDATA,
properties and orbital frame. No inner/pole/B2/fit/sweep is released by this pass.
