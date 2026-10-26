# Compressed CAS12 continuum and scattering-CSF preparation

**Preparation qualifies in both sectors:** target/core8/3/3/0, zero PCOs,
continuum51/34/34/22 and total59/37/37/22. Each sector has990990 raw CSFs;
84942 uncontracted L² CSFs +1410 target-contracted continuum components give
**86352 actual contracted states**. These are measured compressed results,
not substituted equilibrium planning dimensions.

[Compact result](../../progress-compressed-cas12-preparation-result-20261026.json) ·
[Independent public reconstruction](../../verification-compressed-cas12-preparation-public-20261026.json) ·
[Exact execution contract](../../compressed-cas12-preparation-contract-20261026.json)

Source/tests/exact owner froze/pushed at `6bad622` before26October04:21:57UTC
launch; finish04:24:48UTC. Exactly one four-rank integral call10.632730s and
two serial CONGEN calls77.360604/76.412469s; all three native exits0. Original
scientific/command/container/owner/service/resource exits0, no timeout/OOM/retry.
No QC, target solve/reimport, DENPROP, inner solve, S or complex propagation.

Measured owner171.419740s /profiled container170.614204s /202.119108CPU seconds;
kernel peak11489054720bytes =10.700016GiB, allocated output disk6031511552bytes
=5.617283GiB. Four physical cores0–3 /16GiB /one-hour cap. Supervisor metrics
are not scientific-container metrics. All36 input pins and copied original
target/frame/property bytes remain unchanged. All34 source pins match frozen Git.

Six public payloads authenticate83 archived result/source files and36 prerequisites.
Streaming independent reconstruction checks native orbital inventories, final
CONGEN group order/L2 ranges/contraction counts/projection/write markers and
CSF header spin/orbital/electron/dimension fields. Native CSF payloads are fully
hashed without holding multi-GiB arrays in memory. This is **not a CIRMAT phase
audit, native boundary/S-interface qualification or electronic-accuracy result**.
The two generated algebraic CSF files happen to be byte-identical to the saved
equilibrium recipes; physical integrals/Molden/target vectors/properties are
the compressed model's own bytes, never equilibrium substitutes.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-cas12-preparation
mkdir cas12-preparation-replay
cp projects/ukrmol_co/qualification-evidence/compressed-cas12-preparation/artifacts.json cas12-preparation-replay/
uv run --no-sync python projects/ukrmol_co/qualification-evidence/compressed-cas12-preparation/reconstruct-saved.py \
  cas12-preparation-replay/artifacts.json cas12-preparation-replay/verification.json
```

The streaming reader downloads authenticated payloads into its pointer directory
using exclusive creation; run it in a fresh pointer-only directory, or fetch
only the reader/pointer there. The fetch command alone is sufficient to inspect
all original bytes. A separate exact2048-root B1 allocation is now justified by
the measured same-size input; no B2, outgoing pole, fit or sweep is automatic.
All114 holdout S matrices and equilibrium pole exclusion remain protected.
