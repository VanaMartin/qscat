# Fresh-lineage aug-TZ CI-space repairs

Both same-basis restarts of the rejected fresh CAS(10,11)/aug-cc-pVTZ
checkpoint now pass QC. The physical model, forty-component ensemble, tight
orbital/CI tolerances and default 200 CI cycles are unchanged. The numerical
CI search space is enlarged independently to **80/160 vectors**. The exact
rejected parent checkpoint and its failed fresh singlet-A1 convergence flag
remain evidence; its original exit is **one**.

| CI space | Original exit | Wall (seconds) | Kernel memory peak (bytes) | Ground energy (Hartree) | Ground z dipole (a.u.) |
|---|---:|---:|---:|---:|---:|
| 80 | 0 | 641.47 | 449634304 | −112.852303695115 | 0.222676780005 |
| 160 | 0 | 638.23 | 512544768 | −112.852303695116 | 0.222676780002 |

The serial batch takes **1280.74 seconds / 21.35 minutes**, with 32-GiB
containers pinned to CPUs 4–7. Both optimized and fresh CI flags pass, as do
the gradient, spin, Pi and MO-orthogonality gates. Relative to the rejected
parent's recorded spectrum, reoptimization changes a root by at most
**5.45e-8 Hartree**; numerical agreement did not make the original rejection
acceptable.

The two passing restarts satisfy the numerical pair gate:

- All forty roots agree within **1.667e-11 Hartree**.
- Averaged objectives differ by **−2.842e-13 Hartree**.
- Ground z dipoles differ by **−3.230e-12 a.u.**.
- Minimum core/active-subspace singular values are
  **0.9999999999999999 / 0.9999999999999992**.

This demonstrates agreement between CI-space repairs of the same fresh lineage.
Independent competing starts, fixed-orbital lowest-root coverage and UKRmol
all-root/dipole import remain required. In particular, these repairs do not
settle fresh-versus-projected larger-basis qualification.

The companion preserves two completed runs and their frozen batch source,
configurations, seed/final checkpoints, raw QC logs and resources. It also
includes the exact rejected parent, pair comparator/source/command and
acknowledged CPU-slot lease/resumption. The executable verifier reconstructs
all **80 raw final-state energies/spins**, rechecks current QC analysis and
recomputes spectral/dipole/objective pair gates from the targets. Orbital
singular values are measurements retained with the reproducible comparator.
The small read-only pair calculation briefly overlaps the subsequent TZ
import on CPUs 4–7; that overlap is explicit in its execution record.

All **129 payloads** verify by size/digest, with byte-identical repackaging
and a byte-for-byte public fetch. The [original fresh-basis companion](../fresh-basis/README.md)
retains the rejected attempt and its nonzero batch exit.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-augtz-repairs
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-augtz-repairs/fresh-augtz-repairs.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-augtz-repairs.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [continuation](../../CONTINUATION.md) records the finite 96-probe,
600-cycle coverage/residual scan queued for these two checkpoints after the
fresh CAS(10,11) TZ import. The
[input recipe](../../calibration-sa11-fresh-augtz-ci-repairs.json) names both
repairs and the unchanged rejected seed.

The later [coverage companion](../fresh-augtz-coverage/README.md) now passes
both checkpoints' 96 probes/624 eigenpair evaluations in 10.02 minutes /
0.391 GiB. The [64-root aug-TZ import](../fresh-augtz-import/README.md) now also
passes root/dipole/subspace gates in 56.01 minutes / 5.268 GiB. Competing starts
and electronic-model qualification remain open.
