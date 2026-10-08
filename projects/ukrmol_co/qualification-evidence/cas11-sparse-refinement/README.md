# CAS11 native spectral refinement and electronic successors

This companion independently reconstructs **1506 payloads / 477179 resource
profiles** from eight completed components. It preserves both the original
16384-root OOM and the interrupted local watcher's QC-reporting key error.

## Native selected-spectrum controls

Both Pi sectors pass eigenpair, boundary-export, full-grid phase and unchanged
1.8–3.3-requested-eV fixed-window position/width checks at **2048, 4096 and
8192 roots**, against the same 27546-dimensional CAS11 dense reference.

| Roots | Maximum phase error / rad modulo π | Wall / h | Kernel peak / GiB |
|---:|---:|---:|---:|
| 2048 | 5.1e-6 | 2.454 | 7.693 |
| 4096 | 1.0e-7 | 3.335 | 7.365 |
| 8192 | 1.0e-7 | 7.060 | 13.16 |

The original **16384-root B1 SCATCI** exits **137**, reaching exactly its
**24-GiB cgroup limit** after **3309.93 seconds**. Docker records the OOM event;
the original controller exits one and launches no later stage. This is a
preserved resource failure. The separately declared
[40-GiB successor](../../sparse-scattering-cas11-memory-repair-contract.json)
retains the scientific inputs, solver and gates.

The other components are the **96 passing fifty-component equilibrium coverage
probes**, completed stretched CAS11 dense-anchor pipeline, and passing stretched
CAS12 TZ/aug-TZ QC. QC success does not establish independent larger-basis root
coverage or native import. The stretched anchor's fit verdict remains unset.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement
tar -xzf projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement/cas11-sparse-refinement.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement
PYTHONPATH=projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement/state-averaged-evidence/verify-cas11-refinement-bundle.py \
  projects/ukrmol_co/qualification-evidence/cas11-sparse-refinement/state-averaged-evidence
```

The archive contains raw CI/boundary records and the dense oracle, solver
telemetry, phase grids, native outputs, electronic checkpoints, resource records,
historical owner/source snapshots and portable independent reconstruction.
Its deterministic repackaging is byte-identical. Full-model convergence,
resonance extraction/continuity and CAS12 scattering remain separate gates.
