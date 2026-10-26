# A2-only driver pilot — stops at fifth molecular root

**The bounded pilot stops inconclusive; no selector or fresh QC is released.**
Exact non-A1 A2 controls pass both installed drivers against analytic singlet
energies (maximum error1.23e-15Ha), and the small H4 mixed-spin A2 optimization
passes its independent final-integral energy check. Its core overlap changes
to0.9996117274, so the optimizer control is not a fixed-orbital tautology.

The one five-root CO request returns flags `true,true,true,true,false`.
Roots1–4 meet the frozen state gates; root5 fails:

| Fifth-root diagnostic | Measured | Gate |
|---|---:|---:|
| Full physical Hamiltonian residual /Ha |8.79527e-8|1e-9|
| Penalized residual /Ha |1.86017e-7|1e-9|
| Penalized/full action difference /Ha |1.50439e-7|1e-9|
| Restricted/full action difference /Ha |1.30703e-7|1e-9|
| S² action norm |5.59441e-8|1e-9|
| Alpha/beta exchange error |8.38034e-9|1e-9|
| S² expectation |3.94256e-16|absolute error1e-6|

A tiny spin expectation does not override failed CI/action gates. Zero spatial
leakage and CI/MO orthogonality also do not qualify root5. The eight-root request
is **not run**, so requested-count agreement remains unset. There is no larger
space/cycle/alternate-seed repair or retry. The original failed TZ optimization
remains rejected; this is neither a physical TZ rejection nor a basis comparison.

[Independent public replay](../../verification-singlet-a2-pilot-public-20261026.json) ·
[Frozen contract](../../singlet-a2-pilot-contract-20261026.json) ·
[Physics and release boundary](../../../../docs/physics/co-singlet-a2-driver-control.md)

Source/tests/exact owner `dc43f53` freeze/push precedes14:03:12UTC launch;
finish14:03:38UTC26October. Pinned image/PySCF2.11.0, original failed checkpoint,
unchanged160-space/200-cycle/energy/residual/lindep/penalty settings. Four physical
cores0–3 /8GiB /no network. One CO CI request; zero CO QC/native/maps. The declared
synthetic and H4 validations are separate from the zero-CO-QC counter.

Original profiled-command/container/owner/service exits are all**2**; scientific
exit2 and verdict unset. No timeout/OOM. Cost26.085183s profiled /26.423980s owner,
63.130420CPU seconds,585584640bytes (0.545368GiB) kernel peak,131 resource samples.
Sampled run-directory allocated peak33214464bytes excludes source/reference data.
Other host workloads, including ComfyUI, are untouched.

Public replay verifies seven payloads/fifty archived files/twenty frozen sources,
seven failed-run input/four source prerequisites; all twenty sources also match
frozen Git blobs. It independently checks analytic A2 energies, H4 J±K energies,
orbital movement, saved CI/action norms, symmetry/orthogonality, raw log, absence
of eight-root files, telemetry and exits. **CO full-Hamiltonian/S² actions are
saved external contractions, not independently recomputed actions.** No producer
code/PySCF import, eigensolve or optimization occurs during reconstruction.
An initial public-reader HTTP403/exit1 is preserved; only the request User-Agent
is corrected to the existing artifact-fetch convention, public replay exit0.
No scientific bytes/exits are repaired.

The hypothetical1/40-weighted replacement of the five rejected printed A2
energies gives−0.07751764535Ha. This is **unqualified diagnostic substitution**,
including the failed fifth root and the prior rejected sector, not a newly
optimized objective, target, energy improvement or accuracy result. The pilot's
pass-only substitution stage was correctly not reached.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/singlet-a2-pilot
uv run --no-sync python -B projects/ukrmol_co/qualification-evidence/singlet-a2-pilot/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/singlet-a2-pilot/artifacts.json \
  verification-singlet-a2-pilot-fresh.json
```

**Disposition: stop this driver allocation.** Exact controls establish limited
non-A1 feasibility, not a validated molecular route. No selector change, fresh
forty-state CO optimization, coverage/import/properties/native/scattering/B2/fit
or automatic successor. Retain B1 model-position rejection, nine working reference,
2048 roots,114 untouched holdout S matrices, equilibrium exclusion and parked
new50. A separate saved-evidence production-route decision, not another solver
budget increase, is the only bounded next deliverable; no calculation is queued.
