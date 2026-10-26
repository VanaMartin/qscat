# Compressed CAS12/TZ QC prerequisite — stopped at singlet A2

**The QC prerequisite rejects; the DZ/TZ physical comparison is inconclusive.**
Orbital optimization converges in25 macroiterations, final gradient
2.30964e-8 below1e-7, but **singlet-A2 CI is unconverged**. The raw final five
nominal singlet-A2 states print S² values7.78e-5,3.87e-5,.0017931,2.0000126,
1.9952972 instead of zero. The last two are near triplet spin, not valid singlet
ensemble roots. These printed values independently confirm an additional spin
problem; the original code stops on the CI flag **before** its fresh-CI/spin audit.

This is not a failed physical basis-insensitivity test, a timeout, an OOM or a
command-line parsing error. Systemd's `INVALIDARGUMENT` label is its generic
name for exit2; the saved traceback identifies the actual target CI failure.

[Result/resources](../../progress-compressed-cas12-tz-qc-result-20261026.json) ·
[Independent public reconstruction](../../verification-compressed-cas12-tz-qc-public-20261026.json) ·
[Frozen contract](../../compressed-cas12-tz-qc-contract-20261026.json) ·
[Allocation rationale](../../../../docs/physics/co-compressed-cas12-basis-control.md)

Source/tests/exact owner `317ae15` precede26October11:17:37UTC launch;
finish12:39:16UTC. All eight pinned-runtime controls pass with no skips before
the new CO input. Exactly one QC attempt uses the qualified DZ checkpoint,
cc-pVTZ at1.9bohr, active6/3/3/0, unchanged40-component averaging, ordinary200
CI cycles/space160 and150 macroiteration cap. No numerical source, optimizer,
trial space or ensemble was changed after launch.

Original **runner exit1**, **analyzer not run**; scientific/command/container/
owner/service exits**2**. Source integrity and all nine input pins remain
unchanged. Runtime is4898.179626s profiled /4898.480259s scientific /
4899.111255s owner, about1.36h, within forecast/caps. CPU14569.732293s;
kernel peak1253928960bytes (1.167812GiB);24423 telemetry samples reconstruct. The measured
run-directory allocated peak1601536bytes excludes copied cache/source; scratch0.
Systemd23.8M/1.483CPU seconds are supervisor metrics, not the scientific container.
Allocation: four physical cores0–3,24GiB, no network. No timeout/OOM/retry.

Public replay authenticates970 archived files/23 frozen sources/nine root/
seven source prerequisites. The sources also match frozen Git blobs. It
independently parses original diagnostics, the forty raw printed state energies,
weights/spins, optimizer trace, failure traceback, original exits and resource
samples. It checks the1/40-weighted raw energies against the rejected objective
and the printed average S²=1.099930485 against the intended mean1.0. This is
**saved-output/gate reconstruction, not a new CI solve or independent wavefunction
Hamiltonian action**. No QC retry or scientific output repair is performed.

The checkpoint is preserved with SHA256
`2bfbcafabac0fa07ab2af5e31888eabddaa7dabe465a50cde07e386e33f52cf6`.
No successful `target.json`, Molden export or analyzer result exists. A converged
orbital gradient does not qualify an objective containing wrong-spin,
unconverged roots; a later fixed-orbital CI result alone could not retroactively
qualify this optimization. Do not use these energies/checkpoint as a physical
TZ target, basis comparison, property input or production reference.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-cas12-tz-qc
uv run --no-sync python -B projects/ukrmol_co/qualification-evidence/compressed-cas12-tz-qc/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/compressed-cas12-tz-qc/artifacts.json \
  verification-compressed-cas12-tz-qc-fresh.json
```

**Disposition: stop the current TZ route.** The conditional24-probe proposal
is blocked, not released on this checkpoint. No repair/retry, alternative seed,
larger trial space, native import, DENPROP, scattering, pole, B2 or other model
job follows. The prior86.77meV B1 model-position rejection remains; neither model
is qualified. Next is a separately reviewed singlet-A2 correctness/optimizer-
objective decision from saved evidence, not another full QC allocation. All114
holdout S matrices/equilibrium exclusion, nine working reference,2048 baseline
and parked new50 remain protected.
