# Orbital convention does not explain the oxygen-affinity deficit

Saved5Z O/O− replay completed3November2026 at14:13:37UTC. Semicanonical
EA **1.426635554eV**, original **1.426606459eV**: **+0.029095meV**, below
the predeclared0.1meV meaningful-effect criterion. Density/core and CCSD-energy
invariants pass;occupied/virtual spin-Fock residuals fall below1.94e-13Ha.
Original endpoints are not replaced. Absolute accuracy is not established.

Independent reconstruction authenticates24 raw files,14 executed Git source
files from`fb063d9065df7f329a801cd9b4cf4150b6303d93`,the original5Z six input
pins,native triples logs,energy subtraction and saved rotation/density witnesses.
Full double-amplitude transformations are checked by the executed CCSD-energy
invariant,not independently re-evaluated by the saved-array reader.

All original control/worker/science/container/owner/launcher exits0,no OOM/timeout.
Zero SCF/CCSD iterations. Four physical cores0–3,16GiB,owner76.515963s,
worker CPU138.717374s,cumulative kernel peak12,688,596,992bytes=11.817177GiB.
Immediate terminal observation occurred before the first600s periodic check.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/oxygen-affinity-orbitals
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints
uv run --no-sync python -m projects.ukrmol_co.verify_oxygen_affinity \
  projects/ukrmol_co/qualification-evidence/oxygen-affinity-orbitals/evidence.tar.gz \
  projects/ukrmol_co/qualification-evidence/oxygen-affinity-orbitals/collection.json \
  projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/evidence.tar.gz
```

See [the diagnosis](../../../../docs/physics/oxygen-affinity-diagnosis.md).
No molecular/scattering/fit/holdout successor is qualified.
