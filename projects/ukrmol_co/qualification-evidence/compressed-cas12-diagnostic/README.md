# Compressed CAS12 fixed-orbital rejection diagnosis

All **96 diagnostic probes** complete on the two unchanged rejected compressed
QC checkpoints: 5/8 roots, trial spaces 40/80/160, 600 cycles, both spins and
all four irreps. Runtime is **2246.99 seconds / 0.6136 GiB**. The diagnostic
execution exits zero and records the failed scientific gates explicitly.

**84 probes pass; all 12 singlet-A1 probes reject.** Other spin/symmetry sectors
converge with physical residuals below 6.911e-10 Hartree. Singlet A1 retains false
CI flags, physical residuals up to **0.19803 Hartree** and spin contamination.
Its five-root solve includes states with spin-square near **2**, although a
singlet requires zero. Stable energies across trial spaces (differences below
3.695e-13 Hartree) cannot override those failures. Changing the root request
from five to eight changes common returned energies by **1.6561 Hartree**.

The original QC records remain rejected. Passing a saved-orbital diagnostic
would not validate an orbital optimizer that used the rejected CI objective.
The [separate singlet-driver pilot](../../cas12-compressed-singlet-driver-pilot-contract.json)
tests an alternative library driver with analytic and independent Hamiltonian-
action checks, without changing these records or automatically releasing QC.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic/compressed-cas12-diagnostic.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic
PYTHONPATH=projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic/state-averaged-evidence/verify-compressed-diagnostic.py \
  projects/ukrmol_co/qualification-evidence/compressed-cas12-diagnostic/state-averaged-evidence
```

The archive independently reconstructs **344 payloads / 11201 profiles**, all
96 gate decisions, original failed QC flags, unchanged checkpoint/source hashes
and root-count/space differences. Public fetch, reconstruction and byte-identical
repackaging pass. Reconstruction verifies the rejected verdicts; it does not
turn them into converged results.
