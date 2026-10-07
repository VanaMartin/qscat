# Fresh TZ/aug-TZ state-averaged targets

This companion preserves the completed four-entry fresh-RHF-start QC batch at
R=2.1323 bohr: three passing SA-CASSCF targets and one rejected fresh CI audit.
All use ten active electrons, two inactive core orbitals, five roots per
singlet/triplet C2v sector and forty equal-weight ensemble components. The
orbital/CI tolerances and default forty-vector CI trial space match the tight
DZ campaign. These are QC-only records; they have no UKRmol import or scattering.

| Active space | Basis | Original exit | Wall (minutes) | Ground energy (Hartree) | Ground z dipole (a.u.) |
|---|---|---:|---:|---:|---:|
| CAS(10,11) | cc-pVTZ | 0 | 39.75 | −112.901856489893 | 0.143618428763 |
| CAS(10,11) | aug-cc-pVTZ | 1 | 22.80 | Not exported | Not exported |
| CAS(10,12) | cc-pVTZ | 0 | 112.87 | −112.916254208686 | 0.296445434784 |
| CAS(10,12) | aug-cc-pVTZ | 0 | 91.74 | −112.898624623768 | 0.254970756656 |

The serial batch takes **16030.77 seconds / 4.453 hours** on CPUs 4–7, with
32-GiB container limits. Passing-run kernel peaks are 346562560/799301632/
875302912 bytes, below 0.816 GiB each. The failed entry has orbital and optimized
CI convergence, but the fresh singlet-A1 CI flag is false. Its fresh energies
agree within 1.42e-13 Hartree; agreement does not override the failed flag.
The original batch exit remains **one**, and no target is exported for that entry.

These fresh starts expose appreciable electronic-model sensitivity:

| Change | Ground energy change (eV) | Ensemble objective change (eV) | Ground z dipole change (a.u.) |
|---|---:|---:|---:|
| CAS(10,11)→CAS(10,12), cc-pVTZ | −0.391782 | −0.426914 | +0.152827 |
| cc-pVTZ→aug-cc-pVTZ, CAS(10,12) | +0.479725 | −1.795712 | −0.041475 |

The variational quantity being minimized is the **forty-component ensemble
average**, so a lower objective can accompany a higher ground-state energy.
Sorted same-sector spectra change by as much as 3.44/4.99 eV in these two
comparisons; ordinal labels do not establish state identity across models.
Minimum initial-to-final active-subspace overlaps are 0.171594/0.152714/0.036790
for the three passing targets. These are measurements of orbital reorganization,
not evidence that fresh and projected starts share an optimum. Independent
fixed-orbital coverage, competing starts and all-root/dipole imports remain gates.

The public archive preserves the frozen batch source and commands, four run
configurations and resource records, final checkpoints, all convergence
diagnostics, raw QC output and the executable verifier. It reconstructs all
**160 raw final-state energies/spins**, reruns current QC analysis for the three
successes, and checks the rejection without modifying its exit or exporting a
target. All **194 payloads** verify by size/digest, with byte-identical
repackaging and a byte-for-byte public fetch. The main supplement retains its
original aggregate and archive.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-basis
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-basis/fresh-basis-qc.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-basis.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [continuation](../../CONTINUATION.md) records the finite coverage/residual
scans and the [space-80/160 fresh-lineage repairs](../../calibration-sa11-fresh-augtz-ci-repairs.json).

The later [coverage/residual companion](../fresh-coverage/README.md) now passes
all three targets' 144 probes/936 eigenpair evaluations under its 600-cycle
contract. Competing starts and independent imports remain pending.

The rejected aug-TZ CAS(10,11) checkpoint subsequently receives two passing
[space-80/160 fresh-lineage repairs](../fresh-augtz-repairs/README.md). Their
roots/dipoles and core/active subspaces pass the numerical pair gate; independent
coverage, competing starts and import remain pending. The original failure
and this companion's historical aggregate are preserved.
