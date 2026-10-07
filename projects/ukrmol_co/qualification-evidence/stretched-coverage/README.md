# Stretched-DZ CI-repair coverage

This companion preserves the completed **96 fixed-orbital probes** of the
CAS(10,11), R=2.5-bohr, cc-pVDZ 80/160-vector QC repairs. It completes after
the main 46-attempt supplement is frozen. The original controller exits **1**:
both targets fail the existing strict requirement that every ensemble root
converge in every requested probe.

Of 96 probes, **87 fully converge**. All **64 probes at spaces 80/160** fully
converge. Three five-root space-40 probes fail the fifth ensemble-root flag:
triplet A1 on both checkpoints and singlet A2 on the space-160 checkpoint.
Six further space-40 controls fail only extra-root flags. Every probe's first
five energies agree with the recorded target within 9.95e-14 Hartree, including
the failed-flag probes, and all spins agree. Energy agreement does not override
a failed convergence flag or the recorded gate.

The repaired checkpoints agree in roots within 2.63e-12 Hartree, objective
within 7.11e-14 Hartree, dipole within 4.81e-11 a.u. and minimum active-subspace
overlap 0.9999999999999984. Numerical pair agreement passes, while coverage
qualification remains rejected. The worker takes **559.24 seconds / 0.340 GiB**
and resumes the waiting CPU-slot owner.

The later [iteration refinement](../stretched-iteration/README.md) reproduces
all nine failed settings at 200 cycles and repairs them at 600, then passes
complete 48-probe scans of both checkpoints with measured Hamiltonian residuals.
That supports coverage under the refined iteration contract; this original
200-cycle rejection remains preserved.

The same worker measures the matched CAS(10,10) initial/final active-subspace
overlap: minimum **0.9999999998858343**, with core minimum 0.9999999999996424.
These support continuity of the matched control used by the main supplement.
They do not certify active-space, basis or continuum convergence.

The archive preserves raw CASCI spectra/spins, convergence flags, target inputs,
checkpoint copies, sources, image/command/lease records and the original failure.
All 199 payload sizes/digests and raw spectra were checked; repackaging is
byte-identical. The manifest distinguishes the numerical source snapshot from
the committed analysis tools. Together with the main archive, 256 CI probes are
publicly retained; its 46 run and eighteen batch records remain the same.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-coverage
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-coverage/stretched-repair-coverage.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-stretched-repairs.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

See [the main supplement](../README.md) and
[the live continuation](../../CONTINUATION.md) for the pending imports,
basis/start repairs and continuum controls.
