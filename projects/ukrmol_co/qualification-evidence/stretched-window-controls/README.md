# Stretched CAS11 retained-data window and grid controls

The corrected four-rank outer pipeline completes **all 20 stages** at R=2.5
bohr, with both Pi sectors on 0.01–2.0-requested-eV grids at 0.005/0.0025-eV
spacing (399/797 points). It reuses unchanged dense inner-region amplitudes and
channels. Cost is **390.39 seconds / 0.6342 GiB**, excluding prior inner-region
preparation. Earlier serial-under-MPI and completion-wrapper failures remain
preserved alongside the successful fresh successor.

Complete grids, finite/nonnegative channel-summed cross sections, Pi symmetry,
native inputs, original exits and retained amplitude hashes reconstruct from
**310 payloads / 1947 profiles**. Common phases agree with the original
0.01–1.0-eV grid within **2.0e-7 rad modulo π**. The two extended grids agree
at shared printed phase points; complete-data fitted position/full-width
differences are at most **2.205e-6 / 4.293e-5 eV**.

The independent fit diagnostics cover three declared windows
(0.4–1.6, 0.5–1.5 and 0.6–1.4 requested eV), all four polynomial backgrounds,
three separated starting guesses and five interleaved held-point folds. The
48 complete-data fits converge without parameter-bound hits. The native fine
grid emits **MAXFIT**, so its native-fit verdict stays unset despite its printed
candidate. The complete-data independent fits retain all requested points.

**The all-background width-sensitivity gate fails:** the declared constant-to-
cubic family spans **6.72–6.74%**, exceeding the **5%** limit. Position spans
are below 1.058 meV. Terms2–4 alone span about **0.131%** in width, and held-point
residuals improve substantially as background terms are added; this is a
representation diagnostic, not a retroactive removal of the constant-background
control. Full extraction, model selection and geometry/threshold identity remain
unqualified.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-window-controls
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-window-controls/stretched-window-controls.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/stretched-window-controls
PYTHONPATH=projects/ukrmol_co/qualification-evidence/stretched-window-controls/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/stretched-window-controls/state-averaged-evidence/verify-stretched-windows.py \
  projects/ukrmol_co/qualification-evidence/stretched-window-controls/state-averaged-evidence
```

The archive includes original binary retained inputs, complete replay outputs,
both earlier failures, immutable source and independent reconstruction. Public
fetch, reconstruction and deterministic byte-identical repackaging pass.
