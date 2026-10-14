# CO root resolution: sufficient evidence for the next controls

**The user accepts CAS12's passing 2048→4096 comparison as sufficient for now.**
Use the tested 2048-root baseline for matched qualification controls. Stop the
8192-root experiment, preserve its partial evidence and original exits, and
prioritize electronic/channel, state-continuity and extraction questions.

## What the phase numbers mean

| Model / resolution | Measured phase discrepancy, radians modulo π | Reference |
|---|---:|---|
| CAS11, 2048 | 5.1e-6 | Full 27546-root dense result |
| CAS11, 4096 | 1.0e-7 | Full 27546-root dense result |
| CAS11, 8192 | 1.0e-7 | Full 27546-root dense result |
| CAS11, 16384 | 1.0e-8 | Full 27546-root dense result |
| CAS12, 4096 versus 2048 | 1.67e-5 | Matched 2048-root result |

CAS11 has a full-spectrum numerical oracle. CAS12 currently has a matched
two-resolution comparison, so **1.67e-5 is the measured refinement difference,
not an independently measured absolute error of the 4096-root phases**. It
covers both sectors and the full 99-point energy grid. The unchanged phase
budget is 0.05 radians: the difference is about 3000 times smaller.

CAS12's maximum same-window fitted position/full-width changes across background
orders 1–4 are 0.00224/0.04195 meV. Both sectors pass raw exports, residuals,
symmetry, prefix energies, actual native inputs and channel interfaces. The
[completed pair](progress-audited-4096-20261014T0658.json) and
[local reconstruction](verification-audited-4096-20261014.json) independently
support practical root-resolution sufficiency for the next matched controls.
Absolute omitted-spectrum error, electronic-model convergence and pole identity
remain distinct, unset claims; the production sweep remains gated.

## Why stop 8192?

The [CAS11 dense-reference controls](qualification-evidence/cas11-sparse-refinement/README.md)
already passed at 2048 roots. Increasing 4096→8192 more than doubled native
runtime (3.34→7.06 hours) without reducing the reported maximum phase error.
The [16384-root replay](qualification-evidence/cas11-16k-tz-continuation/README.md)
cost 19.65 hours and confirmed an already adequate result. The 10 October
stopping rule recognized that completing a ladder is not a release requirement.

CAS12's new doubling would add empirical confirmation, but would not independently
bound its omitted spectrum. There is no observed root-resolution budget failure
to resolve. The 8192 launch justification did not demonstrate enough marginal
scientific value to warrant the cost; a separate declaration alone does not make
an experiment useful.

## Prospective stopping rule

- Treat the current root-resolution evidence as sufficient for the next controls;
  use 2048 as the economical matched baseline and retain 4096 as completed evidence.
- Launch no further root ladder by default. Require a named new discrepancy,
  model-specific coverage requirement or independent test whose outcome could
  change a scientific decision, with a justified resource budget.
- A change in geometry, basis, active space or explicit channels can require its
  own bounded resolution check; current evidence is not universal qualification.
- Preserve the original 8192 declaration and launch record. Record cancellation
  as an administrative scientific-priority decision, not a solver or budget failure.
- Complete model/channel, continuity/extraction and neutral qualification before
  the frozen restartable production pilot, sweep or held-out potential fitting.

## Recorded cancellation

The [10:34 UTC capture](progress-audited-8192-cancelled-20261014T1034.json)
records the user-approved STOP request at 10:32:24 UTC and owner completion at
10:33:08 UTC. The native stage ends after **2.50016 hours**, with Docker/native
exit 137 and owner exit 1 preserved as administrative termination. It is not a
new OOM or numerical failure. All eight B1 ranks completed Class1/3; Class4,
remaining assembly/solve, eigenpair exports, B2, outer stages and the sequence
reader did not complete. There is no reusable Hamiltonian checkpoint.

Both owner and observer services are inactive; no scientific container remains.
The observer's three immutable captures and original execution record remain
unchanged, with its administrative stop documented separately. Frozen 47-source /
665-input hashes and 242 installed engine artifacts reverify. The final capture
pins 739 immutable remote records, including itself.
