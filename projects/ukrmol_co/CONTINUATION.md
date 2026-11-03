# CO: current handoff

## Oxygen-affinity diagnosis requested

The user now selects a separate reference-correctness investigation: audit the
O/O- affinity and identify any incorrectly configured calculation. This work
is isolated on its own branch, not added to the PR under review. The saved
state audit finds consistent occupations/core/subtraction, but ROHF-to-UHF
conversion does not semicanonicalize the spin Fock blocks. The
[bounded orbital-convention control](../../docs/physics/oxygen-affinity-orbital-convention.md)
tests that specific hypothesis using the saved5Z checkpoints/amplitudes, with
zero SCF/CCSD iterations. Four physical cores/16GiB,1800s science cap,600s active
checks and immediate terminal observation. Next reconstruct the actual result;
no automatic6Z,reference reoptimization,molecular/scattering or holdout successor.

## Goal and disposition

A physically grounded neutral and lowest-doublet-Pi resonance/bound-anion curve,
then a QSCAT two-dimensional potential with excluded-geometry validation.
**No production model, fit or sweep is released.** Low-Pi channel-specific
dynamics is distinct from total thermal attachment near 10 eV.

**Current decision:** stop endpoint basis-insensitivity/production inference and
redirect to neutral dissociation, correlation and reference treatment. The
[settled-results note](../../docs/physics/co-results-and-rejected-routes.md)
records the evidence and rejected routes without the launch chronology.

## Latest completed evidence

- [Matched aug-QZ/5Z atoms](qualification-evidence/atomic-5z-endpoints/README.md):
  all original calculation exits zero; 5Z owner 114.673835 s, kernel high-water
  11.883816 GiB on four physical cores, below the unchanged 16-GiB allocation.
  O/C EA changes are +23.533750/+5.950705 meV. O rejects the 10-meV diagnostic;
  all three fixed-reference offset changes reject 20 meV. No De/D0/CBS claim.
- [Original 5Z startup stop](qualification-evidence/atomic-5z-launch-stop/README.md):
  controls exit 1, outer exits 2, zero integral/SCF/CC calls. A separate exact-source-
  byte control fixes Python-version-dependent AST metadata hashing; no physical
  or numerical gate changed and no successful atom was repeated.
- [Stretched DZ/TZ poles](qualification-evidence/stretched-pole-counterpart/README.md):
  both finite numerical results qualify; position/full-width insensitivity rejects.
- [Aug-TZ preparation](qualification-evidence/stretched-augtz-preparation/README.md)
  qualifies; [unchanged export](qualification-evidence/stretched-augtz-inner-preflight/README.md)
  is blocked before an inner allocation by threshold-crossing grid rows.

No campaign quantum job or monitor is active. Both atomic jobs completed before
their first 600-s periodic check; immediate terminal observations are retained.

## One bounded next deliverable

At most **60 minutes of saved evidence/source review, zero new quantum calls**:
select one stable spin-consistent neutral dissociation/correlation treatment or
an explicitly declared hybrid/literature reference. Require overlap with the
accepted well data, a route to fragment consistency, and one costed small
diagnostic set. If the evidence cannot select a treatment, record the blocker.
Reconcile spin-free/lowest-J conventions before precise external EA scoring.

No 6Z, automatic extrapolation, dense neutral grid, augmented inner/pole work,
fit, sweep or holdout action follows automatically. A new calculation needs
its own tested, frozen, pushed source/input/owner release. Historical contract
flags do not authorize a rerun against current source.

## Standing protections

Preserve accepted calculations, immutable null input ledgers, original failures,
exit statuses and measured costs. Keep all 114 S holdouts untouched, the specified
equilibrium outgoing-pole exclusion, nine as working reference, new50 parked and
the 2048-root baseline. Budgets remain in [production decisions](PRODUCTION_DECISIONS.md).
Review goal convergence every one or two active working days, on resumption,
before a new branch and after a decision-changing result.

Automatic continuation applies within approved bounded stages, with scientific,
resource and scope stops. Active runs target 600-s checks plus immediate terminal
observation; preserve actual gaps against the requested 5–15-minute interval.
There is no perpetual scheduler while idle. Prefer Sadaharu: four physical cores
for small jobs, eight for large and twelve for very large when justified;
paid provisioning remains deferred.

## Evidence, not a live launch ledger

[Run/replay instructions](RUNNING.md), [evidence index](qualification-evidence/README.md)
and [research queue](RESEARCH_PLAN.md) are the restart points. Remote originals
are under `/home/kooza/ukrmol/co-sa-20261006` on Sadaharu; public authenticated
bundles, not transient PIDs or local scratch paths, are the portable record.
Superseded handoffs and figures are in the
[campaign-history archive](qualification-evidence/campaign-history/README.md).

Main retrieval still reports incompatible shared-service dependencies. Last-known
indexed main is `423d5770d3de337571f51e62702ae919bad97e8c`; integrity/current remote
currency cannot be verified through that service. Current source, scoped Git
differences and executed-source evidence govern. No service intervention was made.
