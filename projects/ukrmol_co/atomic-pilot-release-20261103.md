# First four-atom aug-QZ release and monitoring

## Physical decision

User continuation of the reference-first plan selects the smallest missing
endpoint calculation: C(3P), O(3P), O-(2P), C-(4S) in aug-cc-pVQZ, spherical,
nonrelativistic spin-free, one frozen1s orbital per atom. Existing neutral
QZ/5Z data and all scattering data are reused without modification or rerun.
The four finite-method totals and their sums/EAs will determine whether the
reference prescription is usable before choosing a molecular dissociation model.
They cannot by themselves qualify molecular De/D0 or absolute atomic accuracy.

[Exact release](atomic-asymptote-release-20261103.json),
[immutable input ledger](atomic-asymptote-inputs-20261103.json) and
[physics protocol](../../docs/physics/co-atomic-energy-references.md) define the
scope. `atomic_references.calculate` constrains ROHF occupations, checks spin and
majority s/p orbital character, supports internal stability, and converts to UHF
representation without further SCF relaxation before UCCSD(T). Unavailable external
ROHF stability stays null; full correlated L2 purity is not asserted.

Atomic driver, durable owner, runtime controls and release freeze/push before any
quantum call. Six no-quantum controls check state/core/charge/spin inventory,
known-truth energy-zero identities, refusal of unauthorized/expanded releases,
first-timeout exit preservation/no peer, and exact owner allocation. Runtime
source support was checked previously in the same pinned PySCF2.11 image.

## Bounds and autonomous continuation

Four sequential atoms, four distinct physical cores0–3/16GiB,12GB PySCF memory,
minimum24GiB available RAM/10GiB disk. Per-atom science cap1800s,total7200s;
container7500s,durable service7800s. No measured runtime is claimed in advance.
The2hour narrow implementation bound remains separate from compute cost.

The user explicitly requests status checks every5–15minutes and automatic
continuation unless blocked. During this active allocation the owner writes
status every **600seconds**, plus an immediate terminal check. Resource sampling
continues every0.2s through the reused profiler. Actual scheduling gaps must be
reported rather than rounded into compliance. Ten-minute checks apply to active
runs; this creates no perpetual scheduler or unsupported future-session reminders.

Successful atoms advance automatically within these four cases. The first
scientific/input/stability/timeout/OOM/resource failure stops without orbital
retry, peer, method/basis fallback or suppressed original exit. After completion,
automatically read and qualify the saved evidence; do not treat execution0 as an
absolute-EA or molecular-dissociation pass. A later quantum stage needs its own
goal decision, bounded tested contract and freeze, not another permission question
when it remains within approved scope. A scientific blocker or out-of-scope
decision must be reported instead of bypassed.

Per-atom sampled memory values are observations; kernel `memory.peak` is the
cumulative container high-water mark, not an independent peak for each atom.
All raw logs/checkpoints/diagnostics, original command/container/owner/service
exits and source hashes are retained. No image rebuild, paid resources, augmented
inner solve, new neutral grid, fit/sweep or114-holdout analysis is released.
