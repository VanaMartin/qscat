# Explicit forty-state property export on the fifty-component frame

**Later physical comparison,16:49UTC:** the property export remains qualified,
but the [cross-frame retained-subspace screen](ENSEMBLE_FRAME_COMPARISON.md)
rejects both A1 sectors. Selecting five roots per sector does not preserve the
same physical retained manifold across orbital averages. No new properties or
Hamiltonian are generated to repair this rejection.

## Completed and independently read — 21 October,15:38 UTC

**The selected40-state property export passes**, native command/container/owner
**0/0/0**, finishing **15:29:41UTC**. Independent saved-byte reading reconstructs
all40 IDs, five states in every spin/irrep, **676 finite unique dipole/quadrupole
records**, their symmetry indices and the complete sector-to-property map.
Maximum independent-CI energy discrepancy **4.91610e-10 Hartree** passes1e-7;
ground-dipole difference **−1.12928e-11a.u.** passes1e-5. Previously missing
singlet-A2 roots4/5 now have property IDs **39/40**. This establishes the explicit
selection, not physical cross-frame identity or accuracy of every transition moment.

The original and new property files share **38 same-frame states /620 moment
keys**, and all620 values are **exactly equal**. This is a read-only subset
observation, not a retrospective acceptance gate or comparison between different
orbital frames. The raw40-state file SHA256 is
`16ec0905d1487dddb7ee7f7f891e62b3e8157a44da21d58c74d698f99ede9ece`.
Its input deck differs only in the declared `ntgtl` selection; all other native
bytes/tokens stay unchanged. The original six-root seed rejection/exits1/1/1
remain separate.

Measured cost **2232.151s /37.203min**, sampled CPU2231.769s /37.196CPUmin,
kernel high-water **469,364,736 bytes /.43713GiB**,11,150 samples. Native DENPROP
itself takes2231.204s. The run exceeds the30min expectation but stays inside45min
upper /60min hard bounds. Forty selected states did not make this measured run
cheaper than the earlier35.5min global-low50 run; these are not isolated timing
benchmarks. Four reserved cores0–3 /one native thread /4GiB remain unchanged.
Zero orbital optimizations, target diagonalizations or scattering inner solves.

All12 producer sources /13 original pins /10 copied native units verify before/
after. Source freeze **`96b8ed9`** precedes production; the separate independent
reader/tests freeze **`5cbefa0`** before saved reconstruction. Its seven positive/
negative cases plus four export-deck cases pass locally. The
[public evidence](qualification-evidence/fifty-properties/README.md) fresh-verifies
all **8 payloads /22 original generated/owner files**, with9 prerequisite payloads
reauthenticated. Ten unchanged copied input files are authenticated but not
republished; their exact bytes remain in the prerequisite target archive.
The new collection preserves all other generated files, including native
intermediates,93,323,078 uncompressed bytes, SHA256
`d35630ff8ae9ac6ccc7b451f514bdf0a9b79e71ec6575dc767c7856f1522e3d6`.
[Result/resources](progress-fifty-properties-result-20261021.json) and
[independent reconstruction](verification-fifty-properties-public-20261021.json)
retain the full maps, moment counts, hashes and original exits.

No scientific job remains after completion. **Next:** physical cross-frame
channel/subspace correspondence and compatible scattering preparation/cost on
the qualified orbital frame, with no repeated QC optimization. Property IDs and
equal target counts alone do not prove that two orbital averages retain the same
physical state manifold. No new Hamiltonian, phase or width claim is released.

## Earlier launch state

**Status21 October,14:54UTC:** one selected40-state DENPROP calculation is
running, source freeze **`96b8ed9`**, owner PID278257. It starts14:52:28UTC,
expected30min /upper45min /whole-container cap60min. Physical cores0–3 /
one native thread /4GiB; all12 frozen sources and13 original input pins verify.
Native input verification completes; final scientific verdict/exits remain
pending. Persistent outputs:
`prepared/fifty-selected-properties-source-20261021/outputs/properties/` on
Sadaharu. [Launch capture](progress-fifty-property-launch-20261021.json) records
the actual command, hashes and original pending statuses. The completion watcher
is read-only; no other numerical work is submitted by it.

The completed new orbital frame now independently passes48 CI probes and all74
saved native target energy/vector/header/EOF checks. The original seed-frame
preservation rejection is unchanged. The immediate missing input for an ensemble
comparison is the **property selection**, not another target diagonalization.

The previous global-low50 property file contains only **three singlet-A2** roots.
The existing five-A2 reference requires five roots in each spin/irrep. Keeping
the global-low50 export would compare different channel manifolds, even if the
total state count looked adequate.

## One bounded native property calculation

The [contract](fifty-property-export-contract-20261021.json) and
`fifty_property_export.py` require the authenticated independently qualified
frame (`55e89a6`), all ten native input units, the original deck and the completed
qualification result/exits. Only one native card changes:

```text
ntgtl = 7,6,6,8,6,8,3,6,   →   ntgtl = 5,5,5,5,5,5,5,5,
```

All other cards, including `ipol=2` and nuclear contributions, stay unchanged.
Inputs are copied into a new directory from read-only producer mounts. One
DENPROP invocation computes the selected40-state property representation. No
orbitals are optimized, no target roots are rediagonalized and no scattering
Hamiltonian is assembled.

The existing property-inventory reader reconstructs all40 native state IDs and
their explicit spin/irrep/root map, using independent fixed-frame CI energies
with the unchanged1e-7-Hartree budget. It rejects incomplete/duplicate IDs,
wrong sector counts, unordered/nonfinite energies, duplicate/nonfinite moments,
uncovered states or missing dipole/quadrupole ranks. Ground dipole must match
the qualified QC frame within1e-5a.u. All copied native inputs and13 original
pins are authenticated before/after. Inventory coverage is not independent
accuracy validation of every transition moment.

## Cost and release boundary

The previous global-low50 DENPROP stage measured **2131.91s wall /2136.58 sampled
CPU seconds**, .65450GiB sampled peak memory.current /.36610GiB anonymous memory.
It is effectively serial; reserving eight cores does not make it eight-way work.
Forecast **1800s expected /2700s upper /3600s whole-container cap**, with3300s
native timeout, one BLAS thread on reserved cores0–3 and4GiB memory. Preserve all
native/container/owner exits and stage logs. No retries or tolerance changes.

A pass releases **physical cross-frame channel/subspace correspondence and
sparse-preparation/cost design only**. The target-only CIDATA/integral layout
has not been qualified as a direct replacement for full scattering preparation;
the latter must preserve the qualified orbital frame without another QC
optimization and establish compatible continuum integrals/native inputs.
No new scattering solve, phase comparison, width qualification or production
sweep is claimed by this export. The independent model/extraction/continuity/
neutral-reference and restartable-pilot requirements remain.
