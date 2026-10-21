# Explicit matched-manifold property export on the fifty-component frame

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
