# Native continuum-target products: ordering is not a spin convention

This is a source-resolved admission contract, **not a qualified native adapter**.
The [four-block metric](co-full-interior-metric.md) and
[core-sector contractions](co-core-sector-contractions.md) require a coherent
creation-on-target representation in a fixed determinant convention. A passing
one-body self metric or generic spin-lowering test does not determine that map.

The [source audit and hash-bound locators](https://data.qscat.org/co-native-product-phase-audit-20261120/assessment.6d59015f7a73.md)
identify native source `d88fbc9486dccb4d2a433e78b78cf87a656ceb79` and the
inspected subset. Raw native source is not part of the public evidence.

## Three different indices

Keep these indices distinct:

- A **contracted row** is ordered by native target symmetry, target state within
  that symmetry, then continuum spatial orbital. If symmetry t retains n_t
  states and c_t continuum orbitals, its one-based row is
  `sum_(v<t) n_v c_v + (I-1)c_t + j`.
- A **prototype CSF** belongs to a target configuration, not a target eigenstate.
  The inspected Hamiltonian uses pairs of representative continuum-orbital
  prototypes and a stride of two between target configurations. Those pairs
  support same-orbital versus different-orbital matrix-element expansions.
- A **spin component** is an actual continuum alpha/beta creation operator and
  the corresponding target spin projection. It must come from the determinant
  expansion or a validated coupling map, not a prototype's odd/even position.

In particular, interpreting the two adjacent prototypes as alpha and beta
components is unsupported and conflicts with their use as two continuum-orbital
representatives. Spin coupling already resides in each CSF's coefficients.
Do not multiply each prototype by a new Clebsch–Gordan factor on that assumption.

## Coefficient axes and phases

The native target object stores `eigenvectors(configuration, state)`. Its
contracted-Hamiltonian helper accepts the configuration first and state second,
despite ambiguous local parameter documentation. The active call sites, loader
unpacking and array declaration agree on the actual axes. Never transpose target
coefficients based on the helper's prose alone.

The target loader calls `CIRMAT` with set/state selections, target-CI sequence
and `phase_index`, then unpacks the returned vectors into that object. Its
`Phase` pointer is documented as vestigial; that is **not** proof that the
`CIRMAT` phase argument is inactive. Raw CIDATA coefficients cannot be declared
equivalent to loaded coefficients without resolving that boundary.

The subsequently resolved `CIRMAT`/`CPAKPZ` implementation makes the boundary
explicit. With continuum contraction enabled, a selected raw target coefficient
is multiplied by the product of the corresponding target-file phase label and
scattering-CONGEN target phase label. Both are signs. The row transformation is
therefore `C_loaded[m,I] = p_scattering[m] p_target[m] C_raw[m,I]`, with target
selection/order unchanged. It must be applied once; multiplying only one phase
table or independently rephasing the resulting target changes the product basis.
This source rule is not proof that a particular decoded pair uses matching labels.

For the electron-only D2h branch, `MKORBS` gives adjacent spin addresses for
each spatial orbital, spin0 followed by spin1, without the axial-symmetry
fourfold pattern. This admits the interleaved-address conversion only after the
saved header establishes that symmetry branch and electron-only case.

Packed CSFs reconstruct each determinant by substitutions into its reference
occupation. The native manager folds sorting parity into each determinant
coefficient when converting to its ordered bit representation. Conversion to
the project's alpha-then-beta convention must include the complete permutation
once, not repeat a native sort sign or discard the reference-slot order.
The orbital table also supplies spin/spatial labels; numeric address parity
alone is not a source-verified map for an arbitrary primary-model export.

## Smallest discriminating mapping witness

Use one already accepted primary-CAS11 frame and source-matched target/prototype
CSFs under a separate read-only admission. No eigenvector recovery or native
execution is needed merely to test the map. The CAS12-specific target decoder
cannot be silently reused for the thirteen-spatial-orbital primary partition.

1. Pin target selections, raw/loaded coefficient phases and the explicit native
   orbital address/spin map. Resolve `CIRMAT` or export its already loaded target
   coefficients only under a separately authorized saved-reader scope.
2. Independently decode representative singlet and triplet target CSFs and both
   continuum prototypes into canonical sparse determinants. Use enough target
   configurations to expose configuration-dependent phases; a single leading
   determinant or a squared overlap is not sufficient. Bound selected records
   before decoding instead of expanding the entire scattering space.
3. For each prototype, verify exactly one continuum electron in every term and
   remove its explicit spin orbital with signed exterior deletion. Compare the
   residual target functions with the corresponding target spin projections,
   including normalized spin lowering where needed.
4. Reconstruct the full **signed** coupled CSF, term by term, from creation-on-
   target components. Check both representative continuum orbitals, not just
   their norms or a diagonal self Gram. Establish any common/configuration phase
   relation to the coefficient basis actually loaded by the native contraction.
5. Contract the verified map with unchanged target/state coefficients and test
   selected native row orderings. Do not rephase a target while leaving its
   continuum coefficients or L2 coupling unchanged.

As an algebraic diagnostic only, for a doublet M=+1/2 coupled from triplet
targets, write

    g_alpha a†_(c,alpha) T_(1,0) + g_beta a†_(c,beta) T_(1,1).

With unit target projections related by the conventional spin-lowering operator,
total spin raising annihilates the doublet only if
`sqrt(2) g_alpha + g_beta = 0`; unit normalization gives squared magnitudes
1/3 and 2/3. A singlet target instead contributes its M=0 state with continuum
alpha. These relations do **not** supply the native common phase, the native
target-CSF phase matrix, spin-label mapping or native normalization. They are
rejection controls after those conventions are fixed, not a replacement for
direct signed CSF agreement. Moving the continuum creation operator past a
fourteen-electron target has even total parity, but internal spin-orbital sorting
still needs its own determinant-dependent signs.

## Boundary

The initial source subset resolved row order, coefficient axes, prototype role
and sorting behavior but contained only the `CIRMAT` call. Subsequent inspection
of `ci-diag/scatci_routines.f` (SHA256
`cbbaa8e32026f3ca345c2b382ed38cdb4868a391063a963900de59e4a99f9ea6`)
resolves `CIRMAT`/`CPAKPZ` and `MKORBS`; the orbital-table wrapper is pinned at
`63fa977571f307432d6e8d2afe318c6a790ad0b74f176e5bf9336d698d7cc98d`.
These are restricted source identities, not source payloads shipped here.

`native_csf_selection.selected_csf_determinants` decodes bounded selected CSFs
from the six-record int64/float64 electron-only D2h format. It preserves raw
reference-slot substitutions and converts directly to alpha-then-beta canonical
determinants; it does not apply a native sorting sign a second time. Unselected
CSFs remain unexpanded read-only mappings. Explicit phase/group metadata are
returned separately, and `corrected_target_coefficients` applies only the declared
row-phase product. Hand-built format, slot-parity, phase-covariance and resource
controls precede molecular use. This is not a coupled-CSF adapter; the signed
selected-prototype witness above remains an independent admission gate.

Even after that gate, complete interior Q coefficients, a consistently computed
cross-frame one-body metric, spin-changing removal contractions and coherent
PP/PQ/QP/QQ assembly remain necessary. The generated real-CI benchmark establishes
one ingredient's scale only. It does not authorize native recovery, other-R
acquisition, another geometry, refitting or promotion to qscat.
