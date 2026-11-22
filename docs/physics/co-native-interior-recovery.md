# Reconstructing a retained native interior state without changing its model

Full-vector recovery repeats a native Hamiltonian build/eigensolve and needs an
explicit release. Saved integrals, targets and CSFs are reused; no new electronic
prescription is implied. Retaining more coefficients is not retaining more roots.

## Contracted rows and Q tail

The source-resolved native order is target sector, target state, continuum
orbital. The executed `MCONT`, `NOTGT`, `NUMTGT`, target selections and explicit
MO symmetry/continuum flags must agree. Select continuum coordinates in native
MO order, not by sorting orbital energies or treating prototype pairs as spins.
`native_interior.continuum_row_layout` validates this explicit inventory.

The full retained state has continuum rows P followed by uncontracted bound-only
CSFs Q. `assemble_native_interior` maps P through the signed raw-target product
adapter and Q through the already qualified sparse bound-CSF transform. Their
coherent sum is represented by `CoreInteriorExpansion`. Raw target roots carry
no second loaded phase delta. No missing Q is manufactured by the adapter.

A zero-Q input may test a supplied P component but is not a complete physical
state. Such a control must be labelled partial. Physical normalization uses the
complete four-block self metric on the same domain as the cross comparison.

## Recovering the same retained spectral combination

Let old complete native eigenvectors be U (columns), saved continuum rows be
P U, and accepted matching coefficients be c. The desired interior state is U c,
not a nearest inner root. A recovered eigenbasis V can change signs or rotate
within degenerate eigenspaces. Applying c directly to V would change the state.

`recover_retained_state` uses unchanged energy-ordered roots and old/new saved
continuum components. In each prospectively bounded energy cluster, orthogonal
Procrustes gives R such that old continuum rows equal R times recovered rows
(the implementation stores roots as rows). It verifies energy agreement,
identifiability of the visible cluster, orthogonality and the component residual.
The matching coefficients transform as c_new = R.T c_old, without conjugating
this real basis transformation. The complete state is V.T c_new. Complex matching
coefficients are retained literally.

An invisible or rank-deficient cluster is a rejection, not permission to guess a
phase, discard a component, refit c or alter the finite retained space. Energy
proximity alone is insufficient. Root/cluster bounds precede matrix allocation;
algorithmic tolerances do not substitute for a frozen experiment's actual gates.
Numerical native eigenpair residuals and the complete-state norm/boundary checks
remain independent admission requirements.

These operations establish representation and basis correspondence, not unique
continuation through an unsampled interval, electronic accuracy or an extension
of the local CO potential. The [common-ball metric](co-cross-frame-orbital-metric.md)
and [native spin/phase contract](co-native-continuum-product-map.md) remain separate.
