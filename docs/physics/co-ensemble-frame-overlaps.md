# Physical retained-subspace comparison of CO orbital ensembles

## Question and limits

At R=2.1323bohr, compare the original CAS(10,11)/cc-pVDZ forty-component orbital
average with the independently qualified fifty-component average. Both proposed
scattering prescriptions retain the lowest five target roots per spin/irrep.
Equal root counts and property IDs do not establish that they retain comparable
many-electron states. This is a finite **retained-subspace dominance screen**,
not electronic convergence, global state continuity or resonance identification.

Reuse the authenticated old twelve-root CI arrays at the exact original
scattering checkpoint. Generate only eight new fixed-orbital CI sets:13 A1 roots
and8 roots otherwise, per spin, at space160 /600 cycles. These new roots already
have independent coverage; the extra computation supplies the CI arrays that
were not saved by that coverage. No orbitals are optimized, geometries changed,
native engines called or scattering Hamiltonians assembled.

## Physical overlap

Use atomic units, the physical AO cross overlap and the two original orbital
frames, not AO-following Lowdin transport. For core-plus-active MOs A and B,
`M = A† S_AO B`. CI arrays follow ascending bitstring determinants. The
many-electron overlap is the bra-conjugated CI contraction with the product of
alpha- and beta-occupied Slater determinants of M, including the two doubly
occupied frozen core orbitals and core/active cross terms.

The existing `ci_overlap.frozen_core_overlaps` uses the Schur identity. With
core block C and active Schur complement T,

`det M_occupied = det(C) det(T_active-occupied)`,

so the combined spin prefactor is **det(C)²**, not |det(C)|². Root bra coefficients
are conjugated separately. Compare this result with explicit full occupied
minors, independently enumerated by integer bitstrings. Also compare reverse
overlaps with the conjugate transpose, check own-frame CI/MO orthonormality,
and bound physical overlap magnitudes and singular values by1+1e-9. Retain all
raw arrays, numerical differences and checks; no row/projection normalization.

Numerical budgets retain the previous overlap method:1e-10 full-minor and
reciprocity discrepancies,1e-9 own-frame Gram error,1e-12 AO cross-metric
reciprocity. New CI roots retain1e-9-Hartree physical/penalized/penalty-action
norms,1e-6 spin error and1e-7-Hartree agreement with the qualified fixed frame.
The original cached twelve-root CI gates and native first-five agreement must
reconstruct before reuse. Buffers are diagnostics, not new retained channels.

## Subspaces rather than arbitrary root gauges

For each sector save the full12-by13 or12-by8 matrix O. Let A=O[:5,:5]. Since
the two retained root sets are orthonormal, the eigenvalues of A†A are squared
principal overlaps: the smallest is the least captured weight over **every**
normalized combination in either retained subspace. This quantity is invariant
under unitary mixing or phases inside either retained set.

Prospectively require **min(svd(A))²≥0.5 in all eight sectors** for the dominance
screen. The0.5 scale is the earlier single-root weak-overlap scale applied to
the worst retained-subspace direction: more overlap is retained than missing.
It is deliberately not a near-equality criterion or bound on scattering error.
Lost physical overlap can reflect core/active-frame changes and finite root
manifolds; it is not an omitted-state population estimate.

Separately report bidirectional one-to-one assignments of the retained five
roots against all available opposite-frame buffered roots, with the existing
weak-score<0.5 /best-minus-second margin<0.1 /retained-to-extra triggers. Do not
let high buffered roots compete as additional mandatory assignment rows. These
individual-root flags are diagnostics, not a replacement for the invariant
subspace screen: an exact unitary rotation within the retained set must pass
the subspace screen even when individual-root matches are ambiguous.

Mandatory toy controls are retained permutations and unitary mixing (pass), a
retained/buffer exchange (reject), and a uniformly contracted0.7 overlap (reject,
because0.49 is not renormalized to1). These test the decision rule independently
of the CO outcome. Numerical CI/overlap failure stops immediately. A scientific
dominance rejection preserves all eight sector results and exits1, with no new
root-count, space, geometry or tolerance ladder.

## Cost and next decision

The previous three-anchor24-sector /48-matrix study took4322.37s;4298.41s was
inside its sector CI generation/validation. The remaining23.95s also included
metrics, self checks, overlaps and setup, so it is not an isolated overlap timing.
The new frame's already completed48-probe scan took308.66s using symmetry-aware
residual actions. Reusing all old waves and generating only eight new sets is
forecast **300s expected /900s upper /1800s cap**, four threads on reserved cores
0–3 /4GiB. Persist new waves so subsequent work does not repeat these solves.

A pass supports designing a matched **retained-subspace** ensemble contrast;
individual-state and global-identity claims remain limited by their flags.
A rejection means the same ordinal prescription changed the physical retained
subspace too much for this screen; report that measured model change rather than
silently choosing other states. Either outcome informs compatible continuum/
sparse-preparation and resource design. No new inner solve or production sweep
is automatically released.

## Python interfaces and evidence

`ensemble_frame_overlaps.diagnose_overlap(overlap)` accepts an unscaled physical
root-overlap matrix (at least five roots on either side); it returns squared
principal overlaps, both retained/buffer projection lists and both directed
retained-only assignments. Nonfinite or noncontractive matrices reject.
`full_minor_overlaps(left, right, metric, ncore, electrons)` independently builds
full occupied minors using ascending integer-bitstring enumeration; it does not
use the Schur complement. Its complex toy is checked against a scalar Slater
oracle and the existing primitive. No compiled optimization or library promotion
is warranted for this bounded campaign.

`ensemble_frame_comparison` is invoked as a module with `--root`, `--output`
and `--contract`. Its prospective machine contract pins32 prerequisite files and
the two unchanged helper hashes. It authenticates the cached-wave provenance,
both checkpoints/configurations, original native energies and the qualified
new-frame/property proof. Saved old residual gates reconstruct from authenticated
metadata; they are not new residual-action calculations. Every new root receives
fresh physical/penalized/penalty-action checks via the byte-identical
`coverage_reference.state_residual` helper. Raw CI vectors are not renormalized.

Save ten physical AO/MO/self-metric arrays, eight new CI arrays, eight logs/root
records and three overlap routes per sector (Schur, full minors, reverse), plus
numerical/diagnostic JSONs, progressive/final status and any failure record.
State-mixing flags and the original seed rejection remain explicit beside the
invariant dominance verdict. Publish generated artifacts and enough authenticated
prerequisites to reproduce overlaps without another CI solve.
