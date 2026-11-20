# Core-sector CI contractions without a core-block inverse

This is a generated-control and performance-stage implementation for the
[four-block overlap](co-full-interior-metric.md). It does not decode native
target/continuum product phases or supply missing interior-state coefficients.
Native source/addressing and spin-coupling conventions require separate checks.

## Representation and canonical signs

`CoreSectorCI` fixes a spatial partition into core coordinates followed by active
coordinates. It stores explicit occupied-core tuples separately for alpha and
beta spin, active electron counts, and root-by-alpha-string-by-beta-string CI
coefficients. Active strings are ordered by ascending integer bitmask, including
the unique vacuum string. Canonical full determinant order is **all alpha
spatial orbitals, then all beta spatial orbitals**. Coefficients must include
the parity converting any different native ordering into this convention.

For spin s, all occupied coordinates of its fixed core precede its active
coordinates. `remove_spin_orbital` applies coefficient deletion with the index
of that orbital in the complete canonical occupied tuple. Thus beta deletion
includes the number of alpha electrons before the beta block. Core deletion
creates an explicit core-hole sector; active deletion reduces an active count.
The function returns `None` only for an exactly Pauli-absent orbital, not for
missing native data or an unsupported physical sector.

`spin_lower` applies the unnormalized operator sum_p a†_(p,beta) a_(p,alpha)
to a paired, fully occupied frozen core plus active CI. Pauli exclusion makes
core contributions zero only in that paired-core case; other core-hole inputs
are rejected rather than silently dropping terms. This is spin algebra, not
a native Clebsch-Gordan phase convention or automatic state renormalization.

## Occupied-spin minors and batching

For a spin-preserving spatial cross metric S, each full determinant overlap is

    det(S[occupied_alpha_left, occupied_alpha_right])
      * det(S[occupied_beta_left, occupied_beta_right]).

Left/right fixed core tuples can differ, and the active electron counts can
differ as long as the total count in each spin agrees. The corresponding spin
minor matrices can therefore be rectangular. When total spin-resolved counts
disagree, that overlap is exactly zero in this spin-preserving metric. A
spin-changing **removal** overlap need not vanish when its parent target spin
projections differ: compare counts after each removal, not before.

The occupied lists retain every core/active cross entry. `CoreSectorMetric`
forms these full occupied-spin determinants in bounded chunks. It never inverts
the core cross block, so a singular core block does not by itself reject the
observable. For nonsingular well-conditioned cores the result can be checked
against the independent Schur-complement frozen-core formula.

If A and B are the two spin-minor matrices and c_R one right CI tensor, its
transformation is A c_R B^T. Left contraction uses c_L†. No complex conjugation
belongs on B^T; bra conjugation is applied once to the CI coefficients.
`overlaps` processes eight right roots at a time and caches occupied-spin minors
under a checked byte cap. The spatial matrix is copied and made read-only so a
changed caller array cannot invalidate those caches.

`removal_overlaps` assembles one-spin orbital/target removal blocks. Left
removals are grouped by core-hole/active-count sector, right removals are visited
one orbital at a time, and each right root block is transformed once for a
left-sector family. It does not recompute that transformation for every left
orbital in the same sector. It checks the total left-removal allocation before
building arrays; the right arrays and eight-root transforms are bounded scratch.
An explicit `right_spin` additionally supports spin-changing deletion; the
default retains same-spin behavior. The
[coherent four-block consumer](co-core-sector-full-overlap.md) specifies the
crossed deletion spins required by exchange between different continuum spins.

## Scope, controls and limits

The caps are eleven active orbitals, four core coordinates, eight occupied
orbitals per spin, 64 roots per CI sector, at most one million pairs per spin
minor matrix, a 64-MiB minor cache and a 1-GiB left-removal storage budget. A
molecular-size benchmark must report whole-process RSS and CPU/wall time, not
infer them solely from these array caps. Real target coefficients and complex
resonance amplitudes have different memory/cost requirements; a real benchmark
does not establish a complete complex native workload's performance.

Small controls compare with direct full spin-orbital determinants and explicit
fermionic deletion/creation, including complex nonorthogonal metrics, singular
core blocks, nonzero core-hole/active-hole overlaps, singlet/triplet lowering,
spin-changing removal, root batching, cache eviction, vacuum and work/memory
rejection. No orbital/core/state projection or normalization is used to hide a
failed overlap. A generated fourteen-electron benchmark can establish this
kernel's measured representation range, not native state identity or physical
model accuracy.

The coherent PP/PQ/QP/QQ assembly still needs native target-product mapping,
native spin/phases and complete Q/P coefficients in the same declared physical
domain. This module does not read molecular payloads or authorize their
acquisition, an eigensolve, another geometry, Rust optimization or promotion to
qscat. Core-hole sectors must remain explicit in any later consumer.
