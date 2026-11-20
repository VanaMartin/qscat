# Coherent core-sector interior overlap

`core_sector_full_overlap` is a generated-stage molecular-dimension consumer of
the [four-block identity](co-full-interior-metric.md) and
[occupied-spin-minor kernels](co-core-sector-contractions.md). It calls no native
program, electronic solver or orbital integrator and supplies no missing vector.
One-body metrics are dimensionless overlaps on a declared common physical domain.

## Explicit products and one Q root

`CoreInteriorExpansion` contains one bound-only `CoreSectorCI` root Q and up to
two `TargetProducts` families. Each family contains target CI roots at one spin
projection, its explicit continuum spin and a target-by-spatial-continuum complex
amplitude array A. Products are in `a†(continuum) Phi` order. Native target/CSF
phase conversion belongs in the inputs, not a guessed coupling sign inside the
overlap evaluator; follow the
[native map admission contract](co-native-continuum-product-map.md).

For a doublet with total M=+1/2, an M=0 target plus continuum alpha and an M=1
target plus continuum beta contribute coherently. Their target spin-resolved
electron counts differ, but each target-plus-continuum product must have exactly
the Q root's counts. Combine roots with the same continuum spin in one family;
do not separately normalize or square those families. The kernel checks counts,
array shapes, finite coefficients and inherited core/active/root/work caps.

## Four spatial metric blocks and contraction

Let bb be the spatial bound cross metric owned by `CoreSectorMetric`, and let
bc, cb, cc be bound–continuum, continuum–bound and continuum–continuum spatial
cross metrics. Spin labels supply Kronecker deltas; mixed spatial metrics are not
assumed zero across frames. bc and cb are not assumed equal or adjoints within
one directed cross-frame call.

For a left family f with continuum spin s and right family g with spin t,
write O_IJ for target overlaps and

    R_mu,I,nu,J = <D_(mu,t) Phi_I^left | D_(nu,s) Phi_J^right>.

The deletion spin on each target is the **other side's continuum spin**. Thus
s!=t generally requires spin-changing removal, even though the one-body metric
is spin preserving. Extend `removal_overlaps(..., spin=t, right_spin=s)` to
retain those sectors; its default right_spin=spin preserves existing callers.
Parent target projections must not be compared before deletion and used to
drop a nonzero exchange block. Core-hole sectors remain explicit on both sides.

Project channel amplitudes before contracting the removal tensor:

    L_I,nu = sum_a conj(A_left[I,a]) cb[a,nu],
    B_J,mu = sum_b A_right[J,b] bc[mu,b].

This family's PP contribution is

    delta_s,t sum_I,a,J,b conj(A_left[I,a]) O_IJ cc[a,b] A_right[J,b]
      - sum_I,mu,nu,J L_I,nu R_mu,I,nu,J B_J,mu.

Each removal tensor is built under the existing byte/work caps, contracted and
released before moving to the next family pair. No complete continuum-configuration
Gram matrix or cross-core inverse is formed. A real target CI with complex Q/A
does not require converting every target-root tensor into a complex array.

PQ contracts L against `<Phi_I|D_(nu,s)Q_right>`. QP contracts the corresponding
right projection against `<D_(mu,t)Q_left|Phi_J>`. The Q root has one extra
electron, so each deletion lands in the matching target sector; core deletion
must still be included. QQ is the ordinary full occupied-spin-minor overlap.
The result is the same `qq/pq/qp/pp/total` complex-amplitude dictionary as the
small four-block oracle, with `total=qq+pq+qp+pp`, not a sum of sector probabilities.

## Admission and limitations

Small independent full spin-orbital determinant/exterior controls compare every
block, including complex rectangular continuum metrics, singular core blocks,
nonzero different-spin exchange, reversal, target-basis/whole-state phase
covariance and a complete-state norm with nonzero interference. The cross-spin
weighted deletion trace is checked against explicit S-minus on the ket, not
inferred from a same-spin particle-count identity.

A generated 14/15-electron benchmark must separately establish actual time/RSS
for the intended real-target/complex-channel/Q types, root counts and continuum
inventory. Inherited caps alone are not a cost forecast; 48 complex target roots
can exceed the left-removal cap even when real roots fit. Keep Python as the
oracle; there is no promotion to qscat or Rust optimization without its own gate.

A native consumer still needs complete coherent P/Q coefficients, source-matched
native phase/order conversion, four consistently integrated orbital cross blocks
and complete-state self normalization. A passing generated kernel or sampled
native phase map is not a measured molecular endpoint overlap, unsampled branch
continuation, crossing location or electronic accuracy. No acquisition or
eigensolve is authorized by this API.
