# Independent real-axis reconstruction of the saved native CO K matrix

## Scope and boundary conventions

The [bounded combined-exterior control](../../projects/ukrmol_co/COMBINED_EXTERIOR_RESOLUTION.md)
qualifies its defined toy DtN problem. This next method reconstructs real-axis
open-open K from saved zero-Bloch surface sums, channel thresholds and both
multipole matrices. Its output must be compared entry by entry with the saved
native matrix; phase sums and unitarity alone are insufficient. Retain all 298
channels, including the 292 closed ones, and all 2048 exported poles. The
truncated inner spectrum and electronic model are inputs, not newly converged
approximations.

Use the [derived native surface convention](co-native-boundary-conventions.md):
`u(a)=Q(E)u'(a)`, `Q=Σwwᵀ/(E_n−E_total,Hartree)`, a=18 bohr. Incident energy
is the saved requested-eV value times0.0735 Ryd/eV, total Hartree energy is the
first target energy plus half the incident Rydberg energy. With Z=0 the exterior
equation is

`u''=V(r)u`, `V=diag(threshold_Ryd−incident_Ryd)+C2/r²+C3/r³`,

where `C2=diag(l(l+1))+CF1`, `C3=CF2`. The exported coefficients already contain
the Hartree→Rydberg factor. No additional half, radius division, flux factor or
multipole doubling is fitted to native data.

## Source-derived asymptotic sign and flux normalization

The read-only engine's `outer/source/libouter/cfasym.f` has SHA256
`28ef94214b940246938697ff91eaef3523a7deb101215c53c82e1188019725b3`.
`GAILIT`, lines135–162 selects k=sqrt(abs(channel energy)); `ASYWFN`,
lines304–320 divides the Coulomb/Bessel values by sqrt(k) and multiplies their
argument derivatives by sqrt(k). Lines413–419 form the regular/irregular
coupled columns. `BKSHWF`, lines1614–1633 explicitly fixes the neutral limits
`f=sin(kr−lπ/2)/sqrt(k)`, `g=cos(kr−lπ/2)/sqrt(k)`.
The closed branch (`ASYWFN`, lines459–469, and `BKSHWF`, lines1650–1653)
selects a decaying column; its finite nonzero column scaling cannot affect
open-open K. `KMAT_compak`, lines2204–2229 and2260–2269 solves

`[g−Q g', h−Q h'] [K; c] = −(f−Q f')`.

Here h contains all closed decaying columns, and g only the open irregular
ones. Thus the physical solution is `f+gK+hc`, and
`S=(I−iK)^−1(I+iK)`. These literal anchors fix conventions; the Python method
does not reproduce Gailitis recurrences or continued fractions. A source trace
does not yet validate native finite-radius asymptotic errors.

## Independent coupled inverse-power tail

For each reference column j let `s_j=+ik_j` for an open outgoing column, or
`s_j=−κ_j` for a closed decaying column. Write

`U_j(r)=exp(s_j r) Σ[p=0..P] a_p^(j) r^(−p)`, `a_0=e_j`.

For `D=diag(threshold−incident)` substitution gives

`(D−s_j²I)a_n = −2s_j(n−1)a_(n−1)`
`                  +[(n−2)(n−1)I−C2]a_(n−2)−C3 a_(n−3)`.

Outside the threshold-degenerate subspace this determines a_n directly. Within
it the next equation determines

`a_p = ([p(p−1)I−C2]a_(p−1)−C3 a_(p−2))/(2s_j p)`.

Threshold grouping tolerance is an explicit input and its maximum within-group
spread must be recorded. No energy at an exact channel threshold is supported.
Open columns carry the fixed phase/amplitude `exp(−il_jπ/2)/sqrt(k_j)`; their
imaginary and real parts give f and g. Closed columns are rescaled to omit the
common exponential at the outer endpoint, preventing underflow without changing
the permitted solution space. Evaluate their values and analytic derivatives. Series
order and endpoint are independent controls; no tail-cutoff approximation is
silently substituted for this expansion.

## Stable ratio Numerov propagation

On a uniform real grid define `B_i=I−h²V_i/12`, `z_i=B_i u_i`. The matrix
Numerov equation is

`z_(i+1)=(12 B_i^−1−10I)z_i−z_(i−1)`.

Propagate the inward-value ratio outward,

`R_i=z_i z_(i+1)^−1=[12 B_i^−1−10I−R_(i−1)]^−1`.

This avoids propagating exponentially growing closed-channel fundamental
columns. Start with `u(a)=Q`, `u'(a)=I`, and a sixth-order Taylor evaluation at
a+h. Taylor derivatives follow the exact smooth potential:
`u^(n+2)=Σ[j=0..n] binomial(n,j) V^(j) u^(n−j)`.
The initial ratio is `B_0 Q (B_1 u(a+h))^−1`. All matrix divisions are linear
solves; no symmetrization is applied to a measured K.

At the endpoint use the same sixth-order Taylor identity to write
`u(R−h)=A_u u(R)+A_d u'(R)`. Define
`L(X)=[B_(N−1)A_u−R_(N−1)B_N]X+B_(N−1)A_d X'` and solve

`[L(g),L(h)] [K;c]=−L(f)`.

Equilibrate matching columns by their maximum absolute entry before the linear
solve, then undo those scales on the returned coefficients. Record both raw and
equilibrated condition estimates. This makes arbitrary closed-column scaling a
numerically meaningful invariance test without changing the physical references.

Matching values and analytic derivatives is essential: a truncated tail series
does not exactly solve the differential equation at finite R. Substituting its
values at two neighbouring nodes instead would introduce a step-dependent
boundary residual and spoil the fourth-order convergence test.

The expected global spatial error is fourth order. Singular ratio charts,
ill-conditioned solves, nonfinite values and failed controls are recorded and
stop a frozen study without changing the grid. This real-axis method does not
authorize complex-pole searches or general nonzero-Bloch propagation.

## Validation before native comparison

Pure-Python functions live in `projects/ukrmol_co/exterior_numerov.py`:
`tail_columns` evaluates finite-order coupled references;
`numerov_k` returns the unsymmetrized open K and matching diagnostics.
Tests require free centrifugal references and their sign/flux Wronskian,
unequal-threshold coupled constant-coefficient matching, and fourth-order
convergence against independent direct ODE integration on a smooth coupled
inverse-square/inverse-cube finite interval. Closed column rescaling must leave
K invariant. Native inputs are not used for these analytic/numerical tests.

Before a native comparison, freeze exact input hashes, selected energies,
independent step/endpoint/series-order families, K/S/phase/reciprocity budgets and
resource caps. Measure a bounded synthetic 298-channel timing case on Sadaharu
to justify the cost. The first native comparison is a bounded equilibrium-B1
pilot; a pass cannot release the six-case/594-matrix reconstruction by itself.
