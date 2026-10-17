# Physical boundary matching and outgoing poles: a radial shell benchmark

## Purpose and scope

Before using the saved CO boundary data in an independent outer propagator,
test the scattering boundary and outgoing-pole conventions on a self-adjoint
Hamiltonian whose poles have an independent closed-form oracle. This is the
first physical toy-model gate after the [saved-matrix inventory](co-saved-scattering-matrices.md).
The coupled polar tail and native boundary-amplitude normalization require later
gates. No CO file or phase fit enters this experiment.

Lifecycle: design → Python toy → numerical validation here. Optimization is
unwarranted for this small benchmark; promotion is not requested.

## Hamiltonian and exact oracle (derived here)

Use atomic units, electron mass one, reduced radial s wave on the half-line:

\[
 H=-\tfrac12\frac{d^2}{dr^2}+g\delta(r-a),\quad u(0)=0,\quad
 u(a^+)=u(a^-),\quad u'(a^+)-u'(a^-)=2g u(a).
\]

For real positive `g`, the real derivative jump conserves the radial probability
current `Im(u* u')`; the operator has no negative-energy bound state. It has
outgoing resonances: the repulsive shell temporarily confines a cavity wave.
Let `x=r/a`, `q=ka`, `λ=2ga`. Inside, the regular solution normalized by
`u_x(0)=1` is `sin(qx)/q`; at the shell `u_x` jumps by `λu`. Outside write
`u = C(S e^{iqx} − e^{-iqx})/(2i)`. Thus

\[
 F(q)=1+\lambda e^{iq}\frac{\sin q}{q},\qquad S(q)=F(-q)/F(q).
\]

`F(0)=1+λ`, so the apparent zero at `q=0` obtained after multiplication by `q`
is **not a pole**. For nonzero q, `F(q)=0` gives
`e^{2iq}=1−2iq/λ`. With `z=2iq`, put `w=λ−z`; then
`w exp(w)=λ exp(λ)`. The outgoing positive-real-q branches are

\[
 q_n=\frac{i}{2}\{W_{-n}(\lambda e^\lambda)-\lambda\},\quad n=1,2,3,\ldots;
 \qquad E_n=\frac{q_n^2}{2a^2},\quad\Gamma_n=-2\,\mathrm{Im}\,E_n.
\]

The principal Lambert-W branch yields the removable q=0 solution and is excluded.
The benchmark oracle evaluates these functions at 60 and 90 decimal digits;
agreement of the two precisions and the unmultiplied Jost residual are checked.
Lambert W determines the pole independently of the numerical shooting/root solve.
The real-axis oracle uses the exact Jost ratio. Neither uses ODE integration.

## Independent numerical boundary route

Integrate the regular initial-value problem `u_xx = −q²u` from x=0 to 1 with
`(u,u_x)=(0,1)` using DOP853. Apply the exact distributional shell jump, then
integrate the exterior to each prescribed matching radius `t = B/a`.
There is no finite-width approximation to the delta shell.

Keep the homogeneous value/derivative pair rather than its logarithmic derivative:
the latter is singular at wavefunction nodes. For dimensionless Bloch coefficient b,

\[
 d=t u_x-bu,\quad R_b=u/d,\qquad
 S=e^{-2iqt}\frac{d+(b+iqt)u}{d+(b-iqt)u}.
\]

This explicitly defines the **benchmark's** dimensionless R convention,
`u(B)=R_b[B u'(B)−bu(B)]`. Storing `(u,d)` avoids divisions at R poles or zeros;
a common normalization of that pair cancels. The equivalent length-valued
R for b=0 is `u/u' = B R_0`. This derivation does not establish that native
UKRmol boundary amplitudes implement this convention: trace the native
writer/assembler and Hartree/Rydberg factors separately before CO use.

The outgoing shooting residual is `J(q)=exp(iqt)[u_x−iqu]`, which equals F(q)
for the stated regular normalization. Its sign specifies `exp(+ikr)`; searching
in q avoids an ambiguous library square root in energy. Accepted roots must
have `Re q>0`, `Im q<0`, `Re E>0`, and Γ>0. They lie on the continuation of the
positive-k scattering branch into the lower k half-plane. The incoming residual
at the outgoing pole must be nonzero, ruling out an S numerator cancellation.
This is a single-threshold sheet test, not multichannel sheet qualification.

## Frozen validation matrix and stopping rules

`projects/ukrmol_co/physical-pole-contract.json` freezes:

- shell radius a=5 bohr; λ=0 (free), 2, 10, 30;
- real q grid: 101 points from 0.1 to 10.1, plus π/2π/3π wavefunction-node
  controls and q=1e-5, 2e-5, 4e-5 threshold controls (107 distinct points);
- matching radii B/a=1,2,4 and relative ODE tolerances 1e-5,1e-8,1e-11,
  absolute tolerance 0.01 times relative tolerance; **36 real-axis profiles**;
- b=−1,0,2 and common boundary rescaling 7, testing representation invariance;
- nine physical poles (λ=2/10/30 × n=1/2/3), each at all three radii and all
  three ODE tolerances: **81 shooting solves**. Each solve starts once from
  `q=nπ−0.3−0.2i`, without the oracle; root tolerance 1e-9, maxfev=80.

Store every root solver exit/message, residual, oracle discrepancy, and all
three refinement levels. Finest real-axis S error ≤1e-7, phase error modulo π
≤1e-7 rad, unitarity `||S|²−1|` ≤1e-10; finest pole q error ≤1e-7, position
and full-width errors each ≤0.1 meV. Root residual ≤1e-7 and noncancelling
incoming residual magnitude ≥0.1. Bloch/rescaling errors ≤1e-10, maximum
finest radius spread in q ≤1e-7. Real-axis and pole maximum errors must decrease
at least tenfold from coarsest to finest unless coarsest already lies below
1e-10. Intermediate levels remain reported; monotonicity is not assumed near
roundoff. Low-energy `−arg(S)/(2q)` must reproduce `a_s/a=λ/(1+λ)` within 2e-8
at the finest tolerance. Free propagation must reproduce S=1 without shell poles.

Mandatory negatives reverse the shell jump or multiply it by 1/2 at q=1.3,
λ=10: both must differ from the exact S by >0.01. An incoming-boundary residual
evaluated at the outgoing oracle pole must be >0.1. This detects sign/factor
errors; real-axis unitarity alone cannot, since a wrong real shell is also unitary.

Any missing/nonfinite data, failed ODE/root status, or false/unset gate stops
with nonzero command exit after writing the available evidence. No reseeding,
tolerance relaxation, channel addition, native run or automatic successor occurs.
All root statuses are required at every rung; the accuracy budgets apply to the
finest rung and the declared convergence sequence. The physics parameter grid,
initial guesses and numerical gates are not tuned using results of this matrix.

## What a pass releases

A pass qualifies this scalar-shell implementation and the defined boundary
convention in the tested regime. It permits designing the next coupled polar-tail
benchmark, including angular mixing, closed channels, threshold square roots and
tail/outgoing controls. It does not validate a CO pole, infer an uncertainty bound
from the 2048 exported poles, fix the five-A2 label mismatch, qualify a production
sweep or replace the existing extraction rejection.
