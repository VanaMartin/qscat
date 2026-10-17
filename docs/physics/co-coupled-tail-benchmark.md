# Coupled polar-tail and threshold-sheet benchmarks

## Purpose and scope

The [scalar physical pole benchmark](co-physical-pole-benchmark.md) passes.
The next two controls isolate matrix angular tails and unequal target thresholds
before they are combined in a native CO propagator. They are physical radial
Hamiltonians with independent analytic oracles; no native export or fitted phase
is input. Atomic units apply, mass one and a=5 bohr; x=r/a, ε=2a²E.
This remains Python toy/validation work, with no optimization or promotion.

## A. Exact inverse-square dipole tail

On x≥1 impose a hard inner wall u(1)=0 and solve

\[
 -u_{xx}+C u/x^2=q^2 u,\qquad
 C=\begin{pmatrix}2&-2d/\sqrt5\\-2d/\sqrt5&6\end{pmatrix}.
\]

This is the l=1,2, |m|=1 two-wave projection of centrifugal energy plus
`−d cosθ/r²`: the angular recurrence gives
`<1,1|cosθ|2,1>=1/√5`, with zero diagonal elements. Truncation to two waves
defines the model being tested, not an angular-convergence statement about CO.
The origin is excluded by the hard wall; for d=0/0.4/1.2 the matrix is positive.
Both channels have the same threshold, so a constant real orthogonal O diagonalizes
C. The eigenchannel orders are `ν_j=sqrt(c_j+1/4)` and are noninteger when d≠0.
The exact outgoing/incoming matrices are

\[
 H^\pm(x)=O\,\mathrm{diag}\left[
 \sqrt{\pi x/2}\,e^{\pm i(\pi\nu_j/2+\pi/4)}
 H^{(1,2)}_{\nu_j}(qx)\right]O^T
 \mathrm{diag}[e^{\mp i\ell\pi/2}].
\]

They approach `diag(exp[±i(qx−lπ/2)]/sqrt(q))`, fixing unit flux on the real
axis and the physical partial-wave phases. For complex q in the right half-plane
use Hankel1/Hankel2 analytically, not complex conjugation. The hard-wall solution
written `u=H+ S−H−` gives **`S=(H+(1))⁻¹ H−(1)`**. This direct boundary value,
evaluated with independent mpmath eigensystems/Hankel functions at 60/90 digits,
is the oracle. It requires no radial propagation.

The numerical path integrates the full 2×2 fundamental pair U(1)=0,V(1)=I
outward in the original angular basis using DOP853, then solves
`[U, −H+; V, −H+'] [X;Y] = [H−;H−']`, returning **S=−Y**. It uses SciPy
eigensystems/Hankel functions only at the matching boundary. Thus the primary
comparison tests full matrix propagation and matching independently of exact
radial evolution; both routes still share the derived analytic tail model.
Zero-dipole controls reduce to integer-order hard-sphere partial waves.

At a separate finite cutoff B, match the propagated pair to **d=0** centrifugal
Hankel functions. That deliberately drops the dipole tail beyond B. Its errors
against the infinite-tail oracle must be visible and diminish under B/a=8→32→128.
This measures why a free/uncoupled boundary is insufficient, not a convergence
bound for the actual CO multipoles. No claim about CO quadrupole/polarization is made.

## B. Coupled shell with an open/closed threshold

For x>0, u(0)=0, solve a two-internal-channel radial Hamiltonian

\[
 -u_{xx}+T u+\Lambda\delta(x-1)u=\epsilon u,
 \quad T=\mathrm{diag}(0,16),\quad
 \Lambda=\begin{pmatrix}10&c\\c&12\end{pmatrix},\quad c=0,1,3.
\]

T and Λ do not commute when c≠0. The real symmetric shell gives the jump
`V(1+)=V(1−)+ΛU(1)`, conserves current and couples different channel momenta.
This model uses radial s channels to isolate threshold elimination; it is not
the angular tail model A and does not combine unequal thresholds with a dipole tail.
At a real energy below 16 choose `k0=+sqrt(ε)`, `k1=+i sqrt(16−ε)` (decaying);
above 16 choose both positive real. Exact threshold points, where flux-normalized
open waves are singular, are excluded; points ±0.001 either side are retained.

The independent exact shell logarithmic derivative is
`L=diag(k_i cot k_i)+Λ`. Below the second threshold eliminate the closed channel:
`L_eff=Loo−Loc(Lcc−i kc)⁻¹Lco`. Above it use the full L. With open-channel
`P±=diag(exp(±ik)/sqrt(k))`, the exact unit-flux S is

\[
 S=(P^+)^{-1}(L_{\rm eff}-iK_o)^{-1}(L_{\rm eff}+iK_o)P^-.
\]

The oracle uses mpmath trigonometry, Schur elimination and matrix algebra at
60/90 digits. The numerical path integrates U(0)=0,V(0)=I, applies the shell
jump and propagates to the declared matching radius. Its matching system admits
only incident open waves, outgoing open waves and **decaying closed** waves.
Only open columns carry `1/sqrt(k)` flux factors. It does not use the oracle's
logarithmic derivative or Schur complement. Report the linear system condition
number to expose closed-channel growing-solution conditioning.

For an outgoing resonance use q=k0 as root variable, ε=q², and explicitly keep
`k1=+i sqrt(16−q²)` on the continuation below the second threshold. For Re ε<16,
Im ε<0 this branch has Im k1>0 and generally Re k1<0; flipping its sign gives a
different sheet. The numerical Jost determinant is
`det[diag(exp(ik_i B))(V−iK U)]`. Its independent oracle uses exact
`U=diag(sin k/k)` and `V=diag(cos k)+ΛU` at x=1, followed by a high-precision
root solve from a fixed secant pair. For c=0 its first-channel root must reproduce
the earlier independent Lambert-W pole. Changing the sign of k0 at the oracle
pole must give a nonzero determinant, excluding numerator cancellation.

## Frozen matrix and budgets

The JSON contract fixes d/c values, all energies, matching radii, seeds and gates
before either complete study is run:

- Tail: five real q and two complex q, three d, B/a=1/2/4, and ODE rtol
  1e-5/1e-8/1e-11: **189 matrix samples**, 135 real and 54 complex.
- Threshold: twelve real ε, three c, B/a=1/1.5/2, same three tolerances:
  **324 matrix samples**, with one/two open channels either side of ε=16.
- Threshold poles: one first-channel resonance per c at all radii/tolerances:
  **27 numerical root solves**. Each starts once at q=π−0.3−0.2i; no oracle
  root is supplied. Oracle secant seeds are 2.8−0.2i,2.9−0.1i.
- Three finite-tail cutoff controls at d=1.2,q=1.3 and the finest tolerance.

Finest maximum complex S entry error ≤1e-7, real-axis unitarity/reciprocity error
≤1e-8; no unitarity condition is imposed at complex energy. Finest q error and
radius q spread ≤1e-7, position/full-width errors ≤0.1 meV, numerical Jost
determinant residual ≤1e-6. All root statuses must succeed at every rung with
Re q>0,Im q<0, 0<Re ε<16 and positive full width. Incoming-oracle determinant
at each outgoing pole must exceed 1e-3. Finest maximum S/q errors must improve
tenfold from coarse unless already below the 1e-9 roundoff floor; every rung
remains in the report. Both reference precisions must agree within 1e-40.

Mandatory negatives at c=3: selecting a growing rather than decaying closed
channel at ε=8 must change S by >1e-3; omitting relative flux normalization at
ε=16.1 must change S by >1e-3. These wrong real-axis matrices can hide behind
scalar phase/unitarity checks. At d=1.2,q=1.3, omitted-tail error at B/a=8 must
exceed 1e-3, and at 128 be at least fivefold smaller. The cutoff comparison is
reported separately from exact-tail matching accuracy.

All missing/nonfinite results, solver failures and false/unset gates stop the
command after writing available evidence. No retries, reseeding, grid/tolerance
changes, native input reads, exports, inner solves or automatic successors occur.
The previous scalar benchmark and historical native rejections remain separate.

## Interpretation

A pass qualifies these two test models and their stated sheets. It does not
validate a general noncommuting unequal-threshold **and** polar-tail problem,
native amplitude normalization or a CO complex pole. The next implementation
must combine those ingredients with independent controls, derive the native
writer/assembler conventions and reproduce the saved full K matrices before CO
pole searches. Numerical agreement is distinct from model/channel convergence.
