# Noncommuting thresholds and polar tails: exterior boundary operator

## Model and decision

The separate [tail and threshold controls](co-coupled-tail-benchmark.md) pass.
Their combination needs an independent check before native CO propagation.
Test the **exterior Dirichlet-to-Neumann map**, not fitted phases: for prescribed
two-component u(1), return its outgoing derivative u_x(1)=Y(ε)u(1). This is the
boundary operator needed for interior matching and closed-channel elimination.
It has no incoming-wave normalization ambiguity. Scattering amplitudes still
require a subsequent full K-matrix reconstruction gate.

Atomic units, a=5 bohr, x=r/a and ε=2a²E:

\[
 -u_{xx}+[T+C(d)/x^2]u=\epsilon u,\qquad
 T=\operatorname{diag}(0,16),\quad
 C(d)=\begin{pmatrix}2&-2d/\sqrt5\\-2d/\sqrt5&6\end{pmatrix}.
\]

The channels represent two internal thresholds carrying l=1,2, |m|=1 angular
waves with a dipole transition matrix element. This defines a self-adjoint toy
Hamiltonian on the real half-line. For d≠0, **[T,C]≠0**: a constant angular
rotation cannot decouple the radial equation. The tail is uncut in the model.
There is no claim that these two target channels describe CO convergence.

## Outgoing sheet and exterior rotation

Let k0=+sqrt(ε), and below the second threshold use
`k1=+i sqrt(16−ε)`; above it use `+sqrt(ε−16)`. Complex points are continued
from their declared side of the threshold, away from the branch point. At the
two lower-half-plane energies an open outgoing wave grows on the real axis.
Rotate the exterior coordinate from its boundary:

`z(t)=1+(t−1)exp(iθ)`, 1≤t≤L. Require **Im(k_i exp(iθ))>0** for every
asymptotic outgoing/decaying channel. This exposes the outgoing solution as a
decaying BVP. The potential is analytic along these contours; the only inverse-
square singularity z=0 lies outside them. Derivatives returned at t=1 are
converted to **d/dz**, not left as d/dt.

At the finite numerical endpoint impose `u_z=iK u`. This omits asymptotic
inverse-power corrections **at the endpoint**, not the tail throughout the
propagation region. Its effect at z=1 must converge separately under L and θ;
two methods sharing that condition alone do not establish the infinite-domain
answer. Endpoint lengths 40/60/80 and angles 20°/35° are frozen. Threshold-near
points ε=15.9/16.1 are included; exact threshold is excluded.

## Independent numerical paths

1. **Riccati ODE:** set Y(L)=iK and integrate inward
   `dY/dt=exp(iθ)[T+C/z²−εI−Y²]` with DOP853 at rtol 1e-5/1e-8/1e-11,
   atol=0.01rtol. This is a nonlinear 2×2 boundary-operator propagation.
2. **Spectral BVP:** independently discretize the second-order linear equation
   on Chebyshev–Lobatto nodes, with u(1)=I and the same endpoint outgoing Robin
   condition. Solve the full two-channel collocation matrix at N=64/128/256.
   Return `exp(−iθ)D_t u(1)`. This route does not integrate a Riccati equation,
   use ODE output as a boundary condition, or share its radial grid.

Neither method is an exact coupled oracle. Their differential agreement is
combined with three-rung discretization, separate domain/contour controls,
analytic limits and real-axis current checks. Record every linear condition
number, collocation residual, ODE status and full complex Y, including failures.

For **d=0**, the unequal thresholds decouple and the independent analytic
answer is `Y_ii=1/2+k_i H'_(l_i+1/2)(k_i)/H_(l_i+1/2)(k_i)`, with zero
offdiagonals. Evaluate with mpmath at 60/90 digits. For **T=0**, use the exact
noninteger-order Hankel eigensystem from the earlier tail benchmark, checking
the new operator against H+'(1)H+(1)⁻¹. These limits jointly anchor radial
phases, derivative normalization, and both threshold and angular ingredients.

## Frozen matrix and gates

The JSON contract and sources are frozen before the full matrix:

- d=0/0.4/1.2; ε=0.5/8/15.9/16.1/22/30 and 8−0.2i/18−0.2i;
  three endpoints, two angles and three paired accuracy levels:
  **432 coupled-domain cases**, each evaluated by both methods.
- Degenerate-threshold controls: d=0.4/1.2, ε=1/4−0.2i,
  L=80, θ=35°, finest settings: **four additional cases**.
- Finest ODE/BVP Y-entry difference ≤1e-7 across all domains/angles.
  Finest analytic-limit errors ≤1e-7; reference precision change ≤1e-40.
  Each method's finest-versus-middle change ≤1e-7; its maximum coarse-to-middle
  change must shrink tenfold at finest, with a 1e-9 roundoff floor.
- Each method's finest endpoint/contour spread ≤1e-7. Report L=40/60/80
  differences independently; do not call this a rigorous absolute tail bound.
- Real-axis transpose symmetry ≤1e-7; minimum eigenvalue of
  `(Y−Y†)/(2i)` ≥−1e-7. This is outgoing flux positivity, not unitarity of Y.
  Closed-channel algebraic components induced by the open tail are retained.
- All ODE solves succeed, arrays are finite, and all rotated outgoing momenta
  decay. A Riccati singularity is an explicit failed case, without chart switching.
- Wrong open-channel incoming sign in the decoupled exact ε=8 control must
  change Y by >0.001; omitting dipole coupling at d=1.2,ε=8 must do likewise.
  Omitting the `exp(−iθ)` derivative conversion must also differ by >0.001.

False/unset checks stop after saving evidence. No retries, changes to the
declared matrix, fitted poles, native solves or automatic successors. Limit the
Sadaharu container to reserved physical cores 0–3, one BLAS thread, 4 GiB and
300 seconds. Retain the original exit and measured resource samples.

## Interpretation

A pass qualifies this two-channel outgoing **boundary operator** in the tested
real/complex regime, including simultaneous noncommuting thresholds and an
inverse-square tail. It does not qualify S-amplitude normalization, poles of an
interior-plus-exterior system, the CO channel/model truncation, quadrupole or
polarization completeness. Native surface/Bloch/energy conventions require
their own source trace; full native open-open K reproduction remains the next
observable gate before a separately frozen CO pole search.
