# FEM-DVR-ECS: finite-element DVR with exterior complex scaling

**Location:** `qscat.dvr` (grid + kinetic + Hamiltonian/eigen helpers),
`qscat.ecs` (the complex-scaling coordinate map).
**Origin:** ported from eMoScat's `FemDvrEcsGrid.cpp` / `DvrGrid.cpp` /
`KineticEnergy.cpp` / `DiscreteStates.cpp`. The method is Rescigno & McCurdy,
Phys. Rev. A 62, 032706 (2000) -- see
`reference/literature/rescigno-2000-pra62-032706.md`; the design rationale is
recorded at `docs/superpowers/specs/2026-07-21-femdvr-ecs-grid-design.md`.
**Units:** atomic units throughout (energy in Hartree, length in Bohr).

## Method

A Finite-Element Discrete Variable Representation (FEM-DVR) discretizes a
radial coordinate into a sequence of elements, each carrying its own
Gauss-Lobatto-Legendre (GLL) quadrature nodes and Lagrange interpolating
basis functions. Adjacent elements share their boundary node ("bridge"
function), and the two outermost grid points are dropped to enforce a
Dirichlet ($\psi = 0$) boundary condition at both ends. On the
shared-quadrature GLL basis, the kinetic and (diagonal) potential operators
are simple to assemble and the resulting matrix eigenproblem gives bound and
(with ECS) resonance/continuum states directly, without an explicit
numerical integration of the Schrodinger equation.

Exterior Complex Scaling (ECS) extends this to scattering problems: beyond a
pivot radius $R_0$ (chosen to sit exactly on an element boundary), the radial
coordinate is rotated into the complex plane by a fixed angle $\theta$:

$$
z(x) = \begin{cases}
x & x \le R_0 \\
R_0 + (x - R_0)\,e^{i\theta} & x > R_0
\end{cases}
$$

This turns the oscillatory, non-normalizable continuum wavefunctions of the
unscaled problem into square-integrable, decaying functions, which:

- exposes resonances as isolated, $\theta$-stationary complex poles of the
  discretized Hamiltonian, and
- rotates the true continuum spectrum by $\arg(E) \sim -2\theta$ (an
  asymptotic statement, sharper for a pivot far out on the tail — see
  Benchmark 3),

while leaving true bound-state energies (which lie under the rotated
continuum) numerically unchanged for any $\theta$ in a "stable" window.

## Grid construction (`qscat.dvr.grid.FemDvrEcsGrid`)

Given a `GridSpec` (quadrature order `nq` shared by all elements, an ordered
list of `ElementSpec(length, angle_deg)`, and an inner boundary `x_min`):

- Real (unscaled) cumulative element boundaries `ar[i]` are plain running
  sums of element lengths; `R0 = x_min + sum(real element lengths)` is
  computed once in `GridSpec.__post_init__` and is guaranteed to land on an
  element boundary by construction.
- Each element's node positions are built from the reference GLL nodes
  $\xi$ on $(-1, 1)$ (`gll.gll_nodes_weights`), first placed in real
  (unscaled) space, then passed through `qscat.ecs.ecs_map(x, R0, angle_deg)`
  to get the actual (possibly complex) grid point. This is the ECS map's
  only use in the grid: **`ecs_map` is the single source of the coordinate
  transform**, so there is exactly one place in the codebase that defines
  "what counts as an ECS-scaled point."
- The quadrature *weight* Jacobian is independent of how the point is
  computed — it is `hz = 0.5 * length * exp(i*angle_deg)`, the derivative of
  the linear map restricted to that element, and is bridge-summed
  (`+=`) at shared element boundaries.
- The two outermost global points (`x_min`, `x_max`) are dropped (Dirichlet),
  leaving `n = tnel*(nq-1) + 1 - 2` basis functions.
- `element_maps` gives each element's explicit `(local_idx, global_idx)`
  scatter map into the global basis (needed because the Dirichlet drop makes
  the local-to-global offset non-uniform at the first/last element).

**Caveat:** a single ECS angle shared by every element in the complex tail
is the validated configuration (used by all benchmarks below and by
`ecs_map`, which takes one `theta_deg`). A "bent" contour with different
angles on different tail elements is unverified/experimental — `GridSpec`
does not forbid it, but no benchmark exercises it.

## Kinetic-energy assembly (`qscat.dvr.kinetic.kinetic`)

$T = -\dfrac{1}{2\,\mathrm{mass}}\dfrac{d^2}{dz^2}$ is assembled element-by-element:

```text
wze[l]    = hz * wl[l]                      # scaled GLL quadrature weight
dBF[a, l] = dLp[l, a] / hz                  # scaled Lagrange derivative
dBF[a, :] /= sqrt(weights[global_idx(a)])   # normalize by the GLOBAL
                                             # bridge-summed weight, not the
                                             # local element weight
T_local[a, b] = (1/(2*mass)) * sum_l wze[l] * dBF[a, l] * dBF[b, l]
```

`T_local` is computed over all `nq` local indices, then the retained
sub-block (per `element_maps[k]`) is scatter-added into the global matrix;
adjacent elements share one bridge index, so `+=` assembles the
bridge-corner coupling automatically. The classic assembly bug here is
normalizing by the wrong (local vs. global) weight — this is why Benchmark 1
(particle-in-box) is an exact oracle: it is sensitive enough to catch it.

## Diagonal-potential approximation (`qscat.dvr.operators.hamiltonian`)

$H = T + \operatorname{diag}\!\big(V(\text{points})\big)$. Because the DVR
basis functions are (by construction) orthonormal under the same Lobatto
quadrature that builds `T`, the potential-energy matrix element $V_{ij}$ is
well approximated by $V(x_i)\,\delta_{ij}$ — no explicit off-diagonal
quadrature is needed. This is an approximation (not exact quadrature of
$\langle i|V|j\rangle$), and it degrades if $V$ has structure the grid can't
resolve — in particular, a potential discontinuity should land exactly on an
element boundary (see Benchmark 4) rather than inside an element, or the
diagonal approximation smears it.

$H$ is complex-symmetric but non-Hermitian in general (once any element is
ECS-rotated), so `qscat.dvr.operators.eigen` uses the general eigensolver
(`np.linalg.eig`, corresponding to LAPACK `zgeev`) rather than a Hermitian
one, and sorts results by ascending $\mathrm{Re}(E)$.

## Validation benchmarks

All four live in `libs/qscat/tests/test_femdvr_ecs.py`.

1. **B1 — particle in a box** ($\theta = 0$). Exact analytic oracle
   $E_n = n^2\pi^2 / (2\,\mathrm{mass}\,L^2)$. Matched to `rtol <= 1e-6` for
   the first five levels, imaginary parts `~0` to `1e-9`. A companion
   spectral-convergence check (`nq = 4, 5, 6`) confirms the ground-state
   error falls monotonically (exponentially) as quadrature order increases —
   the signature of a correct spectral (FEM-DVR) discretization. This
   benchmark is the primary arbiter of kinetic-assembly correctness because
   it is maximally sensitive to bridge-weight normalization, the Dirichlet
   trim, and scatter bookkeeping errors.
2. **B2 — harmonic oscillator** ($\theta = 0$). Analytic oracle
   $E_n = \omega\,(n + 1/2)$. Matched to `rtol <= 1e-6` for the first five
   levels. Exercises the diagonal-potential DVR approximation on a smooth,
   everywhere-differentiable potential.
3. **B3 — ECS continuum rotation**. For a free particle ($V = 0$) on a grid
   with a real region of length $R_0$ followed by a single complex tail of
   length $L_t$ at angle $\theta$, matching $\psi$ and $\psi'$ at $z = R_0$
   gives the *exact* (not asymptotic) quantization condition
   $k\,(R_0 + L_t e^{i\theta}) = n\pi$, i.e.
   $E_n = n^2\pi^2 / (2\,\mathrm{mass}\,Z_\mathrm{eff}^2)$ with
   $Z_\mathrm{eff} = R_0 + L_t e^{i\theta}$; $\arg(E_n) = -2\arg(Z_\mathrm{eff})$
   for *every* $n$. The textbook "continuum rotates by $-2\theta$" picture is
   the $R_0 \ll L_t$ limit of this ($\arg(Z_\mathrm{eff}) \to \theta$). The
   benchmark grid uses $R_0 / L_t = 0.05$, so $\arg(Z_\mathrm{eff}) = 28.63°$ (vs. the
   $\theta = 30°$ asymptote), giving $\arg(E) = -57.25°$ (vs. the
   $-60°$ asymptote) — well inside a $\pm 5°$ window; the test asserts
   that most mid-spectrum, sizeable-$|E|$ eigenvalues cluster there
   (edge-of-basis "junk" states and near-zero eigenvalues are excluded from
   the selection).

   **Why the grid is deliberately lopsided.** An *equal* real/complex split
   ($R_0 = L_t$) gives $\arg(Z_\mathrm{eff}) = \theta/2$ exactly — since
   $1 + e^{i\theta} = 2\cos(\theta/2)\,e^{i\theta/2}$, independent of the
   common length — so the spectrum clusters at $\arg(E) = -\theta$, not
   $-2\theta$. A benchmark written that way fails its own assertion with
   *zero* eigenvalues in the expected window, and looks like an ECS bug
   while being a test-setup artifact. That it is an artifact was confirmed
   separately: a *uniform*-angle grid (all elements rotated, no real region)
   reproduces the exact scaled-box spectrum
   $E_n = n^2\pi^2 e^{-2i\theta} / (2mL^2)$ to machine precision, so the
   $e^{i\theta}$ Jacobian in `kinetic` is right. Hence `R0/Lt = 0.05`.
4. **B4 — bound-state $\theta$-independence**. A square well $V = -V_0$ on
   $[0, a]$, deep enough to support one bound state. The well edge `a = 3.0`
   is placed exactly on an element boundary (real region length 12, 4
   elements of length 3 each) so the diagonal-potential DVR represents the
   discontinuity cleanly. The bound-state energy is computed at two
   different ECS angles (`theta = 20 deg` and `35 deg`); physically a true
   bound state (lying below the rotated continuum) must not depend on
   $\theta$, and the two energies are required to agree to `< 1e-4`.

## Known limitations / out of scope

- Only a single, contiguous ECS tail at one shared angle is validated;
  `complex_negative` (ECS at the *inner* boundary) and graded/bent tail
  angles from the original eMoScat input format are not implemented.
  (`GridSpec.__post_init__` enforces contiguity of complex elements but not
  angle uniformity — `ecs_map` needs one `theta_deg` when it is used, so a
  per-element angle in `qscat.dvr.grid` is only exercised with a uniform
  tail in tests.)
- The diagonal-potential approximation assumes the potential is smooth
  within each element (or that any discontinuity is placed on an element
  boundary, as in Benchmark 4).
- Resonance identification via a $\theta$-stabilization scan (finding the
  angle range where a complex eigenvalue is stationary) is not automated —
  eMoScat only automated bound-state stability checks (`prec 1e-4`, mirrored
  by Benchmark 4).
