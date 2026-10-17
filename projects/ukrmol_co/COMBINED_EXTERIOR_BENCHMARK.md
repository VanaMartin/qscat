# Combined polar exterior: refinement gate rejects — 17 October 2026

**The prospectively frozen combined boundary-operator gate rejects:**
**15/16 checks pass**, but the spectral middle→finest change is
**7.90991e-7**, above the declared **1e-7** budget. All **432 paired-method
cases and four analytic controls complete**. Original scientific command and
container exits are **1/1**; no grid change, tolerance relaxation or successor
run is used to replace this result.

[Design](../../docs/physics/co-combined-exterior-benchmark.md),
[contract](combined-exterior-contract.json), sources and tests were frozen at
**`a3300b0` before the full matrix**. The independent methods are nonlinear
inward Riccati integration and a second-order linear Chebyshev boundary-value
solve. They share the finite-endpoint Robin condition, so endpoint/angle and
analytic-limit controls are separate requirements. This is a small Python toy
experiment, without native input or a native solve.

## Combined physics and outgoing continuation

The model is `−u_xx+[T+C(d)/x²]u=εu`, x=r/a≥1, a=5 bohr,
`T=diag(0,16)`, `C=[[2,−2d/√5],[−2d/√5,6]]`, d=0/0.4/1.2.
T and C **do not commute** at nonzero d. It retains a genuine inverse-square
transition tail between unequal internal thresholds; a constant basis rotation
cannot reduce it to the earlier degenerate-threshold benchmark.

The observable is the outgoing exterior map **u_x(1)=Y(ε)u(1)**. Six real
energies include **15.9/16.1**, straddling the upper threshold. Two complex
energies **8−0.2i / 18−0.2i** use explicitly continued outgoing momenta.
Exterior rotations **20°/35°** make every outgoing/closed momentum decay.
Endpoint lengths **40/60/80** and three ODE/spectral accuracy levels remain
in the report. A decaying closed channel may have an algebraically forced
component from the open-channel tail; it is not set to zero.

The endpoint condition is the leading outgoing derivative `u_z=iK u`, with
the coupled tail retained everywhere between the inner boundary and endpoint.
Contour/domain stability assesses its endpoint approximation. The returned
derivative includes the `exp(−iθ)` chain-rule conversion. Y is not an S matrix;
its real-axis invariants are transpose symmetry and outgoing flux positivity.

## Numerical outcome and failed refinement diagnosis

| Check | Measured maximum | Frozen budget / verdict |
|---|---:|---|
| Finest Riccati–spectral Y difference | **2.65381e-9** | 1e-7, pass |
| Finest analytic-limit error | **1.59637e-9** | 1e-7, pass |
| ODE middle→finest change | **2.40874e-8** | 1e-7, pass |
| Spectral middle→finest change | **7.90991e-7** | **1e-7, reject** |
| ODE endpoint/angle spread | **1.60161e-9** | 1e-7, pass |
| Spectral endpoint/angle spread | **3.84610e-9** | 1e-7, pass |
| Analytic reference 60→90-digit change | **4.87279e-61** | 1e-40, pass |

Riccati settings rtol **1e-5/1e-8/1e-11** give successive maximum changes
**6.58931e-5 → 2.40874e-8**. Spectral degrees **64/128/256** give
**0.0844746 → 7.90991e-7**. The spectral change improves by far more than
tenfold but still misses the separately required final absolute budget.
The finest reference comparison cannot retrospectively waive that gate.

Exactly **three of 144 physical/domain/angle configurations** exceed the
spectral refinement budget. All have **ε=30,L=80,θ=20°**, one for each dipole.
Their degree-128 method differences are **7.907–7.931e-7**, whereas degree-256
method differences are **8.36e-14–4.08e-11**. The zero-dipole member already
has the defect, anchoring it to the long-domain/high-energy discretization
control rather than requiring noncommuting polar physics to explain it.
This is a source-bound diagnosis, not a repaired verdict or additional run.

The largest collocation infinity-norm condition estimate across all settings
is **2.93239e7**, with maximum linear-system residual **1.016e-9**. Record those
conditioning effects explicitly. Both methods pass finest real-axis transpose
symmetry and flux positivity. Exact uncoupled unequal-threshold and noninteger-
order degenerate-threshold Hankel limits pass. No coupled exact oracle is claimed.

## Negative controls and native normalization progress

At ε=8,d=1.2, dropping the dipole changes Y by **0.180034**; dropping the
rotated-coordinate derivative conversion changes it by **2.14709**. In the
decoupled analytic ε=8 control, reversing the open outgoing sign changes Y by
**5.02831**. All exceed the fixed 0.001 detection threshold.

The independently completed [native convention trace](../../docs/physics/co-native-boundary-conventions.md)
resolves an important factor-of-two ambiguity:

- GBTOlib `UKP_READAMP` **divides surface amplitudes by sqrt(2)** before CI
  contraction. The writer's later sqrt(2) conversion applies only to interior
  RMT sample points, **not the scattering surface**.
- The active native zero-Bloch sum is **Q(E)=Σ wwᵀ/(E_n−E)**, with total
  **Hartree** denominators and **bohr** units, mapping u' to u. Do not insert
  another half or divide by the boundary radius while keeping an unscaled u'.
- Exported multipole coefficients are **already converted to the Rydberg
  radial operator**. Do not double them again. Native requested-eV grids retain
  their historical 0.0735-Ryd/eV conversion.

The [Python sum](native_boundary.py) and two independent analytic tests are
frozen at **`d621541`**. A free Neumann box anchors `Q=tan(ka)/k` and explicitly
tests the finite-pole omitted-tail behavior at **128/512/2048** terms. These
tests do not assert that an uncorrected finite sum is a complete physical R.
Seven unmodified native Fortran sources, their build provenance and three
COPYING files are checksum-preserved. Native full K and general nonzero-Bloch
propagation remain unqualified.

## Reproducibility, resources and next decision

[Twelve public payloads](qualification-evidence/combined-exterior/README.md)
preserve all boundary maps, the original false verdict, both source freezes,
native source/normalization anchors, the exact dependency wheel and resource
records. Host/local/public reconstruction preserves every case and all sixteen
check verdicts, including the scientific exit 1. See
[numerical verification](verification-combined-exterior-20261017.json),
[native source trace](verification-native-boundary-20261017.json) and
[public replay](verification-combined-exterior-public-20261017.json).
Maximum host/public Y-component drift is **4.186e-8** across all rungs and
**3.008e-9** at finest settings, below the 1e-7 matrix budget. Both independently
replayed native analytic normalization tests pass.

Sadaharu cost **11.178 s wall / 0.073 GiB sampled kernel high-water**,
54 samples at 0.2 s, reserved physical cores 0–3, one BLAS thread and 4-GiB
limit. Sampling is not a post-exit peak. The previously pinned mpmath wheel is
reused read-only; its historical scalar-study missing-dependency probe is
explicitly identified rather than asserted as a new failure. Persistent
artifacts: `prepared/combined-exterior-source-20261017/`. Source archive SHA256:
`8a952365752a3b16298cd90a7c2d85752462dde86cb32a42f9e635ebdf9481f2`.
**37 focused tests pass.**

The next decision is a separately justified, prospectively frozen spectral
resolution family addressing the measured ε=30,L=80 discrepancy and its cost.
There is no automatic refinement successor. Only after the combined gate passes
should full saved-native **K**, radius/step/tail and asymptotic normalization
controls be released. CO pole searches remain unset, along with spectral/model,
continuity, extraction and neutral-domain qualification. All earlier failures
and separate scalar/matrix benchmark passes retain their original meanings.
