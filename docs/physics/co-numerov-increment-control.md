# Fixed-grid Numerov floating-point control

## Question and scope

The [two-point attribution](../../projects/ukrmol_co/NATIVE_K_ATTRIBUTION.md)
rejects absolute K gates when halving the step: the finest/native64 discrepancy
is about 7.26e-5 while tail-order sensitivity is about 4.4e-11. This motivates
an algebraic precision control at **the existing grids**, before another
physical or discretisation refinement. It does not establish roundoff as the
cause. All original failures and producer bytes remain evidence.

This experiment implements a separate Python toy function
`increment_numerov_k(q, thresholds, angular, c2, c3, energy, radius, endpoint,
step, order, degeneracy=1e-10)`. It shares the existing tail construction and
unsymmetrized matching, but avoids subtracting nearly equal identity matrices
in initialisation, the ratio recurrence and final boundary matching. This is
one declared floating-point formulation, not a new potential or fitted boundary.
The original `exterior_numerov.numerov_k` is unchanged.

## Algebra, including matrix order

The external equation is `u'' = W(r) u`, with native Rydberg coefficients
`W = diag(thresholds-energy) + c2/r² + c3/r³`, lengths in bohr. Native incident
energy is converted to Hartree only when evaluating the inner surface map Q.
Let `B_n = I-h² W_n/12`, `z_n=B_n u_n` and
`R_n=z_(n-1) z_n⁻¹`. The existing Numerov recurrence is

```
R_(n+1) = [12 B_n⁻¹ - 10 I - R_n]⁻¹.
```

Since `12 B⁻¹ - 10 I = 2 I + h² B⁻¹ W`, define
`D_n = I-R_n`, `C_n = h² solve(B_n,W_n)`. Then

```
X_n = C_n + D_n
D_(n+1) = solve(I + X_n, X_n).
```

The last identity follows from `I-(I+X)⁻¹=(I+X)⁻¹X`; it never moves unrelated
matrices through each other. Q need not be invertible. With `u(a)=Q,u'(a)=I`,
compute the sixth-order Taylor **increment** `delta_u=u(a+h)-Q` directly from
terms of degrees 1–6. Then

```
D_1 = [(B_1-B_0) Q + B_1 delta_u] [B_1 (Q+delta_u)]⁻¹.
```

At the endpoint, let `u(b-h)=T_v u(b)+T_d u'(b)`, using the same original
sixth-order Taylor maps and analytic derivatives of W. The value and derivative
coefficients in the matching equation are

```
A = B_(N-1) (T_v-I) + (B_(N-1)-B_N) + D_N B_N
H = B_(N-1) T_d.
```

`T_v-I` is summed without its identity term; B differences are evaluated as
`-h² (W_left-W_right)/12`, rather than by subtracting two B matrices. The
recurrence, truncation, physical boundary and finite-order tail are unchanged
in exact arithmetic. This combined initialisation/propagation/matching control
cannot by itself localise which floating-point subtraction dominates.

## Analytic qualification before native-data use

1. Check the matrix-order identities on noncommuting symmetric W and boundary
   ratios, and verify singular Q remains supported.
2. Compare short smooth coupled open/open and open/closed problems with an
   independent DOP853 fundamental-solution reference at h=0.2/0.1/0.05. Require
   error reduction at least eightfold per halving and final K error <1e-5.
3. On fixed long free intervals a=18,b=400 at E=0.0735 Ryd, test a scalar
   near-large-K boundary and a rotated two-open/one-closed boundary. Use the
   existing h=0.025/0.0125/0.00625 values. The continuous oracle is direct
   **80-digit** boundary matching of analytic free columns; compare also at
   110 digits (oracle K drift <1e-20). Require candidate finest K error <=1e-6
   and S/phase/invariants <=1e-6. Preserve all original-form errors, even if
   candidate qualification fails. Do not require the original method to fail.
4. Independently evaluate the **discrete scalar Numerov solution** at 80 digits
   through the constant-coefficient trigonometric transfer formula, including
   the same truncated Taylor maps. This distinguishes floating-point drift
   from the discretisation error of the chosen recurrence. Require candidate
   error against this discrete oracle <=1e-7 at each frozen grid.

These are prospective criteria. If qualification is false/unset, stop before
CO data. Local tests are laptop-runnable; measured long-grid evidence is run
on Sadaharu with a fixed 600-s cap, four physical cores and 4 GiB, no engine
invocations. Any successful analytic result qualifies only its tested regime.

The first analytic container exited1/1 before numerical execution because the
pinned tools image lacks mpmath. Preserve that attempt. A fresh dependency-only
successor extracts the [checksum-pinned pure-Python wheel](../../projects/ukrmol_co/numerov-increment-dependencies.json)
to a separate import directory, retaining exactly the numerical source frozen at
`381a201`. No analytic budget or numerical input changes with this supplement.

## Conditional CO experiment

Only after analytic qualification, freeze source, a separate input-hashed
contract and cost estimate before submission. The possible scope is the same
two stretched 1-eV identities, candidate at the three existing spacings with
endpoint400/order8 plus the two existing endpoints200/800 at step0.0125/order8:
**ten independent solves**, no new native invocation or inner solve. Compare
against saved original-form rows and native64, preserving separate 1e-6
K/S/phase/invariant vetoes and all one-factor failures. No choice of a favourable
row or relaxed K budget is allowed. No further rungs or automatic successor.
