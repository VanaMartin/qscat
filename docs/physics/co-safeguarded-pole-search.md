# Safeguarded local boundary-root search

## Scope and lifecycle

This is a project Python safeguard for the existing finite-boundary pole
construction, not a global resonance finder or a new electronic model. It
addresses unrestricted local steps without replacing the outgoing map,
Hamiltonian, energy sheet, native roots, contour/residue construction or
physical-character checks. A finite search region remains a hard allocation,
not evidence that a pole exists inside it. No automatic box expansion, alternate
seed, precision ladder or new geometry follows from a failure.

Lifecycle: mathematical design, NumPy/SciPy toy controls, bounded application
only after runtime controls and a separate immutable release. There is no
measured optimization need and no Rust or qscat promotion.

## Unchanged analytic root function and units

With incident energy z in Hartree, obtain Q(E_ground+z) in bohr and the outgoing
Dirichlet-to-Neumann Y(2z) in inverse bohr. Form M=I-QY as in the
[finite-pole construction](co-compressed-pole-reference.md). At the initial
point take the last left singular vector b and last row c of V^H; freeze both
complex vectors for the entire solve. Define alpha(z) by

    [ M(z)  b ] [ x     ] = [ 0 ],
    [ c     0 ] [ alpha ]   [ 1 ].

The border is constant, so alpha is meromorphic wherever the map is analytic.
No energy-dependent Euclidean normalization, conjugated-energy operation,
singular-value minimization or finite-data surrogate becomes the root function.
Singular values remain diagnostics and acceptance tests, not an analytic scalar.

`safeguarded_search(evaluate, seeds, box, limit, *, max_step=0.025,
max_backtracks=8, armijo=1e-4, observer=None)` takes the existing Q/Y callback.
The two finite seeds specify an initial center and a positive real finite-
difference increment h. The callback is never invoked outside the original
closed box. The per-solve cap counts every unique callback point, including
derivative probes and rejected trials; a cache prevents exact repeats.

## Fixed, bounded iteration

At each current center z, evaluate alpha and the original null diagnostics.
Probe at z+h; if that would leave the box, use z-h. If neither is inside, stop
without a callback. Recompute the forward/backward derivative from this local
probe at every accepted center; a previous long secant chord is not retained.
Zero/nonfinite derivatives stop rather than trigger an alternate search.

The full Newton proposal is p=-alpha/alpha'. Clip its norm to the prescribed
maximum step, with clipping fraction k<=1. Try the fractions t=1,1/2,...,2^-B
of that clipped proposal, where B is the frozen maximum-backtrack count. Reject
an out-of-box trial algebraically without evaluating its map. An evaluated trial
must satisfy

    |alpha(z+t*k*p)| <= (1-armijo*t*k) |alpha(z)|.

This is a sufficient-decrease test for the analytic scalar, not proof of a
physical resonance or a global convergence theorem. Record every proposal,
out-of-box rejection, evaluated value, decrease rejection and accepted center.
If no prescribed trial passes, stop as inconclusive. Do not add probes, enlarge
the region, reset the border or silently accept a worsening step. An exception
from the model/ODE callback is a hard stop, not a reason to try smaller steps.

A finite probe/trial can return a candidate only if the original simultaneous
tests pass: |alpha|<=1e-10, sigma_min(M)<=1e-8, and the second-smallest singular
value>=1e-4 when present. A small alpha with a multiple null, coordinate scaling
artifact or no null is not success. An exact accepted initial point needs no
derivative or trial. No step-size/stagnation tolerance substitutes for a null.

## What this does and does not establish

Fresh derivatives and rejected worsening steps address a specific local-solver
weakness. They do not guarantee convergence from an arbitrary interpolation
seed, locate every pole in a box, distinguish all basins or validate a complex
exterior map. A monotone scalar can still approach a non-root limit or box edge;
the evaluation/backtracking bounds then stop explicitly with no pole accepted.

Every candidate still needs the original three-setting energy stability,
interior resolvent contour, isolated/simple winding, open-channel signal,
derivative/contour residue and angular-character/pair gates. The same-model
position/full-width budget remains0.1meV. Native inputs and all protected
validation exclusions are unchanged. Passing analytic controls authorizes no
CO calculation by itself.

## Required controls

- The exact positive-weight one-level model and independent contour/residue
  tests, including a nonlocal seed and decreasing finite-difference increments.
- A holomorphic polynomial whose full Newton step worsens the residual, proving
  that backtracking occurs and accepted centers decrease rather than merely
  reporting a final root.
- Roots near a box edge, ensuring no out-of-box callback, and a root outside
  the box or a constant/non-root function, requiring bounded failure.
- Multiple/cancelling or invisible features: a candidate/null alone cannot
  bypass the existing simple-root and physical-residue witnesses.
- Callback errors, nonfinite matrices, invalid controls, insufficient local
  probe room, singular borders and exhausted unique-evaluation caps.
- Historical local search/default dispatch remains unchanged unless the new
  search is selected explicitly in a separately frozen application contract.

The runner selects it only with `root_search: safeguarded-newton-v1` and complete
`safeguard_controls` containing `max_step`, `max_backtracks` and `armijo`. Unknown
selectors, missing controls or controls on the legacy default reject. Absence of
the selector preserves the original frozen secant path. The new synthetic
control module also runs in the frozen runtime before any CO root use.

## Exact saved-map reuse

A new application may supply explicit `reused_evaluations` for already captured
root maps. Reuse requires the exact same complex energy, numerical setting,
amplitude mode, sector/case, sheet reference, image, pinned four-native-file
inventory, outgoing-map and surface-map
source bytes, total energy and complete finite298-channel arrays. The capture's
own byte digest must be an input pin. No nearest-energy or cross-frame reuse is
allowed; mismatched declared metadata rejects instead of silently recomputing.
Reused maps remain part of the evaluation cap, are copied byte-for-byte and are
marked `new_ode_call: false`; their original diagnostics describe the earlier
computation, not the successor's elapsed cost. Fresh calls are marked true.
Contour/amplitude maps are outside this narrow reuse path. Actual application
image/runtime identity and the accepted parent evidence must also be frozen.
