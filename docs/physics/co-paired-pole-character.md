# Compressed-sector counterpart and open-residue character

## Decision and scope

The [compressed B1 result](../../projects/ukrmol_co/COMPRESSED_POLE_REFERENCE.md)
is one visible pole of a finite CAS11/nine-A2 reference, not yet a Pi assignment.
Before extending its geometry curve, ask whether its **open-channel residue**
has the intended angular character and whether the saved B2 model has a matching
counterpart. Degenerate B1/B2 poles can also have higher odd angular projection.
Do not use degeneracy, a phase-fit width, or an inner R-matrix eigenvalue alone.

This is design → small project Python classifier → synthetic controls → one
staged application. No solver port, Rust optimization or qscat promotion.
The [contract](../../projects/ukrmol_co/paired-pole-character-contract-20261023.json)
freezes every threshold **before any actual CO residue-character score**.
First screen saved B1 evidence without a propagation. Only a B1 pass makes the
bounded B2 calculation useful; rejection/inconclusive stops before any B2 map.
No new sector/anchor beyond this counterpart, geometry, electronic calculation,
root/channel ladder, potential fit or factory export belongs to this experiment.

## Native labels, not array-index assumptions

The native channel record is `(channel, target, l, m)` followed by its threshold.
Authenticated `Postprocessing_module.F90::get_channel_info`, lines765–778,
obtains `lchl`, `mchl` and target indices; `::write_channel_info`, lines1082–1083,
writes those exact columns. Authenticated `outerio.f::READTC`, lines1235–1239
and1304–1305, names/reads orbital angular momentum LCHL and its M quantum number.
Both source hashes are pinned in the contract; these are native-source locators,
not published-paper citations. Existing saved-matrix metadata independently
records the same labels.

CO lies along z; the neutral ground target is the existing spin-singlet A1
ground-state assignment (native target label `[1,0,1,0]`). All six open channels
attach to target1. For that Sigma target, projectile `|m|=1` identifies the Pi
open-channel subspace; `|m|=3` identifies higher angular character in these sectors.
This use of projectile m would **not** suffice for an excited target with its
own nonzero angular projection. Reject incompatible metadata, not reinterpret it.

Verified open labels:

```
B1: (1,1,1,1), (2,1,2,1), (3,1,3,1), (4,1,3,3), (5,1,4,1), (6,1,4,3)
B2: (1,1,1,-1),(2,1,2,-1),(3,1,3,-3),(4,1,3,-1),(5,1,4,-3),(6,1,4,-1)
```

Thus Pi indices are B1 `[0,1,2,4]`, B2 `[0,1,3,5]`. Pair labels by
`(target,l,|m|)`, never by channel ordinal or the sign of m. Sum across **all l**;
neither p-wave dominance nor single-l sufficiency follows from a Pi score.

## Phase-free necessary character test

Let R be the **S-matrix pole residue** in the saved flux-normalized open basis,
in Hartree, and P the orthogonal diagonal projector onto target1, `|m|=1`.
With the ordinary Frobenius diagnostic norm define

    eta_out = ||(I-P) R||_F^2 / ||R||_F^2,
    eta_in  = ||R (I-P)||_F^2 / ||R||_F^2,
    eta = max(eta_out, eta_in).

No eigenvector gauge, Takagi factor or arbitrary inter-sector phase alignment
is needed. R itself still comes from the **bilinear**, non-conjugated pole
construction; the conjugated Frobenius norm is only a diagnostic of its support.
For a rank-one reciprocal residue `R=c g g.T`, both fractions reduce to the
Euclidean support fraction of g outside P. These are **not probabilities,
partial widths, branching fractions or an interior-state Lambda expectation**.
They can miss closed-channel structure; even zero leakage does not prove exact
Lambda conservation of the finite Gaussian/C2v electronic model.

Validity guards: finite nonzero square R with compatible verified labels,
positive full width, `||R||_F/Gamma > 1e-6`, singular-value ratio
`sigma2/sigma1 <= 1e-3` and relative transpose defect `<=1e-3`. A failed signal,
rank or symmetry guard is inconclusive, not a fabricated angular rejection.
Malformed/unverified labels are an input failure, also qualification unset.

Prospective per-sector policy:

|Measured leakage eta|Disposition|
|---|---|
|`eta <= 1e-4`|pass open-residue Pi-character condition|
|`eta >= 1e-2`|reject this pole's proposed Pi open-residue assignment|
|between|inconclusive; no tolerance change or automatic refinement|

These are declared decision tolerances, not universal symmetry constants or
electronic-error bounds. The pass permits at most1% leakage **amplitude norm**;
the reject threshold corresponds to10%. The two-order gap avoids a forced binary
classification. Saved residue-construction errors near1e-7 relative are much
smaller than these policy scales, but that does not qualify electronic accuracy.

Use all four existing estimates:16-node residue,32-node residue, and both fixed
derivative residues. Require all four dispositions identical, maximum absolute
leakage change `<=1e-5`, canonical profile total-variation change `<=1e-4`, and
relative norm change `<=1e-3`. Otherwise inconclusive. Never infer zero error
from identical root coordinates or rename numerical stability model convergence.

The diagnostic is invariant under nonzero complex residue scaling, channel
permutations with labels, channel phase/sign changes, and orthonormal changes
within the Pi/complement subspaces. Mixing those subspaces without transforming
P is a physical change, not a gauge invariance to demand.

## What inter-sector comparison can establish

After both sector screens pass and both poles have the frozen contour/residue
witnesses, compare position and **full width** at `1e-4 eV` (0.1meV) each,
using canonical `qscat.units`. Compare phase-free canonical channel profiles

    p_i = (sum_j |R_ij|^2 + sum_j |R_ji|^2) / (2 ||R||_F^2)

after matching `(target,l,|m|)`. Use total variation `0.5 sum_i |p_i(B1)-p_i(B2)|`
and symmetric norm mismatch `|N1-N2|/max(N1,N2)`: each passes at `<=1e-3`, rejects
at `>=1e-2`, otherwise inconclusive. These profiles are invariant to phases and
reordering, **not** arbitrary mixing of different l labels. Report every l and
|m| contribution; do not introduce a post-hoc dominant-l gate.

Profile/norm agreement discards relative channel phases and is only a necessary
counterpart consistency check, not equivalence of residues/full S, an exact
axial-symmetry proof or geometry continuity. A pass supports the proposed
finite-reference open-Pi pair. It returns priority to a decision-changing
**matched electronic-model pole contrast**, not another general exterior study.
It does not select a production model or qualify a CO curve/potential.

## Minimal B2 application, only after the saved B1 screen passes

Reuse every B1 map/witness; do not rerun its pole. B2 retains all298 channels /
2048 roots and the same native exterior. Before complex work, check B2 full S at
development indices60 and80, using only final setting120bohr /40deg /order10 /
rtol1e-12. Require entry, unitarity, reciprocity and Wronskian errors `<=1e-6`,
flux defect `<=1e-7`, and eigenphase-sum error modulo pi `<=1e-6 rad`.
These are minimum B2 normalization/consumer checks, not another convergence
campaign. Neither index is a route holdout; no reserved S value is analyzed.

Use the already known B1 pole as B2 seed; the second seed is +2e-5Ha.
The fixed B2 box is +/-0.002Ha in each coordinate about that B1 pole. A box
exit/cap or no candidate is inconclusive for this local counterpart search,
**not proof that B2 has no pole**. No alternate seed or expanded region.
Use the original three map settings and root caps12/8/8. Reuse the existing
bordered search, projected interior contour and derivative/open-residue routines;
do not build another solver. Retain the fixed radius0.001Ha,32 contour points
(nested16), eight probes /seed314159, derivative steps1e-5 and5e-6Ha, winding+1,
and every original pole witness gate. Merely matching a boundary zero is not a pass.

At most **67 new maps**:2 real +28 root +1 center +4 derivative +32 contour.
Eight physical cores0–7 /one BLAS thread per worker /at most8 contour workers /
8GiB. Target600s /upper1200s /hard cap1800s;90s/map. The previous47 maps cost181s,
with1.94GiB cgroup peak; this motivates allocation, not guaranteed completion.
Saved B1 screen needs no map, four physical cores,1GiB and60s hard cap.

Implementation/diagnosis is limited to this active working cycle, at most two
hours additional effort after the native-label reconnaissance. The deliverable
is this frozen test and one staged result, not a general assignment package.
If the B2 adapter/source cannot be frozen within that bound, stop at the saved
B1 result and hand off the exact missing execution step. No automatic repair,
successor, extra quadrature, root/channel/precision ladder, geometry or fit.

## Controls and execution boundary

Before actual scores: pure Pi at l=1 **and** l=3; pure higher-|m| doublet with
an identical pole (must reject despite degeneracy); mixed-grey support; rank-two,
zero-signal and asymmetric residues; wrong/missing/duplicate labels; permutation,
complex phase/scaling and block-unitary controls; profile-mismatched counterpart;
and perturbation across a decision boundary (must be inconclusive).

Freeze/push design, classifier, controls and contract before the saved B1 score.
The saved-only CLI cannot propagate or search. B2 requires its separately frozen
adapter and exact owner/runtime before use; the current B1 source freeze alone
never authorizes it. Preserve original exits, every partial artifact and source/
input hashes. Exit0 means the stage's stated finite condition passes, exit1 a
resolved rejection, exit2 inconclusive/input/runtime stop. Report the paired
verdict as unset until **both** sectors and all counterpart gates pass.

Protect all114 route holdout S values and the equilibrium pole-target exclusion.
New50 stays parked. Electronic accuracy, omitted spectrum, independent exterior,
geometry continuity, neutral correlation, training-node count, asymptote and
full-scattering embedding remain separate, unresolved production gates.
