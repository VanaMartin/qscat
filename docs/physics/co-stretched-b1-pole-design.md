# First finite stretched B1 pole: own-DZ design

## Decision served and evidence boundary

The finite CAS12 DZ/TZ comparison at R=2.5bohr is intended to test basis-plus-
optimized-frame sensitivity before any pole-targeted potential is fitted. Both
own64-target/400-channel/B1/2048 exports qualify. The separately authorized
20degree real-axis studies now qualify both models' fixed development60/80
identities, including full flux/Wronskian and DZ80/TZ80 endpoint stability.
[Independent completion evidence](../../projects/ukrmol_co/qualification-evidence/stretched-real-completion/README.md)
retains the failed40degree guard and zero-map capsule failure separately.

No stretched outgoing pole has yet been measured. Do not transport a compressed
pole identity, promote native RESON fits to pole targets, infer exact Lambda from
C2v labels, or call real-axis consistency an independent exterior solution. The
next decision is one **own-DZ B1 finite pole attempt**, not both bases or a sweep.

## Saved-only seed reconnaissance

[Reproducible saved reader](../../projects/ukrmol_co/stretched_pole_seed_review.py)
authenticates the completed real allocation, then reads only each own native
99-point K export. For symmetric K, the eigenphase sum is

    phi(E_Ry) = unwrap(2 sum_j atan(lambda_j(K)))/2.

A single isolated Breit–Wigner phase without background would have maximum
dphi/dE=2/Gamma_Ry. This motivates a seed, not a width estimate with controlled
background/sampling uncertainty. The measured gradient maxima are:

|Own model|Native index|E /Ry|dphi/dE /Ry^-1|Heuristic full width /eV|
|---|---:|---:|---:|---:|
|DZ|16|0.06615|95.3085|0.285508|
|TZ|18|0.07350|82.5575|0.329605|

The corresponding grid energies are0.900017/1.000018eV. Maximum adjacent phase
changes0.362957/0.313654rad resolve this sampled rising feature, but cannot rule
out unsampled narrow features or other poles. These are not an outgoing-pole or
basis-sensitivity result. The original114 route S holdouts are different data
and remain unanalyzed. The lowest own closed thresholds are0.371824846Ry (DZ)
and0.370734394Ry (TZ), above the entire native grid.

Choose DZ seeds z=(0.033075,-0.00525) and(0.033095,-0.00525) **incident Hartree**,
with fixed rectangle Re(z)=[0.0225,0.0525], Im(z)=[-0.010,-0.002]Hartree. The seed
real part is half the native gradient-peak energy; Im(z) is minus one quarter
of the heuristic full width in Ry. No compressed result enters these choices.
The rectangle is deliberately broader than a grid-bin diagnostic, not a fitted
confidence interval. A seed/box exit is inconclusive with no reseeding or rescue.

## Sheet, units and contour admissibility

Use the native radial equation at energy2z (Ry), threshold-relative momenta on
the six-ground-channel outgoing sheet; keep reference0.06615Ry fixed. The inner
surface sum uses total Hartree energy E_ground+z, never the native radial energy
as an inner denominator. The whole rectangle lies below the first closed target.

For open k=sqrt(2z), exposure requires theta > -arg(z)/2. Over the whole box,
the largest required angle is half atan(0.010/0.0225), about11.98degrees. Both
15 and20degrees expose it; no branch point is crossed. Closed momenta at these
energies decay on both contours. The much larger closed/open asymptotic rate
separation here supports testing these angles, not a coupled tail-error bound.
Every root, contour and derivative evaluation must remain inside the frozen box.

Before complex maps, compare own native DZ index16 real S at120bohr/15degrees
and120bohr/20degrees, order10/rtol1e-12, with the original1e-6 S/phase/invariant
and1e-7 flux gates. Failure stops before root work. These are new near-feature
checks; previously qualified DZ80 is not repeated. Off-axis amplitudes must be
explicitly continued canonical Jost rows; never call real-axis `physical_s`
with a complex boundary energy or enforce normalization by projection.

## Frozen numerical and character witness plan

[Unreleased proposal](../../projects/ukrmol_co/stretched-dz-pole-proposal-20261102.json).
Reuse existing `pole_reference` bordered-root, projected resolvent, nested moments
and derivative-residue algorithms; add only an explicit own64/400 adapter and
prerequisite/owner routing. No new eigensolver or generic solver framework.

Three root settings separate endpoint and contour dependence:90/15,120/15,
120/20 bohr/degrees, all order10/rtol1e-12. Root evaluation caps16/12/12. Require
the same-root position **and full width** to agree within0.0001eV across settings;
the physical model comparison later retains0.05eV position/5% full-width with
0.001eV floor/0.05rad phase budgets. Do not use model budgets as numerical gates.

Use one fixed0.001Hartree circle,32 contour nodes with the16-node rule reusing
even nodes, eight projected interior columns/seed314159, and derivative steps
1e-5/5e-6Hartree. Retain the existing unit winding, isolation/singular-value,
interior/open signal/rank, nested-rule, root/moment agreement and1e-3 derivative-
residue gates. Capture every Q/Y/A and projected resolvent/open residue. These
are different constructions sharing the same native model/exterior ODE, not
independent exterior engines or an omitted-spectrum bound.

Require all four open-residue witnesses to pass the existing physical-Pi support
screen on native ground-target labels. Verified B1 open channels have(l,|m|)=
(1,1),(2,1),(3,1),(3,3),(4,1),(4,3). Non-|m|=1 leakage <=1e-4 passes, >=1e-2
rejects; intermediate is inconclusive. Rank/transpose ratios<=1e-3, signal/full-
width>1e-6, and witness leakage/profile/norm changes<=1e-5/1e-4/1e-3 are required.
Residue support is not a probability, exact Lambda proof, cross-basis wavefunction
overlap or global geometry continuity. Numerical/isolation failure is inconclusive;
wrong character rejects this candidate without trying another seed or B2.

## Cost, stop and following decision

At most79 maps: two near-feature real checks,40 root calls,32 contour nodes,
four derivative calls and one amplitude-centre call. Eight physical cores/8GiB/
one library thread per worker/eight contour workers; persistent no-network source-
read-only owner. The measured400-channel real maps take roughly9–19s each;
forecast600–2400s, hard7200s science/7800s owner,90s/50000 RHS per map. Reserve
at least16GiB available RAM/2GiB free persistent disk; measure actual costs.
Implementation60minutes; source/tests/exact adapter/owner/contract freeze and
push before any root/near-feature map, with source/input closure gates. Check
at least every60s. No dependency/image build, paid host or native/electronic repeat.

After public real-completion and seed receipts, explicit authorization is still
required. Any failure preserves original exits/partial captures and stops without
retry, expanded box, alternate seed, angle/order/tolerance/RHS ladder or extra
quadrature. A full pass supplies one locally isolated physically visible pole of
this finite DZ model only. Then separately decide a TZ counterpart using its own
data/character/local identity. No automatic peer, B2, geometry, production-model,
curve, fit, factory target or114-holdout analysis follows.

## Separate execution authorization —2 November2026

The user explicitly authorizes this bounded attempt after the public real
completion and source-only design. The historical proposal remains unchanged;
the [separate execution contract](../../projects/ukrmol_co/stretched-dz-pole-execution-20261102.json)
binds its bytes and all seven own-native/public prerequisites. The adapter
populates the strict own64/400 model state consumed by the existing
`compressed_pole_reference.calculate` function, never its compressed loader or
experiment owner. Bordered search, contour/residue/character functions and Jost
propagation remain byte-identical. Two real checks execute before complex maps.
The new eight-core owner verifies complete immutable source/mount/resource caps.
Concurrent contour failures cancel not-yet-running work; already running bounded
maps finish as preserved partial evidence, never as a retry or replacement.

### Observed adapter stop; numerical attempt unstarted

[Independent public evidence](../../projects/ukrmol_co/qualification-evidence/stretched-dz-pole-stop/README.md)
preserves the authorized `d04d3b0` allocation, all original exits2 and zero maps.
The own numerical model has no nested `interface` dictionary. The adapter asks
for that absent metadata after its own-array checks, producing the recorded
`KeyError('interface')` before near-feature guards or complex use. The26 passing
controls did not test the successful own-model return schema. This is an adapter
implementation failure, not scientific evidence against the pole or model.

The six native labels independently reproduce from the unchanged raw interface.
The [unreleased correction proposal](../../projects/ukrmol_co/stretched-dz-pole-adapter-repair-proposal-20261102.json)
uses the existing native metadata reader and requires actual-return-schema and
positive saved-native-load controls before a separately frozen fresh release.
No numerical method, seed, box, gate or cap is to change. Do not alter/restart
the stopped capsule or infer correction permission from the prior attempt.

### Separately authorized adapter correction

The user explicitly authorizes the metadata-only correction after the public
zero-map stop. A [separate immutable release](../../projects/ukrmol_co/stretched-dz-pole-adapter-release-20261102.json)
inherits every original numerical setting and binds the failed-stop reading.
Labels are read through existing `read_interface` on the same pinned native
bytes; the numerical model's eleven-field return schema is not changed. Positive
actual-schema and wrong-label controls pass, and the corrected adapter loads
the five hash-pinned saved native/receipt files (64/400/2048,99 K references)
successfully before freeze, with zero ODE/native calls. Original algorithms,
seeds/box/gates/caps and old capsule/exits2 stay unchanged. Only contract-path
routing/provenance change in the eight-core owner; no scientific fallback.

## Qualified finite own-DZ result —2 November2026

[Independent fresh-public evidence](../../projects/ukrmol_co/qualification-evidence/stretched-dz-pole/README.md)
qualifies the corrected `fbd16ee` allocation:48 maps/all original exits0,
incident pole0.03359318355534015−0.00490773471261789iHartree, position0.914117eV,
full width0.267093eV. Two own-native near-feature real checks pass, all three
root settings accept the same coordinate with sigma_min about3e-11, contour
windings16/32 are1/1, and root–moment full-width disagreement is1.58e-11eV.
Recorded zero coordinate differences are not zero uncertainty: residual,
contour and derivative gates remain required. Derivative residue changes
5.04e-6/1.26e-6 pass1e-3. All four Pi-support witnesses pass, reference
non-|m|=1 leakage2.03e-11, profile/norm changes1.34e-10/4.70e-6.

Eight fresh public payloads authenticate79 raw files/48 executed Git sources/
eight inputs/nine source prerequisites/2044 samples/14 observations. Independent
NumPy reconstruction uses saved native sums/Jost matrices, no new ODE/native
call. Profile409.613145s/CPU1194.174534s/3.002216GiB kernel peak. The original
adapter/40degree/capsule failures remain preserved. No job/monitor remains active.

This removes the uncontrolled phase-width diagnostic for **one local pole of
the finite own-DZ model**. It does not release an electronic model, absolute
CO accuracy, full-S embedding, omitted channels/spectrum, a Pi pair, global
identity, a width curve, fit or production sweep. The114 S holdouts stay untouched.

## Prospective own-TZ counterpart: separate decision

The [unreleased TZ proposal](../../projects/ukrmol_co/stretched-tz-pole-proposal-20261102.json)
inherits the DZ numerical allocation unchanged but selects its own native
phase-gradient index18, E=0.0735Ry, seeds(0.03675,−0.00606) and(0.03677,−0.00606)
incident Hartree. Its own diagnostic width0.02422555099992847Ry motivates the
imaginary seed; neither the DZ pole nor a compressed state is transported.
Box Re=[0.025,0.056], Im=[−0.011,−0.0025]Hartree lies below the first closed
threshold0.370734394Ry; its worst open exposure angle is below15degrees.
Near-feature own-TZ index18 real guards at120/15 and120/20 precede complex maps.

Same strict own64/400/2048 array/receipt/label gates, positive saved-native adapter
load, root/contour/derivative/character witnesses and79-map/8core/8GiB/7200s science/
7800s owner caps apply. Reuse the public DZ pole and both real prerequisites,
no successful calculation repeat. Explicit permission and a separate exact
source/contract/tests/owner freeze are needed; any failure stops without rescue.

A pass would supply a second finite local Pi-support pole, not cross-basis
wavefunction identity. Review open-label/residue compatibility and isolation on
saved evidence before interpreting a counterpart shift; a different character
or ambiguous identity is inconclusive, not a basis rejection. The eventual
model comparison retains0.05eV position and5% full width with0.001eV floor;
the numerical budget stays0.1meV. The comparison changes atomic basis plus its
optimized orbital frame, not isolated basis-size convergence. No automatic
B2/geometry/curve/fit/sweep or holdout analysis follows even from two pole passes.

### Separate own-TZ authorization

The user explicitly authorizes the own-TZ counterpart under the
[new execution contract](../../projects/ukrmol_co/stretched-tz-pole-execution-20261102.json).
The strict shared driver accepts only the two declared finite own64/400 B1
models; it selects the byte-pinned own-TZ proposal, receipt and native arrays,
never substitutes DZ inputs. All original numerical algorithms/gates/caps are
inherited unchanged. A qualified public DZ pole is a prerequisite, not its seed.
Positive actual own-TZ saved-native loading (including its99 K references) and
new receipt/seed/unit controls precede freeze. Next independently reconstruct
the actual pole-or-stop outcome, then review saved necessary counterpart
evidence before a basis-plus-frame shift interpretation. No numerical successor.
