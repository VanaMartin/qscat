# Necessary representability screen for direct CO scattering fitting

## Decision and finite scope

Before building a direct-S fitter, ask whether the **bare isotropic scalar
potential family**, embedded in the fixed physical open-channel basis, can meet
a full-matrix scattering accuracy requirement. This is a necessary-condition
screen, not an inverse fit or a width extraction. It uses saved CAS11/nine-A2
data at1.9/2.1323/2.5bohr, B1/B2, with no new native calculations.

`FlexibleDiatomicModel.v_int(r,R)` is angular-independent and its existing
factory consumes pole/neutral/coupling targets, not native multichannel S.
An isotropic potential with free asymptotic channels produces diagonal S in the
physical angular basis. Relax that family generously to **any diagonal unitary
matrix independently at each energy**, without requiring one common radial
potential or equal phases for repeated l values. Rejecting even this larger
class rules out direct full-S fitting by the bare scalar family. It does **not**
rule out scalar resonance reduction, tail-dressed scalar models, anisotropic or
coupled models, or fitting selected observables instead of the full matrix.

The alternative extraction route currently lacks a physical background/pole
restriction that resolves its known identifiability failures. Merely adding
polynomial orders is ruled out by existing evidence. This screen is preferred
first because it can eliminate an incompatible direct-fitter specification from
saved data without a new propagator, potential optimization or complex pole search.

## Metric and bounds

All S matrices are dimensionless and flux normalized. Retain the actual six
ground-target channels and their saved real-harmonic basis. Do not diagonalize S
and call its eigenchannels an isotropic physical-channel representation.
For unitary S,D define the relative eigenphase distance

    d(S,D) = 1/2 max_j |Arg lambda_j(S† D)|,
    ||S-D||_2 = 2 sin d(S,D).

Adopt a **new prospective full-S representability requirement** d≤0.05rad,
equivalently tau=2sin(0.05). Its angular scale matches the campaign phase budget,
but it is a stronger/different observable than the eigenphase-sum gate; do not
retroactively fail previous calculations under this new definition.

For every diagonal D and basis vector e_j,

    ||S-D||_2 >= ||(S-D)e_j||_2
               >= sqrt(sum_{i != j}|S_ij|²).

Thus L=max_j sqrt(sum_{i!=j}|S_ij|²) is a rigorous lower bound even for arbitrary
complex diagonal matrices. Choose D0_jj=S_jj/|S_jj| (1 when S_jj=0); then
U=||S-D0||_2 is a constructive upper bound for the relaxed diagonal-unitary class.
No optimizer is required. Both bounds are invariant to independent channel
rephasings/permutations; arbitrary channel mixing is not allowed.

- L>tau+1e-10: **reject** the bare diagonal full-S family at that energy.
- U≤tau−1e-10: a diagonal unitary representation exists at that energy; this
  does not establish a common potential, accurate poles or inverse uniqueness.
- Otherwise: **inconclusive**, with no automatic fitting/refinement successor.

Require finite6x6 matrices, unitarity/reciprocity defects≤1e-10, rtol0. Reuse the
previously validated native K reader and Cayley transform, while checking an
independent K-eigendecomposition reconstruction within1e-10. This is saved-data
algebra, not independent outer propagation or absolute electronic accuracy.

## Data separation and controls

Pin the original input manifest and completed saved-matrix report. Its aggregate
five/nine gate remains false: only the six qualified **nine-A2** cases are used.
Authenticate all123 original inputs before/after. Preserve original failed gates.
Each case has99 energies0.1–5 requested-eV; native conversion remains
0.0735Ryd/requested-eV. Verify actual units/decks, thresholds and physical labels
against the qualified case metadata.

Evaluate only zero-based indices i%5!=4:80 development points per case,
480 total. The remaining19/case,114 total, are reserved for a possible later
route test and are **not analyzed here**. They were used in historical studies,
so they are prospective holdouts for this route only, not globally untouched
data. Reading/checksumming complete files and checking energy-grid completeness
does not use reserved S values for this representation decision.

Tests first: diagonal unitary matrices including known scalar-pole factors;
orthogonally mixed negative controls; independent relative-eigenphase/chord
identity; invariance to rephasing/permutation; malformed/nonunitary input rejection.
Wrong-width diagonal factors must also pass the representation screen while
width qualification stays unset: this explicitly prevents treating representation
as identification. These algebraic factors are not certified causal CO potentials.

## Execution and stopping

Lifecycle: design → small project-only Python metric → analytic/differential tests
→ one measured saved-only application. No Rust or qscat promotion is warranted.
Freeze/push before execution. Expected60s, upper180s, hard cap300s, four physical
cores0–3 reserved, one BLAS thread,4GiB cap. Implementation/diagnosis is bounded
to this working cycle; a format/runtime failure returns to the goal review,
not an automatic repair campaign. Preserve source/input hashes, all matrices
actually analyzed, bounds/candidate phases, resource records and original exits.

All six cases are reported. Any development-point rejection makes the aggregate
bare-family verdict false and scientific exit1; an inconclusive aggregate exits2;
only all constructive admissions exit0. A failed numerical/input gate stops with
an error record rather than manufacturing a scientific rejection.

**Production consequence:** reject → do not build the bare scalar full-S fitter;
name the missing tail/coupling embedding or explicitly different reduced target.
Admit/inconclusive → return to route selection before a separately bounded common-
potential/identifiability test. No widths, new geometry sweep, model-convergence
claim, native propagation, inner Hamiltonian or reopening of new50 is released.
