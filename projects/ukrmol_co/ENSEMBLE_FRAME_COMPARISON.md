# Physical orbital-ensemble frame comparison

## Outcome — 21 October 2026,16:49 UTC

**The finite retained-subspace dominance screen rejects.** Both A1 sectors fail
the prospectively frozen minimum squared principal-overlap criterion0.5; the
other six sectors pass. All numerical gates pass, all eight sectors complete,
and original command/container/owner **1/1/1** are intentional scientific exits,
not a CI crash or timeout. No retained set or tolerance is changed after seeing
this outcome. The original seed-preservation rejection is a different result and
remains unchanged.

| Sector | Minimum squared principal overlap | Dominance screen |
|---|---:|---|
| Singlet A1 |0.0692964|Reject|
| Singlet B1 |0.7839631|Pass|
| Singlet B2 |0.7839631|Pass|
| Singlet A2 |0.9859719|Pass|
| Triplet A1 |8.20e-24, numerically near zero|Reject|
| Triplet B1 |0.5866026|Pass|
| Triplet B2 |0.5866026|Pass|
| Triplet A2 |0.9817204|Pass|

These values compare the **physical many-electron retained-five subspaces**,
including both frozen cores and core/active cross terms. They are invariant to
phases/unitary mixing within each retained set and are not normalized to discard
physical loss. The0.5 scale is a dominance screen, not a scattering error bound,
near-identity criterion or proof of electronic convergence. In particular, a
passing0.5866 direction is not a nearly unchanged state.

## Retained/omitted exchanges are resolved, not just mixed root labels

- Old **singlet-A1 root5** maps to **new root6**, squared overlap **0.9074664**.
  Its total projection into the new retained five is only **0.0794180**.
- Old **triplet-A1 root5** maps to **new root6**, squared overlap **0.9837919**.
  Its total projection into the new retained five is **8.22470e-6**.
- New **triplet-A1 root3** has only **0.3885705** total projection into all twelve
  old buffered roots. Its best/assigned old root10 scores **0.3629442**. It is
  therefore not established as a clean match to one available old buffered root.
- New **singlet-A1 root5** is also weak in the old buffer: total projection
  **0.2790905**. The bidirectional retained-five assignment scores its forced
  old-root5 match0.0787591; the row-best old root is3 instead. A constrained
  one-to-one assignment is not an identity certificate.

Across80 directed retained-root assignments: **2 weak /0 ambiguous /3
retained-to-extra triggers**. Buffered roots are diagnostic alternatives, not
new retained channels. Both directions assign only five mandatory rows; buffered
rows never compete for those assignments. Thus an arbitrary retained gauge
rotation alone cannot explain the A1 rejection. No buffer expansion, changed
selection, property export or new scattering calculation follows automatically.

## Numerical evidence and measured cost

Reuse the eight authenticated original twelve-root CI arrays at the exact old
scattering checkpoint. Generate only eight fixed-new-frame CI sets:13 A1 roots
and8 roots otherwise per spin, space160 /600-cycle cap, unchanged1e-12 energy /
1e-22 linear-dependence settings. The byte-identical qualified
`coverage_reference.state_residual` helper is reused without edits. All74 new
roots pass physical/penalized/penalty-action, spin, Gram and fixed-frame energy
gates. The old saved residual gates reconstruct from authenticated metadata;
their actions are not rerun. All new CI arrays persist for reuse.

| Numerical check | Maximum observed discrepancy |
|---|---:|
| New physical CI residual /Hartree |6.98481e-10|
| New penalized CI residual /Hartree |9.79196e-10|
| New spin-penalty action /Hartree |6.11475e-10|
| Schur versus full occupied minors |6.66134e-16|
| Forward versus reverse conjugate transpose |3.10862e-15|
| Fresh public-wave full-minor versus saved overlap |2.55351e-14|

Own-frame AO/MO/CI physical Grams, physical contractivity and metric reciprocity
also pass. Retained permutation/unitary mixing controls pass; retained-buffer
exchange and uniformly contracted0.7 overlaps reject. The analytic Hubbard
exact/perturbed residual controls pass. Complex full-Slater toy/oracle and cache
corruption controls are covered by21 focused tests (32 with related coverage/
binary tests, plus one local PySCF skip exercised remotely).

Source freeze **`869a588`** precedes the run **16:42:39–16:43:43UTC**. Cost
**63.503s**, sampled CPU **170.640s**, sampled kernel peak **529,281,024 bytes /
0.49293GiB**,315 samples. Eight new CI sets take57.660s in aggregate, overlap/
self/numerical checks4.920s; this is not an isolated implementation benchmark.
Four threads on reserved physical cores0–3 /4GiB. Forecast300s expected /900s
upper /1800s cap is not exceeded. All12 frozen sources /32 original input hashes
verify before/after. **Zero orbital optimizations, native engine invocations or
scattering inner solves.**

The [prospective design](../../docs/physics/co-ensemble-frame-overlaps.md) and
[machine contract](ensemble-frame-overlap-contract-20261021.json) remain frozen.
[Result/resources](progress-ensemble-frame-result-20261021.json),
[fresh numerical replay](verification-ensemble-frame-public-20261021.json) and
[public evidence/reproduction](qualification-evidence/ensemble-frame/README.md)
preserve the false decision/exits1/1/1. All11 public payloads fresh-verify;
18 prerequisite payloads reauthenticate. The collection preserves85 files /
108,765,937 uncompressed bytes. All32 pinned inputs map to existing public
prerequisite members; no old wavefunctions are republished. Saved-only replay
(`e066879`) takes8.990s on the Mac with **no CI/native solve**, using full occupied
minors instead of the Schur contraction and Hermitian Gram eigenvalues instead
of squared SVD. It shares the already validated full-minor cross-check routine;
it does not independently regenerate AO integrals or CI residual actions.

## Production decision

The new frame and explicit40-state properties remain independently qualified.
What fails is **the proposed clean same-physical-manifold ensemble contrast**
using five lowest roots in each sector on both frames. This is already a
decision-changing model discrepancy, not a reason to repeat QC or tune a solver.

Keep the original forty-component frame as the existing baseline. Distinguish:

1. A **combined prescription-sensitivity** comparison can intentionally compare
   lowest-five manifolds on both frames, but it varies retained physical channels
   as well as orbitals. It must not be reported as an isolated orbital-average
   effect or used to infer a channel-error budget.
2. A **matched-manifold** contrast needs a prospectively defined physical
   selection/projection prescription and adequate evidence on both frames.
   Replacing each new root5 by6 is not sufficient: new triplet-A1 root3 remains
   weakly represented in the available old buffer. No automatic root-count ladder
   or stronger identity claim is justified by these data.

The technical step common to either route is an immutable-frame continuum/native
preparation design, not another optimization. The target-only native layout is
not the scattering layout: the saved new `fort.16` is44,278 bytes versus the old
full-continuum12,983,436 bytes, and target CONGEN records differ by48 bytes. These
size observations do not establish conversion compatibility. `run.py`'s
`--target-initial-checkpoint` projects/reoptimizes; it is not an immutable import.
Qualify a read-only Molden/checkpoint reuse adapter and compatible integral/target
layout with explicit overlap/map checks **before** a new native preparation.

Old full-continuum integral-stage time11.020s and new target-only0.114s are merely
orientation; neither forecasts the whole preparation or scattering job. The old
dense full pipeline4.423h /25.201GiB is not a sparse/new-map cost estimate.
Select the scientific contrast, freeze the finite preparation stages and measure
their cost before requesting a sparse2048-root Hamiltonian. No dense scattering
production, phase/width claim, ensemble convergence or production sweep is
released. No scientific job remains running.
