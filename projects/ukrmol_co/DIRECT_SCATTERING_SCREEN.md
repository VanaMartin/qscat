# Direct fitting route: bare scalar full-S family rejected

## Decision — 23 October 2026,10:12 UTC

**Do not build a full-native-S fitter around the existing bare scalar radial
potential.** All three saved CAS11/nine-A2 anchors, both Pi sectors, reject even
the larger family of arbitrary diagonal unitary matrices independently at each
energy. Off-diagonal scattering alone is too large; changing Gaussian depths,
ranges or radial shells cannot remove that structural mismatch.

This eliminates a route before implementing its fitter. It does **not** reject
the effective two-dimensional potential goal, the pole-targeted factory,
anisotropic/coupled potentials or a scalar model embedded in a known coupled
exterior. It does not choose the production electronic model.

## Measured necessary-condition result

The prospective metric is the maximum relative eigenphase of two unitary S
matrices, not their eigenphase sum. Required distance≤0.05rad corresponds to
spectral chord error≤2sin(0.05)=0.0999583. The off-diagonal column norm supplies
a rigorous lower bound for **every** diagonal model. See the
[derivation](../../docs/physics/co-direct-scattering-representability.md).

| Anchor / R,bohr | Rejected development points per sector /80 | Inconclusive /80 | Largest lower bound on relative eigenphase distance |
|---|---:|---:|---:|
|Compressed /1.900|52|8|0.491729rad|
|Equilibrium /2.1323|64|0|0.523269rad|
|Stretched /2.500|80|0|0.439228rad|

B1/B2 give matching counts. Across six cases: **392 rejects /16 inconclusive
/72 constructive pointwise admissions**,480 total. Admission only means some
diagonal unitary fits that one matrix, not that any shared physical potential
or correct width has been found. The114 prospective route holdout S matrices
are **not analyzed** or captured. They are not globally unused historical data.

The full-S criterion is newly declared for this feasibility screen. Rejections
do not retroactively change earlier eigenphase-sum gates or diagnose an error
in native scattering. Known polar-exterior coupling can contribute to the
mismatch; the screen does not separate exterior and short-range sources and
cannot reject a tail-dressed scalar model.

## Scope, resources and reproduction

Source/tests/design/contract freeze **`8d0cbc8`** precedes execution. Completion
**10:08:45UTC**, original scientific/command/container/owner **1/1/1/1**: an
intentional scientific rejection, not exception, timeout or OOM. All six cases
complete. Container wall **0.87157s**, Python **0.20946s**, process peak RSS
66,364KiB (**.06329GiB**); no sampled cgroup peak/aggregate CPU measurement is
claimed. Four physical cores0–3 reserved, one BLAS thread,4GiB cap.

All123 original inputs, two prerequisite pins and12 frozen source entries
authenticate. The old saved-matrix aggregate false/five-A2 unset outcomes remain;
only its six qualified nine-A2 cases are used. **Zero native calls, propagations,
potential fits, pole searches, QC, target solves, properties, new geometries or
inner Hamiltonians.**

[Fresh public replay](verification-direct-scattering-screen-public-20261023.json)
authenticates **6 payloads /24 files /12 source entries**, downloads two prior
payloads and verifies all123 inputs. Independent record decoding builds S from
K eigenvectors instead of the producer's Cayley solve; a Hermitian residual-Gram
eigenvalue replaces its singular-value norm. All480 matrices/decisions reproduce,
maximum bound difference **1.26565e-14**. No producer numerical module or native
program runs; shared saved data/NumPy are not an independent scattering engine.
Original exits stay1. [Fetch/replay](qualification-evidence/direct-scattering-screen/README.md).

Persistent Sadaharu source/output:
`prepared/direct-scattering-screen-source-20261023/`; publication mirror
`prepared/publication-direct-scattering-screen-20261023/`. No job remains.

## Cycle disposition and single missing ingredient

**Stop the bare full-S fitter; return to the physical target definition.**
The scalar-potential goal needs a justified **channel-to-effective-model map**:
either a channel-aware/polar-tail embedding defining which native matrices a
forward model predicts, or an explicitly reduced resonant target with a
controlled relation to physical coupled scattering. Choosing an eigenchannel
at each energy, subtracting arbitrary background or fitting an eigenphase sum
does not itself provide that map or an identifiable width.

Keep the pole-targeted factory as the integration baseline, not a newly qualified
CO solution. Restricted extraction still lacks the independent physical
background/pole information needed to overcome earlier counterexamples. Neither
a qualified width curve nor a ready direct-S shortcut exists yet. Do not
auto-implement a coupled fitter or restart general complex propagation in
response to this rejection.

**Next bounded deliverable:** specify one channel-to-effective-model target map
and falsifiable validation observable from saved anchors, including retained
background/tail physics and deliberate reduction. If this cannot be specified
within one working cycle, report the modeling assumption needing a decision
instead of launching another ladder. Later numerical work needs a separate
finite contract. New50 stays parked; model, continuity, neutral-reference and
production-pilot gates remain open.
