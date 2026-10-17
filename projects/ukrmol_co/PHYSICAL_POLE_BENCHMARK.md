# Physical boundary/pole benchmark passes — 17 October 2026

**The first physical benchmark passes all nineteen frozen checks.** Independent
numerical boundary matching reproduces a real finite-range radial shell's exact
scattering and nine outgoing resonance poles. All **36 real-axis profiles /
81 shooting solves** are retained, with no reseeding or changed tolerance.
Original scientific command/container exits are **0/0**.

The previous local phase counterexamples showed why a flexible analytic/unitary
fit cannot establish a width. This benchmark instead starts from a self-adjoint
Hamiltonian, including its regular-origin condition, real derivative jump and
outgoing-wave boundary. It tests physical conventions before applying them to CO.

## Independent mathematics and frozen scope

The [design](../../docs/physics/co-physical-pole-benchmark.md),
[contract](physical-pole-contract.json), Python implementation and three tests
were frozen at **`3d94f49` before the full validation matrix**. The Hamiltonian
is `−½d²/dr² + gδ(r−a)` on r>0 with u(0)=0, a=5 bohr and λ=2ga=0/2/10/30.
It conserves current across its real shell and has a physical outgoing
resonance series. This scalar s wave is a diagnostic model, not a scalar-p-wave
approximation to the six-open-channel polar CO problem.

Two independent routes are compared:

- **Oracle:** exact Jost ratio `S=F(−q)/F(q)`,
  `F(q)=1+λ exp(iq) sin(q)/q`, q=ka. Its pole positions come from the derived
  closed form `q_n=i[W_{−n}(λ exp λ)−λ]/2`. Computation at 60/90 decimal digits
  agrees within **5.33e-61**, and the unmultiplied Jost residual stays below
  **5.55e-60**. The principal-W branch's apparent q=0 root is removable and
  explicitly rejected: `F(0)=1+λ`.
- **Numerical route:** DOP853 integration of the regular radial equation, an
  exact shell derivative jump, and value/derivative matching to free waves.
  Root solves use only that outgoing residual and fixed seeds
  **`q=nπ−0.3−0.2i`**; they receive no oracle pole. Physical roots must have
  positive real q and negative imaginary q, positive real energy and width,
  and a nonzero incoming residual so that numerator cancellation cannot mimic
  an S pole.

All three matching radii B/a=1/2/4 and ODE relative tolerances 1e-5/1e-8/1e-11
are used. Each real profile contains **107 points**, including π/2π/3π
wavefunction-node controls and three threshold points: **3,852 real samples**.
The homogeneous pair `(u, B u'−b u)` avoids logarithmic-derivative divisions at
nodes. Three Bloch coefficients and a common pair rescaling check the defined
boundary representation without changing the physical Hamiltonian.

## Convergence and measured accuracy

| ODE relative tolerance | Maximum complex S error | Maximum pole q error |
|---|---:|---:|
| 1e-5 | 2.17349e-4 | 5.67359e-6 |
| 1e-8 | 2.57166e-7 | 6.58875e-9 |
| 1e-11 | **2.58327e-10** | **6.49474e-12** |

These are errors against the independent oracle, not successive differences.
Each level is retained; errors decrease strongly under refinement. Finest
maximum errors and the prospectively frozen budgets are:

| Check | Measured | Budget |
|---|---:|---:|
| Phase modulo π | 1.292e-10 rad | 1e-7 rad |
| Pole position | 6.054e-11 eV | 1e-4 eV (0.1 meV) |
| Full width | 1.125e-10 eV | 1e-4 eV (0.1 meV) |
| Outgoing shooting residual | 4.090e-12 | 1e-7 |
| Pole q spread across matching radii | 6.202e-13 | 1e-7 |
| Scattering-length ratio versus λ/(1+λ) | 3.951e-11 | 2e-8 |

Across **all** real-axis refinement levels, unitarity error `||S|²−1|` is at
most **8.89e-16** and Bloch/rescaling discrepancy **5.98e-16**, both below
1e-10. The free control gives S=1. The smallest finest incoming residual at an
outgoing pole is **0.1988**, exceeding the 0.1 noncancellation threshold. Every
root solver succeeds at every rung; the full report retains each status/message,
evaluation count, seed, recovered pole and exact counterpart.

Mandatory negative controls change the real-axis S at q=1.3, λ=10 by
**0.51199** for a reversed shell jump and **0.22569** for a missing factor two,
both above the required 0.01. A wrong real shell can still be unitary: this is
why oracle agreement and boundary/jump conventions are separate tests from
flux conservation.

## Reproduction, resources and next gate

[Public evidence](qualification-evidence/physical-pole/README.md) preserves nine
content-addressed payloads, including **ten source files**, the independent
high-precision dependency wheel and host/local reports. Fresh public replay
reproduces all nineteen checks, all 36 profiles and all 81 successful pole
solves. Maximum host/public q drift is **3.80e-12**, and maximum position/width
drift **3.60e-11 eV**. See [quoted verification](verification-physical-pole-20261017.json)
and [public replay verification](verification-physical-pole-public-20261017.json).

Sadaharu cost: **6.722 s wall / 0.051 GiB sampled kernel high-water**, reserved
physical cores 0–3, one BLAS thread, 4-GiB memory limit. There are 32 samples at
0.2-second intervals; the high-water observation is not a post-exit peak.
The initial dependency probe's original exit **1** (`mpmath` absent) remains
recorded. Its pinned pure-Python wheel is mounted read-only without modifying
the engine or pinned image. Source archive SHA256:
`999a1ca1a2be34b86d5a50f3ae6019b3dc1ec394dec54f69753f5743b93bcb05`.
Persistent files/logs are in `prepared/physical-pole-source-20261017/`.
**28 focused tests pass**. No native file is read or inner solve launched.

**Next:** design a coupled polar-tail benchmark with independently checkable
matrix scattering, open/closed thresholds and outgoing continuation, then
qualify its real-axis and complex-energy numerics. Separately derive the
actual native boundary-amplitude normalization, Bloch and Hartree/Rydberg
conventions from the writer/assembler before consuming CO boundary data. Full
native K-matrix reproduction and radius/step/tail convergence remain required
before a native complex-pole search. The scalar pass is not that qualification.

The five-A2 label reading remains unset with its original aggregate false/exit
1. Compressed/equilibrium extraction fidelity rejection, stretched channel
sensitivity rejection, active-space/model uncertainty, continuity and neutral
domain limitations are preserved. Production remains gated.
