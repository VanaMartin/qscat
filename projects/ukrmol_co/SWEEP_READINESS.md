# CO: gates before a production sweep

**No production model, potential fit or sweep is released.** Engine and
finite-method successes are not substitutes for these physical gates.

| Gate | Current evidence | Still required |
|---|---|---|
| Reproducible target and engine | Common state-averaged orbitals, independent roots/dipoles, CI/spin/Pi checks; source-built serial/MPI controls | All used model/geometries pass their own frozen checks |
| Affordable numerical scattering | Qualified selected-spectrum exports and 2048-root working baseline with refinement evidence | Numerical/extraction uncertainty below the physical comparison budget over the chosen domain |
| Common atomic/neutral zero | QZ/5Z finite atoms complete; both finite-basis diagnostics reject overall | Explicit correlation/reference treatment and remaining alignment error; no automatic 6Z/CBS inference |
| Neutral dissociation | 52 points over 1.9–2.5 bohr pass relative-basis/interpolation gates; four RHF attempts at 3/4 bohr fail stability | Stable spin-consistent continuation, overlap with accepted well data, fragment/size-consistency and tail checks |
| Electronic model | Compressed CAS11/CAS12 and stretched DZ/TZ pole sensitivity reject; channel/ensemble effects measured | One justified model with basis/active-space/channel/ensemble uncertainty, not selection by literature proximity |
| Continuous resonance/bound branch | Local finite DZ/TZ poles at 2.5 bohr and retained-target correspondence qualify | Global resonance identity, consistent thresholds and bound continuation; unresolved width is not zero width |
| Restartable geometry pilot | No production prescription frozen | Exact basis/core/active/ensemble/channels/roots/continuum/energy/extraction/domain, held-outs, restart checks and measured cost forecast |
| Predictive 2-D potential | No released fit | Freeze training once, validate excluded geometry/energy scattering and QSCAT discretization; nuclear check before attachment claims |

## Quantitative boundaries

Physical budgets stay **50 meV position**, **5% full width with a 1-meV floor**
and **0.05 rad in the specified phase metric**. Same-model numerical pole
budgets stay 0.1 meV and 1e-6 rad. Study-specific representation, flux, K and S
checks remain distinct. No lower-accuracy exploratory tier is adopted.

The neutral relative-energy pilot's 20-meV basis and 1-meV interpolation budgets
do not qualify absolute atomic alignment. Atomic 10-meV EA / 20-meV fixed-R
offset diagnostics reject overall; offsets at R=2.1323 are not De/D0 or initial-v
thresholds. Full correlated term purity, CBS and absolute accuracy remain open.

All **114 S holdouts** and the specified equilibrium outgoing-pole exclusion
remain untouched. Keep nine as working reference, new50 parked and the accepted
2048-root baseline. No silent threshold-grid truncation, root/basis/channel
escalation or augmented inner allocation follows from this checklist.

## Immediate disposition

One bounded saved neutral/correlation/reference-treatment decision, zero new
quantum calls, under [production decisions](PRODUCTION_DECISIONS.md). No active
job or monitor. [Settled results and blind alleys](../../docs/physics/co-results-and-rejected-routes.md)
explain the scientific stops; [the evidence index](qualification-evidence/README.md)
preserves the reproducible records. Historical queue items live in the optional
[history archive](qualification-evidence/campaign-history/README.md), not here.
