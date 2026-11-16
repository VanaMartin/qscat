# CO: established results and rejected routes

## Scope

This experiment supplies a reproducible UKRmol+ fixed-nuclei workflow and
finite-model diagnostics, **not a production CO potential**. Its intended target
is a neutral curve and the lowest doublet-Pi resonance/bound-anion continuation,
on a consistent atomic zero, followed by independent scattering validation of
a QSCAT two-dimensional model. Low-Pi vibrational excitation or channel-specific
vibrationally assisted attachment does not establish total thermal CO attachment
near 10 eV, where other resonance families matter.

Lengths and total energies are in bohr and Hartree. Presentation energies are
in eV; widths are full widths, with a pole written as position minus half-width
times the imaginary unit. A native phase-fit candidate is not an outgoing pole.

## What is established

| Contribution | Verified scope | What it does not establish |
|---|---|---|
| State-averaged target builder | Common multi-spin/multi-irrep PySCF orbitals, independent UKRmol roots/dipoles, explicit occupations and spin/Pi/CI gates | Electronic-model accuracy or continuity merely from converged orbitals |
| Selected-spectrum scattering | Native sparse boundary export and dense-oracle comparisons; 2048-root working baseline with higher-root controls | A universal root count or converged active space/basis/channel prescription |
| Independent outer-region controls | Native normalization, K/S/phase reconstruction, exact shell and coupled-tail controls, captured-operand replay | Uniform real/complex convergence outside the tested domains |
| Extended residual kernel | Optional Rust/PyO3 and C long-double dot products with NumPy differential oracle | Portable extended arithmetic on every host or a replacement for model qualification |
| Correlated neutral well | 52 frozen-core RHF/CCSD(T) QZ/5Z points, 1.9–2.5 bohr | Molecular dissociation, absolute atomic alignment or unrestricted extrapolation |
| Atomic endpoints | C(3P), O(3P), O-(2P), C-(4S), matched aug-QZ/5Z ROHF-reference UCCSD(T), one frozen 1s pair per atom | Absolute accuracy, CBS convergence, full correlated orbital-term purity or a molecular size-consistent limit |
| Stretched finite poles | Own DZ and TZ CAS12 B1 poles at 2.5 bohr, independently reconstructed with local Pi/residue controls | Cross-basis unique N+1 identity, a continuous width curve or a production model |

The target implementation and conventions are in
[state-averaged targets](co-state-averaged-target.md) and
[electronic qualification](co-electronic-qualification.md). Numerical evidence
and fetch/replay instructions are indexed under
[qualification evidence](../../projects/ukrmol_co/qualification-evidence/README.md).

## Rejected inferences and parked branches

"Rejected" refers to the named prescription or inference, not to all methods of
that class. "Blocked" means that a prerequisite failed; it is not a physical
model comparison. No failed record is retrospectively promoted by a later repair.

| Route or shortcut | Evidence | Settled disposition |
|---|---|---|
| Choose SEP by agreement with one literature value | Changing TZ virtual space 39→60 shifts equilibrium fitted position about 0.294 eV | Reject agreement-as-convergence; retain SEP as a diagnostic |
| Treat ground-only CASSCF as the multi-state target | Different orbital objective and excited spectrum; basis/external-space sensitivity persists | Not the qualified state-averaged prescription |
| Repair every orbital solution by increasing CI space | Some final orbitals have converged energies but failed root flags/spins; compressed CAS12/TZ and its A2 pilot fail required CI/action gates | Stop those releases; no full-QC retry ladder or claim of a physical TZ rejection |
| Adopt the equilibrium 50-component ensemble automatically | Physical retained A1 subspaces change, including root5→root6 exchanges | New50 parked; changing the ensemble is not simply adding scattering channels |
| Five A2 channels are adequate at every anchor | Five→nine phase changes: compressed 0.0115944, equilibrium 0.0163665, stretched 0.0679062 rad | Stretched rejects 0.05 rad; nine is the working reference, not full channel convergence |
| More accurate phases identify a unique fitted width | Endpoint-isolated analytic audit finds 94 known ±10% wrong-width representations passing all unused phase points within 0.001 rad | Reject phase prediction alone as width qualification; these toy backgrounds are not proven causal CO Hamiltonians |
| Fit the bare saved full S matrix directly | Saved-data representability screen rejects the unqualified shortcut | No bare full-S potential-fitting release |
| Infer outgoing-map accuracy from S/phase alone | Original 40-degree real guard fails flux/Wronskian checks despite S/phase agreement | Original stop remains; later separately qualified 20-degree prescription does not erase it |
| CAS11/CAS12 compressed poles are insensitive | Position changes 86.774662 meV, beyond 50 meV; full-width change 16.066613 meV passes | Reject model insensitivity; no automatic B2 extension |
| Stretched DZ/TZ poles are insensitive | Position changes 103.369940 meV; full width changes 44.171679 meV (16.54%) | Reject both physical budgets; target correspondence is necessary, not N+1 wavefunction identity |
| Continue aug-TZ with the unchanged fixed-six-open export | Triplet thresholds at 4.941440507 eV cross requested 4.95/5.00-eV rows | Block before inner allocation; retain qualified properties and 87,198-state preparation; no silent grid truncation |
| Extend restricted neutral CCSD(T) to dissociation | All four RHF attempts at 3/4 bohr in aug-TZ/QZ fail external stability before CCSD(T) | Reject that continuation; select a stable spin-consistent correlation/reference treatment |
| Declare atomic endpoints basis-insensitive after QZ→5Z | O EA shifts 23.533750 meV against 10 meV; all fixed-reference offsets exceed 20 meV | Reject finite-basis insensitivity; no automatic 6Z or CBS extrapolation |
| Four-node cubic Er/cubic log(full width) is accurate throughout the local interval | Fixed midpoint prediction passes; fixed compressed-gap errors29.36meV/32.28meV exceed5meV/13.89meV, while both outgoing pairs qualify | Reject the candidate's representation accuracy, not the six valid pairs. Reuse saved rows for separately designed successor/validation; never relabel a validation-driven refit as prospective. [Evidence](https://data.qscat.org/co-compressed-gap-r20621125-outgoing-20261113/assessment.e4e8336303f5.md) |
| Six-node PCHIP meets the frozen position budget | Both prospective position errors (+20.135/−14.527 meV) exceed 5 meV; both full-width tests pass their own budgets | Reject the candidate's stated position accuracy, not its eight valid local pairs. Old validation may become development only for a separately frozen recipe with genuinely new validation; original predictions/scores remain unchanged. [Frozen-score evidence](https://data.qscat.org/co-successor-r198105625-outgoing-result-20261114/assessment.12500e437981.md) |

Detailed controls remain in [phase identifiability](co-phase-identifiability.md),
[direct-S representability](co-direct-scattering-representability.md),
[compressed active-space comparison](co-compressed-cas12-model-control.md),
[stretched pole design](co-stretched-b1-pole-design.md), and
[diffuse-flexibility limits](co-stretched-diffuse-flexibility.md).
These are measured limits, not invitations to reopen supporting campaigns.

## Neutral and atomic energy alignment

The near-well QZ→5Z relative-energy change is at most 17.2433 meV; cubic held-point
error is about 0.114 meV. Both pass their 20/1-meV pilot budgets. Each curve uses
its own molecular energy at R=2.1323 as zero; small relative errors do not imply
an accurate absolute molecular-to-atomic offset.

Atomic calculations use occupation-constrained ROHF followed by explicit
ROHF-to-UHF conversion and UCCSD(T), **without further UHF optimization**.
Reference internal stability is checked; external ROHF stability is unavailable,
not passed. State, reference-spin and occupation diagnostics do not prove full
correlated orbital-angular-momentum purity.

| Quantity | Aug-QZ | Aug-5Z | 5Z−QZ |
|---|---:|---:|---:|
| O electron affinity / eV | 1.403072709 | 1.426606459 | +23.533750 meV |
| C electron affinity / eV | 1.235049235 | 1.240999940 | +5.950705 meV |
| C+O minus matching-basis molecular reference / eV | 11.115011294 | 11.167927472 | +52.916179 meV |
| C+O- minus matching-basis molecular reference / eV | 9.711938585 | 9.741321013 | +29.382428 meV |
| C-+O minus matching-basis molecular reference / eV | 9.879962059 | 9.926927532 | +46.965473 meV |

The predeclared EA diagnostic is 10 meV; the fixed-reference-offset diagnostic
is 20 meV. Carbon's EA passes its finite difference; oxygen and all three offsets
reject. **These offsets are not De, D0 or initial-v attachment thresholds.**
The molecular reference is not an optimized minimum or a vibrational level.
Spin-free calculated terms must not be scored as lowest-J experimental thresholds
without sourced corrections to both neutral and anion energies.

QZ costs 17.913795 s / 1.943348 GiB; 5Z costs 114.673835 s / 11.883816 GiB,
each on four physical cores with all original calculation exits zero. Container
kernel peaks are lifetime high-water marks, not independent per-atom allocations.
The separate initial 5Z startup stops before any integral/SCF/CC call because
Python versions serialize AST metadata differently. Exact calculation-source-byte
equality fixes that control in a separate release; it does not requalify the
original failure. All numerical, physical and resource gates stay unchanged.

See [atomic states, totals and zero conventions](co-atomic-energy-references.md)
and the [public matched-basis reconstruction](../../projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/README.md).

## No-repeat and reopening conditions

### Target coverage versus pole identity

Physical target-metric coverage is complete for the original four-geometry
CAS11/nine-A2/2048 path: 16 missing sector matrices qualified and eight final-gap
matrices were reused. Direct interior retained/buffer exchanges and weak
compressed correspondence leave unique N+1 continuation unqualified. Finite-step
core motion can attenuate overlaps; low target overlap alone does not prove a
wrong local pole or a discontinuous Hamiltonian. The then-clean final gap had
only two actual pole nodes, so it did not qualify a narrower fit at that
assessment. These historical findings do not invalidate later qualified pairs
or the completed prospective local representation. Reuse accepted QC, native
inputs, poles and metrics; do not change retained roots to hide exchanges or
repeat successful work merely for another pass. Reopen correspondence only when
a separately selected observable justifies a pole/N+1-specific or interior-
sensitivity test with a frozen outcome-to-decision rule. Scoped interpolation
acceptance remains distinct from this unresolved identity claim.

[Completed target-metric batch](https://data.qscat.org/co-common-frame-target-overlaps-20261113/results.tar.8161225e1209.gz),
SHA256 `8161225e1209da295bd61fccf97ecfe517907aee4a8560e9280bf6f02dfa57b0`.
The later [accepted local representation](#local-representations-accepted-scope)
does not settle this identity question.

## Local representations: accepted scope

**Scoped local representations are accepted over 1.9–2.5 real bohr only.**
The primary neutral is the saved degree-6 aug-cc-pV5Z frozen-two-core
RHF/CCSD(T) Chebyshev representation, zeroed to its own molecular total at
R=2.1323 bohr, not the minimum or fragments; QZ and MR results remain separate.
The CAS11/cc-pVDZ/nine-A2/48-target/2048-root resonance representation evaluates
unchanged not-a-knot cubics of Er and log(full Gamma), using eight development
nodes. Two separately frozen prospective tests pass in both B1/B2 sectors:
maximum observed position/full-width errors are 3.345/1.171 meV, within their
5-meV and max(1% predicted full width, 1 meV) gates. Extraction retains its
separate 0.1-meV budgets. Original cubic and PCHIP accuracy rejections remain
rejected; their qualified poles, frozen predictions and scores are retained.
These selected-point tests establish neither a uniform bound, electronic/global
accuracy, unique N+1 continuation nor full-S/nuclear predictions.
With z=Er−i Gamma/2, Vres_hybrid=V0_5Z+Er_native is explicitly mixed-method,
not a correlated anion total. Reuse the saved delivery/evaluator; stop compute.
Reopen only for a separately selected observable and bounded release, never to
rescue failed scores or fill an unapproved global model.

- [Evaluator usage](../../projects/ukrmol_co/README.md#using-the-local-representations)
  and [mathematical/API conventions](co-effective-resonance-target.md#saved-local-evaluator-not-a-factory-model).
- [Immutable local delivery](https://data.qscat.org/co-local-representation-delivery-20261115T2020/delivery-capsule.tar.437ea9e58c35.gz),
  SHA256 `437ea9e58c351b2a940bd2f57be05f2c09d34c3bb0c906735bd273b419c53555`.
  It retains the original candidate's historical false acceptance flag alongside
  the separate post-validation acceptance sidecar.
- [Pre-compression documentation archive](https://data.qscat.org/co-local-representation-docs-archive-20261115/documents.tar.d29aa7fc1fea.gz),
  SHA256 `d29aa7fc1feadb05d4f9507968b49845bf303a565932a11e140f13af982b3d35`;
  [fresh-download authentication](https://data.qscat.org/co-local-representation-docs-archive-20261115/verification/authentication.b9de58929e6b.json).
  This preserves the seven exact original pages, including local-only prose;
  original scientific inputs, failures, checkpoints and costs remain in the
  [external evidence catalogue](../../projects/ukrmol_co/experiments.md).

## Limits that survive cleanup

- Physical budgets remain 50 meV in position, 5% full width with a 1-meV floor,
  and 0.05 rad in the specified phase metric. Same-model numerical pole budgets
  remain 0.1 meV and 1e-6 rad. Per-study gates are not interchangeable.
- An empty phase fit or unresolved width is not a bound anion or zero width.
- Native pole energies are relative to their own target threshold. Combining
  them with an independently correlated neutral is a **declared hybrid**, not
  the calculated absolute anion energy of either method.
- All 114 reserved S-matrix holdouts remain untouched; the specified equilibrium
  outgoing-pole exclusion remains reserved. No fit or sweep is released.
- Accepted calculations, immutable inputs, source pins and original failures
  remain available. Launch receipts, host PIDs and repeated status figures are
  historical evidence, not current execution authority.

Global neutral dissociation/correlation, resonance-to-bound joining and independent
physical validation remain unresolved before a production potential exists.
They are distinct from the completed scoped local representation milestone.
