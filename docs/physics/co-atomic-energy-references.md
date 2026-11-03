# CO atomic energy references and finite-basis endpoint controls

## Purpose and present status

The independent near-equilibrium neutral calculation is molecular frozen-core
RHF/CCSD(T), not a dissociation reference. Its accepted aug-QZ/5Z data span
1.9–2.5bohr. All four attempted3/4bohr RHF references failed external stability
before CCSD(T). At the initial review,no qualified atomic C/O/O-/C- totals or
fragment sums were found. The immutable pre-execution fields remain unset in the
[input/reference ledger](../../projects/ukrmol_co/atomic-asymptote-inputs-20261103.json).

**Completed first pilot:** [independent public reconstruction](../../projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/README.md)
now supplies all four finite-method aug-QZ totals and fragment sums. EAs are
O1.403072709/C1.235049235eV,not an absolute-accuracy or molecular-dissociation
qualification. Original exits0,17.913795s owner/1.943348GiB container high-water.
No supported bound-state continuation or practical matching radius follows from
atomic sums alone. Original source/protocol inputs and their earlier missing-
result status remain preserved in the executed capsule.

**Completed matched5Z control:** [public reconstruction](../../projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/README.md)
supplies5Z EAs O1.426606459/C1.240999940eV. Same calculation-source bytes,atomic
states,one frozen1s and ROHF-reference UCCSD(T),without further UHF optimization.
All original exits0;114.673835s owner/11.883816GiB/four cores. The original
zero-quantum AST-control failure is preserved separately;the corrected release
requires exact source-byte equality and passes eight actual-image controls before
science. No numerical/physical gate changed.

The predeclared **finite-basis**10meV EA and20meV atom-to-fixed-molecular-reference
offset diagnostics reject insensitivity: O EA changes23.533750meV (C5.950705meV),
and C+O/C+O-/C-+O offsets change52.916179/29.382428/46.965473meV. References
are each basis's molecular energy atR=2.1323,not optimized minima or nuclear
levels;these offsets are **not De,D0 or initial-v DA thresholds**. Small relative
well-shape errors do not qualify absolute atomic alignment. Two bases give no
CBS bound;stop endpoint-convergence inference,no6Z/automatic extrapolation.
Next select a neutral/correlation/reference treatment with a demonstrated route
to molecular fragment consistency;absolute accuracy and correlated L2 remain open.

The first endpoint pilot targeted four free atoms in **aug-cc-pVQZ**, matching one
existing molecular basis. It does not repeat the neutral well grid and does not
claim that restricted molecular CCSD(T) dissociates correctly. A second basis is
conditional on the first actual endpoint/EA results, not an automatic ladder;
the user-selected matched5Z control above completes that bounded comparison.
This note fixes the protocol; no new atom execution is released until implementation,
state checks, exact source/image/input pins, tests and a durable costed owner are
frozen and pushed.

## Atomic states and correlation convention

| Species | Charge | Electrons | Term targeted | 2S | Nominal valence configuration |
|---|---:|---:|---|---:|---|
| C |0|6|3P|2|2s2 2p2|
| O |0|8|3P|2|2s2 2p4|
| O- |-1|9|2P|1|2s2 2p5|
| C- |-1|7|4S|3|2s2 2p3|

These are the ground-fragment terms in Dora2016 p.4, Table3; its lowest-doublet-Pi
curve correlates with C(3P)+O-(2P) in p.5, Fig.1/discussion. This is the physical
target for our continuation, **not** an already demonstrated assignment of our
scattering branch. Neutral molecular X1Sigma+ requires a total-singlet coupling of
two triplet atoms, not two closed-shell fragments.

Use nonrelativistic spin-free Hamiltonian, spherical Gaussian functions, isolated
atoms without ghost functions, and one frozen occupied1s spatial orbital per
atom. Keep its mean-field energy in the total. This matches the molecular
valence-only correlation count, not an all-electron/core-correlated reference.

Proposed implementation: a constrained high-spin **ROHF reference followed by
PySCF UCCSD(T) on its explicit UHF representation**, with no additional UHF
optimization. Record the actual class/conversion, orbital energies and triples
convention. Do not call it closed-shell restricted CCSD(T), or assume that another
package's differently defined open-shell CCSD(T) is an exact numerical oracle.
The [software-only runtime inspection](../../projects/ukrmol_co/atomic-runtime-preflight-20261103.json)
confirms this route's source-level support in existing PySCF2.11.0: UCCSD calls
`to_uhf` for a non-UHF reference and conversion documents unchanged reference
energy/wavefunction. It also establishes that **external ROHF stability is not
available**. Run supported internal checks; external status stays null/unavailable,
never a manufactured pass. This inspection performs no atom, integral, SCF or CC
calculation and does not validate the numerical pilot. No image rebuild is needed.

Enforce occupation in an atomic D2h subgroup. The doubly occupied1s and2s orbitals
have Ag symmetry. For the three p irreps, high-spin occupations are respectively
`[1,1,0]`, `[2,1,1]`, `[2,2,1]`, `[1,1,1]` for C/O/O-/C-. Equivalent orientation
choices are degenerate; freeze one assignment rather than selecting it by energy
agreement. These constraints support the intended valence term, but multiplicity
and point-group counts alone are not a full correlated L2 proof. Report orbital
character/state diagnostics and any unresolved term mixing rather than claiming
exact term purity from the species label.

Minimum runtime checks: electron/charge/2S and irrep counts, overlap-rank/linear-
dependence diagnostics, intended occupied core/valence character, reference S2,
finite converged SCF and coupled-cluster iterations, finite triples correction,
and immutable input/source hashes. Record available reference stability diagnostics
with their supported meaning; a symmetry-constrained atomic reference is not the
same problem as the rejected molecular RHF singlet. Convergence proves a finite-
method energy, not physical EA accuracy or molecular size consistency.

## Energy ledger and zero conversion

For each basis and **one fixed method/core/reference convention**, compute

\[
\begin{aligned}
A_0 &= E(C)+E(O), & A_O &= E(C)+E(O^-), & A_C &= E(C^-)+E(O),\\
\mathrm{EA}_O &= E(O)-E(O^-), & \mathrm{EA}_C &= E(C)-E(C^-).
\end{aligned}
\]

Retain all totals in Hartree and the algebraic identities `AO-A0=-EA_O`,
`AC-A0=-EA_C`, `AC-AO=EA_O-EA_C`. Convert to physical eV only through
`qscat.units`. The first report gives these as **finite-method predictions**.
Dora2016 p.4, Table3 gives rounded external bindings1.461/1.262eV. Their ordering
and approximate0.199eV separation are plausibility checks, not meV-precision
acceptance data or a basis-convergence proof. More precise isotope/fine-structure
references require their own explicitly sourced convention before scoring errors.

The current molecular pilot's zero is its own energy atR=2.1323, not its fitted
minimum. Report `A0-E_CO(2.1323)` as an **atomic-to-reference offset**; do not label
it De. Only after locating a supported minimum is `De=A0-E_CO(Rmin)` meaningful;
`D0=De-ZPE` additionally needs the physical initial nuclear level. Atoms alone
cannot certify either molecular dissociation depth.

For an eventual common-zero effective model,

\[
V_0(R)=E_{\rm CO,ref}(R)-A_0,\qquad
U_-(R)=V_0(R)+z_{\rm native}(R),\qquad
z_{\rm native}=E_{\rm res}-i\Gamma/2.
\]

If native resonance differences are combined with a separately correlated neutral,
label the result a hybrid reference construction; it is not the absolute anion
energy from either electronic method. Spin-free term energies and lowest-J
experimental thresholds also differ: any correction must include neutral **and**
anion fine-structure shifts. An empirical EA alignment, if later selected, stays
separate from the raw theoretical ledger and cannot hide a wrong neutral depth.

## What the atomic pilot can decide

It supplies missing finite-model endpoints and tests whether their charge/spin/
correlation conventions are usable. Gross or unresolved EA/state discrepancies
trigger one explicit model/reference decision, not automatic reoptimization or
basis escalation. The exact quantitative external gate must await a comparable
spin-free reference; successful job execution is not an absolute-accuracy verdict.

Proposed first release envelope: four sequential atoms, four physical Sadaharu
cores,16GiB container/12GB PySCF memory, minimum24GiB available RAM/10GiB disk;
30minutes per-atom and2hours total science caps, one source/owner implementation
cycle capped at2hours. These are conservative **planning caps**, not measured
runtime or an execution contract. Profile actual cost and retain any first failure.
No fallbacks, repeated molecular QC, extra basis, native inner solve or sweep.

The next molecular work remains separate: stable spin-consistent dissociation
treatment, overlap with accepted near-equilibrium data, bound-anion assignment,
size-consistency/BSSE and finite-R tail residuals. Do not extrapolate the accepted
spline or impose Gamma=0 merely because a width becomes unresolved.

## References

- [Near-equilibrium molecular evidence](../../projects/ukrmol_co/qualification-evidence/neutral-near-equilibrium/README.md).
- [Molecular implementation](../../projects/ukrmol_co/neutral.py), `calculate`.
- [Dora2016 reference note](../../reference/literature/dora-2016-epjd70-197.md),
  pp.4–5, Table3/Fig.1, for asymptotic terms and published bound-region scope.
- [Effective resonance target](co-effective-resonance-target.md), for the native
  threshold, hybrid neutral substitution and independent validation boundary.
