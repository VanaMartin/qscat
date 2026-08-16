# reference/literature

Local copies of the papers and theses the qModeling models are built from.
**Nothing in this directory except this README is committed** (see `.gitignore`) —
fetch what you need from the URLs below.

Like the rest of `reference/`, these are read-only: they are cited and read, never
imported or built.

## Sources

| File | Reference |
|---|---|
| `vana-2017-thesis.pdf` | M. Váňa, *A model of resonant collisions of electrons with molecules and molecular ions*, doctoral thesis, Charles University, Prague, 2017. <https://dspace.cuni.cz/handle/20.500.11956/92902> |
| `vana-houfek-2017-pra95-022714.pdf` | M. Váňa, K. Houfek, *Time-dependent formulation of the two-dimensional model of resonant electron collisions with diatomic molecules and interpretation of the vibrational excitation cross sections*, Phys. Rev. A **95**, 022714 (2017). <https://doi.org/10.1103/PhysRevA.95.022714> |
| `sieg.pdf` | D. Hvizdoš, M. Váňa, K. Houfek, C. H. Greene, T. N. Rescigno, C. W. McCurdy, R. Čurík, *Dissociative recombination by frame transformation to Siegert pseudostates: a comparison with a numerically solvable model*, Phys. Rev. A **97**, 022704 (2018). <https://arxiv.org/abs/1710.10333> |
| `hvizdos-2016-thesis.pdf` | D. Hvizdoš, *Two-dimensional model of dissociative recombination*, master's thesis, Charles University, Prague, 2016. <https://dspace.cuni.cz/handle/20.500.11956/96080> |

```bash
mkdir -p reference/literature
curl -L -o reference/literature/vana-2017-thesis.pdf \
  "https://dspace.cuni.cz/bitstream/handle/20.500.11956/92902/140060325.pdf?sequence=1"
```

A plain-text extraction (`*.txt`, greppable) can be regenerated with:

```bash
uv run --with pypdf python -c "
from pypdf import PdfReader
r = PdfReader('reference/literature/vana-2017-thesis.pdf')
open('reference/literature/vana-2017-thesis.txt', 'w').write(
    ''.join(f'\n===== PAGE {i+1} =====\n' + (p.extract_text() or '')
            for i, p in enumerate(r.pages)))
"
```

## What the thesis covers

The published counterpart of `reference/eMoScat`, covering every model in
`qscat.model`: the N₂/NO/F₂-like electron-molecule models (ch. 3) and the H₂⁺
dissociative-recombination model (ch. 4).

- Table 1.1 / Table 1.2 — model constants. The N₂/NO/F₂ values match
  `qscat.model.library` exactly (verified 2026-08-15).
- §1.4.4 Eq. 1.60 — `σ = 4π³/k² |T|²`, the normalization the TI oracle uses.
- §1.5 Eq. 1.63/1.65 — the local complex potential
  `V_res(R) = E_res(R) − (i/2)Γ(R)` and `H_LCP`. qscat's `Vd` is the thesis's
  `E_res`.
- §1.6 Eq. 1.69/1.70 — the discrete state `φ_d` and the projection `Ψ_d(R,t)`.
- Tables 2.1/2.2 — the FEM-DVR-ECS grids. These **differ** from the eMoScat JSON
  decks locked by `validation/diatomic/test_da_grid.py` (e.g. F₂ nuclear
  `nq = 20`, θ = 35° here vs order 12, θ = 25° in the deck). Both are plausible
  discretizations of the same model, and quantities computed on them are
  grid-convergent, so the two are treated as a convergence check on each other
  rather than one superseding the other.
- §3.4 — the boomerang / quasi-bound-state interpretation of the cross-section
  structure.
- Bibliography — the canonical citations for this work.

## What Váňa & Houfek (2017) covers

**The peer-reviewed publication of this work's electron–molecule side** — the
time-dependent formulation of the Houfek/Rescigno/McCurdy 2-D model, and the
published counterpart of the thesis's ch. 3. Prefer it over the thesis when
citing anything it contains.

- §II — the 2-D model and its parameters.
- §IV Eqs. (40)/(41) — the local complex potential,
  `V_res(R) = E_res(R) − (i/2)Γ(R)`, and the statement that **the imaginary part
  is nonzero only where `V0(R) < E_res(R)`**. This is the peer-reviewed source
  for the Γ-support condition; the thesis §1.5 repeats it.
- §VIII — the quasibound-state interpretation: each narrow cross-section peak is
  a quasibound vibrational state of the anion (an eigenstate in `V_res(R)`),
  elastic boomerang maxima sit at roughly those energies while the VE 0→1 maxima
  are displaced, and the NO lifetime estimates (first VE 0→1 peak forming at
  t > 10 000; lowest state > 30 000 a.u.).

## What Hvizdoš (2016) covers

**The first time-independent solution of the H₂⁺ model** — the master's thesis
this repo's `dr_cross_section` descends from, supervised by Houfek.

- §1.2.1 and Table 1.1 — the `e⁻ + H₂⁺` model parameters. All match
  `qscat.model.H2P`, **including µ = 918.076** — a third published source
  agreeing with Váňa 2017 Table 1.2 and Hvizdoš et al. (2018) §II A against
  eMoScat's 918.25.
- §1.3 — FEM-DVR, exterior complex scaling, and the driven-equation solution
  (solved in place of Lippmann–Schwinger), which is what `qscat.core.driven`
  and `qscat.core.dissociation` implement.
- §2.1-2.3 — DR cross sections, convergence tests, and the interpretation of
  their structures via the Rydberg potential curves `V_n(R)`.

## What Hvizdoš et al. (2018) covers

The exact-2D benchmark for H₂⁺ dissociative recombination, and the reference for
the frame-transformation / Siegert-pseudostate approximation tested against it.

- §II A Eqs. (4)-(8) — the H₂⁺ model Hamiltonian and Hamilton's coupling
  potential. Matches `qscat.model.library.H2P` **except** the reduced mass: the
  paper and thesis Table 1.2 both give `M = 918.076` (= mₚ/2), qscat has 918.25.
- §II Table I — a third H₂⁺ FEM-DVR-ECS grid parametrization (`nq = 6`, θ = 20°
  on both coordinates).
- §II — the hard ECS bending-angle bounds for this model: nuclear θ < π/8,
  electronic θ < π/4, or `V(R,r)` diverges.
- Appendix A Eqs. (A3)-(A5) — Siegert **pseudostates**: outgoing-wave boundary
  condition at finite `a`, with the surface-corrected orthogonality relation.
  Distinct from qscat's ECS complex-scaled resonance eigenstates; see the
  terminology table in the resonance-levels spec.
