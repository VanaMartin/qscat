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
| `sieg.pdf` | D. Hvizdoš, M. Váňa, K. Houfek, C. H. Greene, T. N. Rescigno, C. W. McCurdy, R. Čurík, *Dissociative recombination by frame transformation to Siegert pseudostates: a comparison with a numerically solvable model*, Phys. Rev. A **97**, 022704 (2018). <https://arxiv.org/abs/1710.10333> |

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
