# reference/literature

The papers and theses qModeling is built on, and the **reference notes** that
carry what the repository actually relies on from each.

**The PDFs are local copies only** — `*.pdf` and the `*.txt` extractions are
gitignored, so they do not travel with a clone. **The `*.md` notes are tracked**,
because they are what the repository cites. A note carries every published fact
the repo depends on, each anchored to a page plus an equation, table or figure,
so a reader can find it in the source in one jump.

Like the rest of `reference/`, the sources are read, never imported or built.

To write or update a note, use the `mastering-references` skill
(`.claude/skills/mastering-references/`).

## Sources

| Note | Source | Role here |
|---|---|---|
| [`rescigno-2000-pra62-032706.md`](rescigno-2000-pra62-032706.md) | T. N. Rescigno, C. W. McCurdy, *Numerical grid methods for quantum-mechanical scattering problems*, Phys. Rev. A **62**, 032706 (2000). [DOI](https://doi.org/10.1103/PhysRevA.62.032706) | The FEM-DVR + ECS method itself — the foundation of `qscat.dvr` and `qscat.ecs`. |
| [`tannor-1993-jcp98-3884.md`](tannor-1993-jcp98-3884.md) | D. J. Tannor, D. E. Weeks, *Wave packet correlation function formulation of scattering theory: the quantum analog of classical S-matrix theory*, J. Chem. Phys. **98**, 3884 (1993). [DOI](https://doi.org/10.1063/1.464016) | **The Tannor-Weeks method** — propagate-and-correlate energy extraction, the `eta_incident`/`eta_outgoing` normalization factors, and the "outgoing normalization isn't 1" fact behind the elastic free-reference subtraction. |
| [`houfek-2006-pra73-032721.md`](houfek-2006-pra73-032721.md) | K. Houfek, T. N. Rescigno, C. W. McCurdy, *Numerically solvable model for resonant collisions of electrons with diatomic molecules*, Phys. Rev. A **73**, 032721 (2006). [DOI](https://doi.org/10.1103/PhysRevA.73.032721) | **The 2-D model this whole repository implements.** N₂/NO constants. |
| [`houfek-2008-pra77-012710.md`](houfek-2008-pra77-012710.md) | K. Houfek, T. N. Rescigno, C. W. McCurdy, *Probing the nonlocal approximation to resonant collisions of electrons with diatomic molecules*, Phys. Rev. A **77**, 012710 (2008). [DOI](https://doi.org/10.1103/PhysRevA.77.012710) | The discrete-state choice, nonlocal vs LCP — and **F₂'s published constants**. |
| [`vana-houfek-2017-pra95-022714.md`](vana-houfek-2017-pra95-022714.md) | M. Váňa, K. Houfek, *Time-dependent formulation of the two-dimensional model…*, Phys. Rev. A **95**, 022714 (2017). [DOI](https://doi.org/10.1103/PhysRevA.95.022714) | **This work's own publication.** The LCP curve, the Γ-support condition, the quasibound-state interpretation. |
| [`vana-2017-thesis.md`](vana-2017-thesis.md) | M. Váňa, *A model of resonant collisions of electrons with molecules and molecular ions*, doctoral thesis, Charles University, Prague, 2017. [dspace](https://dspace.cuni.cz/handle/20.500.11956/92902) | The broadest single source: all four models, grids, the discrete-state projection. |
| [`hvizdos-2016-thesis.md`](hvizdos-2016-thesis.md) | D. Hvizdoš, *Two-dimensional model of dissociative recombination*, master's thesis, Charles University, Prague, 2016. [dspace](https://dspace.cuni.cz/handle/20.500.11956/96080) | **The first time-independent solution of the H₂⁺ model**, which `dr_cross_section` descends from. |
| [`hvizdos-2018-pra97-022704.md`](hvizdos-2018-pra97-022704.md) | D. Hvizdoš, M. Váňa, K. Houfek, C. H. Greene, T. N. Rescigno, C. W. McCurdy, R. Čurík, *Dissociative recombination by frame transformation to Siegert pseudostates…*, Phys. Rev. A **97**, 022704 (2018). [arXiv](https://arxiv.org/abs/1710.10333) | The Siegert-**pseudostate** construction qModeling's ECS resonance eigenstates are explicitly *not*. H₂⁺ parameters and ECS angle bounds. |
| [`vandijk-2007-pre75-036707.md`](vandijk-2007-pre75-036707.md) | W. van Dijk, F. M. Toyama, *Accurate numerical solutions of the time-dependent Schrödinger equation*, Phys. Rev. E **75**, 036707 (2007). [DOI](https://doi.org/10.1103/PhysRevE.75.036707) | The origin of the order-N diagonal-Padé time propagator (`qscat.evolution.make_pade_stepper`); order 1 is Crank-Nicolson. |
| [`mccurdy-1991-pra43-5980.md`](mccurdy-1991-pra43-5980.md) | C. W. McCurdy, C. K. Stroud, M. K. Wisinski, *Solving the time-dependent Schrödinger equation using complex-coordinate contours*, Phys. Rev. A **43**, 5980 (1991). [DOI](https://doi.org/10.1103/PhysRevA.43.5980) | Origin of ECS + Crank-Nicolson propagation on complex-symmetric matrices — the ancestor of `qscat.evolution`. Carries a third, distinct ECS angle bound. |
| [`mccurdy-1991-cpc63-323.md`](mccurdy-1991-cpc63-323.md) | C. W. McCurdy, C. K. Stroud, *Eliminating wavepacket reflection from grid boundaries using complex coordinate contours*, Comput. Phys. Commun. **63**, 323 (1991). [DOI](https://doi.org/10.1016/0010-4655(91)90259-N) | Precursor: where ECS-as-absorber for time propagation was first demonstrated. Background, not the repo's cited source. |
| [`mccurdy-2004-jpb37-r137.md`](mccurdy-2004-jpb37-r137.md) | C. W. McCurdy, M. Baertschy, T. N. Rescigno, *Solving the three-body Coulomb breakup problem using exterior complex scaling*, J. Phys. B **37**, R137 (2004). [DOI](https://doi.org/10.1088/0953-4075/37/17/R01) | The general ECS contour formalism our sharp contour specializes from, and the Wronskian amplitude extraction behind `td_extractors.Flux`. |
| [`gertitschke-1993-jpb26-2927.md`](gertitschke-1993-jpb26-2927.md) | P. L. Gertitschke, W. Domcke, *Systematically improved local complex potential approximation for the dynamics of electron-molecule collision complexes*, J. Phys. B **26**, 2927 (1993). [DOI](https://doi.org/10.1088/0953-4075/26/17/024) | Corroborating evidence on where LCP holds and fails. **Not** the Phys. Rep. 208 review — that is `domcke-1991-physrep208-97.md`. |
| [`domcke-1991-physrep208-97.md`](domcke-1991-physrep208-97.md) | W. Domcke, *Theory of resonance and threshold effects in electron-molecule collisions: the projection-operator approach*, Phys. Rep. **208**, 97 (1991). PII `0370-1573(91)90125-6` | **The canonical nonlocal resonance model.** The nuclear equation `qscat.core.nrm` solves, the LCP limit derived from it, and the Eq. (4.14) coupling PRA 77 disputes. |
| [`gertitschke-1993-pra47-1031.md`](gertitschke-1993-pra47-1031.md) | P. L. Gertitschke, W. Domcke, *Time-dependent wave-packet description of dissociative electron attachment*, Phys. Rev. A **47**, 1031 (1993). [DOI](https://doi.org/10.1103/PhysRevA.47.1031) | The time-dependent nonlocal treatment of DA. Nothing cites it yet — ingested for the planned TD-NRM sub-project. Quantifies the LCP's e+H₂ failure (14×/23×) and its wave-packet-splitting mechanism. |
| [`formanek-2010-aip1281-667.md`](formanek-2010-aip1281-667.md) | M. Formánek, M. Váňa, K. Houfek, *Comparison of the Chebyshev Method and the Generalized Crank-Nicholson Method for Time Propagation in Quantum Mechanics*, AIP Conf. Proc. **1281**, 667 (2010). [DOI](https://doi.org/10.1063/1.3498565) | Background: the authors' own efficiency comparison, context for why qscat propagates with Padé rather than Chebyshev. |

## Fetching the sources

The notes stand on their own; fetch a PDF only when you need the full text.

```bash
mkdir -p reference/literature
curl -L -o reference/literature/vana-2017-thesis.pdf \
  "https://dspace.cuni.cz/bitstream/handle/20.500.11956/92902/140060325.pdf?sequence=1"
curl -L -o reference/literature/hvizdos-2016-thesis.pdf \
  "https://dspace.cuni.cz/bitstream/handle/20.500.11956/96080/150040279.pdf"
```

The APS articles are paywalled; obtain them through an institutional
subscription and save them under the note's stem with a `.pdf` extension.

To regenerate a text extraction for searching, see
`.claude/skills/mastering-references/references/extracting.md`.

## Pagination caveats

Locators in each note use the **printed** page, which is not always the
extractor's page index. Each note states its own mapping; two are worth knowing
up front:

- **The theses** carry unnumbered front matter — Hvizdoš 2016 runs at
  extractor = printed + 8, Váňa 2017 at + 8.
- **`hvizdos-2018-pra97-022704.pdf` is a preprint copy**, self-paginated 1–26
  rather than the journal's `022704-N`. Its note's locators follow the preprint;
  a reader holding the published version must translate.
- **`mccurdy-1991-pra43-5980.pdf` and `gertitschke-1993-pra47-1031.pdf` are
  scans with no text layer.** Their notes were written from rendered page
  images; there is no `.txt` extraction to grep.
- **`domcke-1991-physrep208-97.pdf` has an OCR text layer that garbles
  equations.** The `.txt` is usable for locating an equation by number but not
  for reading it; every equation in that note was verified against a rendered
  page image. It is an offprint, so extractor page N = printed page `96 + N`.

## Parity with the code

Every note states, for each published constant, whether it was checked against
the repository and what the check found. Three findings so far came out of that
discipline:

- **H₂⁺ reduced mass.** Váňa 2017 Table 1.2, Hvizdoš 2016 Table 1.1 and
  Hvizdoš et al. 2018 §II A all give **918.076** = mₚ/2; eMoScat's deck carried
  918.25. The repo was corrected to 918.076.
- **F₂'s provenance.** Its constants are published in Houfek et al. (2008)
  Table I, not merely transcribed from the eMoScat deck.
- **Grid decks.** Váňa 2017 Tables 2.1/2.2 differ from the eMoScat decks locked
  by `validation/diatomic/test_da_grid.py`. Both are plausible discretizations
  of the same model and quantities on them are grid-convergent, so the two serve
  as a convergence check on each other.

## Known gaps

Conventions the repository leans on that no note here anchors to a page.

- **The c-product.** `qscat.linalg.c_product` — the bilinear, non-conjugated
  pairing `sum_i a_i b_i` — is used for every T-matrix projection, correlation
  function, state normalization and overlap in `qscat.core`. What the notes DO
  anchor is that `H` is complex **symmetric** under ECS (McCurdy et al. 1991,
  Phys. Rev. A 43, p. 5985, Eqs. (21)-(32); the CPC 63 note carries the explicit
  "conjugate -> transpose" for the propagator). The step from there to the
  *inner product* is inferred, not cited: no note says in so many words that the
  pairing is `int psi^2` without a conjugate.

  The inference is sound and independently corroborated — the Hermitian
  convention gives a negative sigma on N2's S-matrix, and eMoScat's `cblas_zdotc`
  agrees with the c-product to 3.4e-12 only because it zeroes every channel
  function on the scaled tail. But it is the convention the repo relies on most
  heavily with the weakest citation, and a reader cannot follow it to a page.

  **Wanted:** a source that states the bilinear pairing for complex-scaled /
  non-Hermitian quantum mechanics directly. Moiseyev's *Non-Hermitian Quantum
  Mechanics* (or his Phys. Rep. **302**, 212 (1998) review) is the obvious
  candidate; drop the PDF in and write the note per the `mastering-references`
  skill. Do not write the note without the source — every locator has to be
  checked against a page, and this file exists because that discipline has
  already caught three real discrepancies.
