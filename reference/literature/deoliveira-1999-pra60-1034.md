# de Oliveira et al., Phys. Rev. A 60, 1034 (1999) — atomic electron affinities

**Source:** `reference/literature/deoliveira-1999-pra60-1034.pdf` (gitignored) ·
[DOI](https://doi.org/10.1103/PhysRevA.60.1034) ·
[published-edition copy](https://repository.uantwerpen.be/docman/irua/873320/7956.pdf)
**Pagination:** printed 1034–1045; extraction page N = 1033 + N. The headers
on extraction pages 3/4 read 1036/1037; the final page reads 1045. Table I and
the oxygen discussion were also checked in rendered published pages, not
translated from the differently paginated arXiv preprint.
PDF SHA256: `b7c22b61a8b03201170dcac0acffa50717f4f391c9144c09037d07961d08f583`.

## Why this repository cares

The oxygen-affinity diagnosis must distinguish finite valence-basis CCSD(T)
from a complete ground-threshold prediction. This source quantifies missing
basis and higher-excitation effects rather than attributing their sum to an
incorrect occupation or orbital conversion. Its corrections are external
benchmarks, not calculated corrections adopted by this repository.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| ROHF references are used unless stated otherwise | p. 1036, II B | Atomic-method comparison |
| Valence CCSD(T), core correlation, spin-orbit, scalar relativistic and higher-excitation terms are distinct contributions | p. 1036, Table I and II B | Oxygen reference error budget |
| For O, extrapolation beyond AV6Z still adds about 0.016 eV | p. 1037, oxygen paragraph, left-column bottom continuing at right-column top | Reject assuming a 5Z absolute result is converged |
| O's roughly 0.012-eV higher-excitation correction is mostly quadruples and higher; CCSDT−CCSD(T) alone is about 0.002 eV | p. 1037, right-column top | Correlation-control choice |
| FCI for O is performed in AVDZ; CCSDT triples calibration uses AVQZ | p. 1036, II B, right-column top | Limits of a small-basis FCI control |

## Equations

None implemented from this source. Its extrapolation expressions on p. 1036,
II B are not used to produce a repository CBS number.

## Parameters and numeric values

All entries below are eV from **p. 1036, Table I, O row**. The CCSD(T) column
is a *correlation contribution*, not the entire affinity.

| Contribution | Published value |
|---|---:|
| SCF limit | −0.53902 |
| Valence CCSD(T) correlation limit | +1.99391 |
| Core correlation | +0.00173 |
| Spin-orbit | −0.00222 |
| Darwin + mass-velocity | −0.00588 |
| FCI/higher-excitation calibration | +0.01223 |
| Best calculated total | 1.46075 |
| Then-experimental value | 1.461122(3) |

Parity check: reconstructed saved 5Z O/O− logs give SCF −0.538835931 eV
and CCSD(T) correlation increment +1.965442390 eV. Neither is a published
limit. The difference to the published valence correlation limit is
28.467610 meV; this is a cross-method/basis benchmark comparison, **not** a
measured correction or a CBS error bound. No table correction is installed in
the atomic driver. Published rounded components sum to 1.46075 eV.

## Findings and limits

O has exceptionally slow valence-basis convergence and appreciable excitations
beyond perturbative triples (p. 1037, oxygen paragraph). The higher-excitation
estimate uses different finite bases for CCSDT and FCI (p. 1036, II B), so it
does not establish a rigorous uncertainty bound for our UCCSD(T) convention.
The experimental number is the value cited in 1999, not a newly verified
isotope/field-resolved experimental recommendation.

## Terminology map

| Paper | Repository |
|---|---|
| AVnZ | aug-cc-pVnZ |
| n-particle / FCI correction | correlation beyond CCSD(T), not basis augmentation |
| DMV | scalar Darwin + mass-velocity correction, not spin-orbit |

## Not used here

DFT comparisons, second-row recommendations and the paper's extrapolated
numbers as substitutes for our own calculations. No empirical EA alignment.
