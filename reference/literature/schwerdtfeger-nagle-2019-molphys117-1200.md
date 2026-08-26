# Schwerdtfeger & Nagle, Mol. Phys. 117, 1200 (2019) — 2018 Table of static dipole polarizabilities of the neutral elements (2023 update)

**Source:** `reference/literature/schwerdtfeger-nagle-2019-molphys117-1200.pdf` (gitignored) · <https://doi.org/10.1080/00268976.2018.1535143> · the file is the authors' maintained *2023 Table of Polarizabilities* (last update 17 July 2023), <https://ctcp.massey.ac.nz/2023Tablepol.pdf>
**Pagination:** the PDF's own page numbers, 1 … (p. 1 the periodic-table figure, p. 3 onward Table 1 by increasing `Z`).

## Why this repository cares

The potential factory ends every anion curve at the atom + ion energy `−EA` and, at finite `R`, at the ion–atom polarisation tail `−α_d/(2R⁴)` (`projects/potential_factory/target.py::polarisation_tail`). That tail needs the neutral partner atom's static dipole polarisability `α_d`, an ATOMIC constant with no place in a molecular figure; this table is the recommended-value source for it, with a stated uncertainty. Without it the O₂ target's long-range form would be a guess.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| Recommended `α_d(O, ³P) = 5.3 ± 0.2` a.u. | p. 6, Table 1, `Z = 8` (also p. 1, Figure 1) | `validation/factory/targets/o2.py` (`ALPHA_D_O`) |
| Recommended `α_d(N, ⁴S) = 7.4 ± 0.2` a.u. | p. 5, Table 1, `Z = 7` | not yet used (an N₂-like target's tail) |
| Recommended `α_d(F, ²P) = 3.74 ± 0.08` a.u. | p. 6, Table 1, `Z = 9` | not yet used (an F₂-like target's tail) |
| `1 a.u. = 0.1481847113 Å³` | p. 2, Table 1 caption | unit sanity only |

## Equations

None. The table lists values; the tail form `−α_d/(2R⁴)` the repository attaches to them is the textbook charge–induced-dipole interaction, not a result of this source.

## Parameters and numeric values

| atom | `α_d` (a.u.) | state | locator |
|---|---|---|---|
| N | 7.4 ± 0.2 | ⁴S | p. 5, Table 1 |
| O | 5.3 ± 0.2 | ³P | p. 6, Table 1 |
| F | 3.74 ± 0.08 | ²P | p. 6, Table 1 |

Checked against the repo (`grep -n "ALPHA_D_O" validation/factory/targets/o2.py`): `ALPHA_D_O = 5.3` — **matches Table 1 exactly (verified 2026-08-25)**. N and F are not in the code.

## Findings and limits

- The recommended values are `M_L`-averaged scalar polarisabilities unless a state symmetry is given (p. 2, Table 1 caption); the O value's ±0.2 spread covers CCSD(T) (5.24 ± 0.04, 5.21 ± 0.07) and the experimental 5.2 ± 0.4 (p. 6).
- The tail's uncertainty from `α_d` alone is ~4 % of `−α_d/(2R⁴)`: at O₂'s `R_inf = 14` bohr that term is 0.07 mHa, so the uncertainty is far below the factory's 1 mHa gate — `R_inf` is what decides the tail's accuracy, not `α_d`.

## Terminology map

| table | qModeling |
|---|---|
| `α_D` (static scalar dipole polarisability, a.u.) | `alpha_d` in `polarisation_tail(alpha_d)` |

## Not used here

Every other element, the per-reference entries behind each recommendation, hyperpolarisabilities, and the method abbreviations of Table 1's caption.
