# O₂ curves from Alt & Houfek (2021), Fig. 2

Source: V. Alt, K. Houfek, *Resonant collisions of electrons with O₂ via the
lowest-lying ²Π_g state of O₂⁻*, Phys. Rev. A **103**, 032829 (2021), Fig. 2
(p. 032829-3): the ³Σ_g⁻ curve of O₂ and the ²Π_g curve of O₂⁻ (CASSCF/MRCI,
aug-cc-pVQZ, the anion curve shifted to the experimental EA(O) = 1.4611 eV —
their §III A and Table I), the real part of the resonance energy below the
crossing, and the local width Γ(R) (plotted ×2) from their R-matrix fits.

**These are not digitised by eye.** The published PDF embeds the figure as
vector graphics — every curve is a thin filled polyline outline of the
original stroke — so `validation/factory/extract_fig2.py` reads the paths
directly, calibrates the axes from the ten `R` tick marks (1.5–6.0 bohr) and
the eight `E` tick marks (−6–+1 eV; residuals 3e-4 bohr and 1.4e-3 eV), samples
each outline segment every 0.2 pt, and takes the median of the outline's two
edges per 0.6-pt bin (in `x` where the curve is flat, in `y` inside each steep
run). The outline's half-width is ~0.3 pt, so the **vertical precision is about
0.02 eV** on the flat parts and the figure's own polyline resolution (~0.09 bohr
between vertices) elsewhere. Re-run the extractor with the paper's PDF at
`reference/literature/alt-houfek-2021-pra103-032829.pdf` (gitignored).

| file | curve in Fig. 2 | columns | points |
|---|---|---|---|
| `v0.csv` | ³Σ_g⁻ (black, full), zero at O(³P)+O(³P) | `R_bohr, E_eV` | 469, R 1.74–6.0 |
| `v_ion_bound.csv` | ²Π_g where bound, R ≥ 2.289 (blue, full) | `R_bohr, E_eV` | 279 |
| `e_res_dashed.csv` | Re of the resonance energy, R < 2.289 (green, dashed) | `R_bohr, E_eV` | 104 |
| `gamma_x2.csv` | Γ(R) × 2 (black, dash-dotted); zero beyond 2.41 | `R_bohr, E_eV` | 51 |

Sanity of the extraction: the `v0` minimum is −5.259 eV at 2.25 bohr against
Table I's D₀ = 5.159 eV (calc.) plus the 0.098 eV zero-point energy; `v0(6 bohr)`
= −0.01 eV; the anion/neutral crossing is at 2.289 bohr; Γ falls from 0.67 eV at
1.8 bohr to zero at 2.3–2.4 bohr. The ladder of the extracted `v0` gives
`G(1)−G(0)` ≈ 1607 cm⁻¹, ~3 % above the spectroscopic 1551 cm⁻¹ — that is a
property of the paper's MRCI curve (and of the extraction floor near the
minimum), and the factory's T0 check therefore compares against the curve's
own ladder, not the spectroscopic constant.

`v_ion_bound` at 6 bohr is still 0.2 eV below the O(³P)+O⁻(²P) asymptote the
figure marks (−1.461 eV); the factory's electron-affinity constraint uses
Table I's EA(O), not this curve's last point.

## Fig. 1 — the spin–orbit splitting of the ²Π_g curve

`so_split.csv` (`R_bohr, split_meV`, 426 points, 1–16 bohr) is the paper's
own Gaussian fit of the total spin–orbit splitting Δ_SO(R) of the O₂⁻ ²Π_g
curve (Fig. 1, p. 032829-3; 19.7 meV at 2.1 bohr, a 12.2 meV plateau beyond
9 bohr), vector-extracted by `validation/factory/extract_fig1.py` (median of
the outline's edges per 0.4-pt bin, ~0.05 meV). It is what the paper used to
build the ²Π_{1/2} and ²Π_{3/2} curves at ∓Δ_SO(R)/2 around ²Π_g (§III A),
and what `targets/o2.py::o2_target(so=±1)` uses for the same purpose. Only
the fit curve is taken, not the MOLPRO points it was fitted to.

## Fig. 5 — the paper's VE cross sections (theory only)

`fig5_ve_0{v}_{nrm,lcp}.csv` (`E_eV, sigma_a0^2`), v' = 0…5, are the
paper's own nonlocal-model (NRM, blue) and local-complex-potential (LCP,
dashed green) vibrational-excitation cross sections from Fig. 5
(p. 032829-6), extracted the same way by `validation/factory/extract_fig5.py`:
each of the six panels is calibrated from its tick rectangles (the labels
are glyph outlines, so their values are fixed in the extractor after one
reading of the rendered page), and the curve is the outline's **upper
envelope** per 0.15-pt bin minus half the stroke width — a comb of meV-wide
peaks is exactly where a median-of-edges centreline fails, while the
envelope keeps every peak's height to 0.25 pt (0.3–0.7 % of a panel's range)
at an energy resolution of one bin, 0.9–1.4 meV. Their only use is the
theory-vs-theory overlay against the factory model's spin–orbit-resolved
exact 2-D cross section (⅓ ²Π_{1/2} + ⅓ ²Π_{3/2}, the same composition as
the paper's curves); nothing is fitted to them, and the experimental traces
of Figs. 7–9 are not extracted.
