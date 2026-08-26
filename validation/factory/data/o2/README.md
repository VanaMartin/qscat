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
