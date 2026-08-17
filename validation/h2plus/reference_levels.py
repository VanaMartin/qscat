"""Published Born-Oppenheimer level positions `omega_i^j` for the H2+ model.

The Rydberg curves `E_Ry_j(R)` and the vibrational levels they support are the
quasi-bound states the DR cross-section peaks are conventionally assigned to
(M. Vana, doctoral thesis, Charles University 2017, Fig. 4.3 and Table 4.1).
These values are the author's own computed level table for this model,
transcribed here as a golden dataset so the repository has an external anchor
for `rydberg_levels`.

**Frame:** electron energy, `E = e_tot - eps[0]` with `eps[0] = -0.0976049` Ha
the cation `v=0` vibrational threshold -- the same convention the published
cross-section windows use. `LEVELS[j][i]` is the thesis's `omega_i^j`: `j` the
Rydberg curve, `i` the vibrational level within it (Fig. 4.3's caption, p. 64,
states the convention explicitly).

**Why the rows are ragged.** Curve `Ry_0` is a much shallower well than the
rest and supports only 5 distinct bound vibrational levels; past that its
published entries repeat -1.2873228 to the printed precision, which is the
accumulation point rather than further levels. Rows here are truncated to the
distinct levels, so a consumer never compares against a repeated continuum
value. Curves `Ry_1`-`Ry_4` are given to 12 levels, well past the range the
published DR windows cover.

**The residual against these values is a known model correction, not error.**
Agreement is <=4e-6 Ha with a mean of 2.4e-6. Substituting eMoScat's reduced
mass `918.25` for this repository's `918.076` drops the mean to **1.1e-7 Ha**,
a 23x improvement -- so the published table was computed with `918.25` and
essentially the whole residual is that one constant. This repository keeps
`918.076` deliberately: Vana 2017 Table 1.2, Hvizdos 2016 Table 1.1 and Hvizdos
et al. 2018 Sec. II A all give `918.076 = m_p/2`, which eMoScat's deck
contradicts (see `reference/literature/README.md`). The gate is therefore set to
accommodate the correction rather than to chase the deck, and the right reading
of a <=4e-6 Ha residual here is "reproduces the published levels to ~1e-7 Ha at
matched constants", not "differs by 4e-6".

**This is a position anchor, not a cross-section anchor**, and that distinction
is the point. Pointwise `sigma_DR(E)` agreement is a badly conditioned test of
this model -- the resonances are ~2e-5 Ha wide and the published sweep's own
1e-5 Ha sampling leaves most of them unresolved, so a few-uHa position
difference moves sigma by tens of percent. These levels test the same physics
where it is well conditioned. See docs/physics/h2plus-dr.md.
"""

from __future__ import annotations

__all__ = ["EPS0", "LEVELS"]

# Cation v=0 vibrational threshold (Ha) -- the zero of the electron-energy frame.
EPS0 = -0.0976049

LEVELS: tuple[tuple[float, ...], ...] = (
    # Ry_0 -- 5 distinct levels; the published row then repeats its accumulation point.
    (-1.291420887, -1.289030648, -1.287807313, -1.287341647, -1.287322820),
    # Ry_1
    (
        -0.063815746,
        -0.055327952,
        -0.047711287,
        -0.041159597,
        -0.037588895,
        -0.036104099,
        -0.034672673,
        -0.032841958,
        -0.031338706,
        -0.029873085,
        -0.028677777,
        -0.027841636,
    ),
    # Ry_2
    (
        -0.034494591,
        -0.025125141,
        -0.016351116,
        -0.008194319,
        -0.000678205,
        0.006173677,
        0.012335044,
        0.017765761,
        0.022412739,
        0.026310312,
        0.029702414,
        0.032780401,
    ),
    # Ry_3
    (
        -0.021612490,
        -0.012012692,
        -0.002963487,
        0.005527095,
        0.013450896,
        0.020800405,
        0.027568702,
        0.033747663,
        0.039326335,
        0.044298293,
        0.048687224,
        0.052562073,
    ),
    # Ry_4
    (
        -0.014806782,
        -0.005116457,
        0.004038586,
        0.012654363,
        0.020726921,
        0.028252736,
        0.035228773,
        0.041651823,
        0.047517732,
        0.052823422,
        0.057575344,
        0.061798852,
    ),
)
