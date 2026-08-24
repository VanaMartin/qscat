"""Task 5: the six C5/C4 anchors, run through the exact 2-D solver, compared
three ways (`validation.n2.exact2d.compute_exact2d_results`).

Two families of checks, mirroring `projects/n2_ti_cross_section/test_cross_section.py`
and `validation/n2/cross_section.py`'s own split:

- INTERNAL (model-independent): sigma is real and >=0 at every anchor.
- V4, the GATE: the 4 GATED anchors (VE channels clear of their own
  threshold) must agree with Houfek's independent `CSVE.V00.J00` almost
  exactly -- this is a differential-oracle check of THIS solver against an
  independent implementation of the *same* model, not a loose cross-model
  bound like the LCP's factor-of-3 (`reference.ANCHOR_FACTOR`).

Tolerance derivation (measured, not picked in advance): running
`compute_exact2d_results()` gives, at the 4 GATED anchors,
`|ratio_exact_vs_houfek - 1|` = 1.149e-06 (E=0.2, v'=1), 8.430e-06 (E=0.2,
v'=2), 1.665e-05 (E=0.2, v'=3), 2.672e-04 (E=0.1, v'=1) -- the last of
these, the largest, sets the scale. `GATED_RTOL = 1e-3` sits ~3.7x above
that largest measured deviation: enough headroom for run-to-run solver/BLAS
variation while remaining a real, tight differential-oracle gate, nowhere
close to the LCP's factor-of-3 cross-model band.

V5 (the scientific deliverable, `ratio_lcp_vs_exact`) and V6 (the two
LCP NOTEs re-examined with a real oracle) are not separately gated here --
they are not reported here -- but this file DOES check the qualitative V5 claim that
motivates the whole sub-project: once gated, the exact solver is at least
as trustworthy a comparison point as Houfek itself, so the LCP should never
look *better* against Houfek than the exact solver does.
"""

from __future__ import annotations

import numpy as np

from validation.n2.exact2d import GATED_RTOL, compute_exact2d_results

# GATED_RTOL is defined once in `validation.n2.exact2d` (see its module-level
# docstring for the derivation from the measured deviations) and imported
# here and by `validation/n2/experiment.py` (Group E) so the two cannot
# drift apart.


def test_sigma_real_and_nonnegative() -> None:
    for r in compute_exact2d_results():
        assert np.isfinite(r.sigma_exact)
        assert r.sigma_exact >= 0.0


def test_gated_anchors_match_houfek_within_derived_tolerance() -> None:
    """V4, the gate: the exact 2-D solver vs. Houfek's independent
    `CSVE.V00.J00`, at the 4 anchors clear of their own vibrational
    threshold. See module docstring for how `GATED_RTOL` was derived from
    the measured deviations, not chosen in advance.
    """
    results = compute_exact2d_results()
    gated = [r for r in results if r.gated]
    assert len(gated) == 4  # (0.2,1), (0.2,2), (0.2,3), (0.1,1)

    for r in gated:
        dev = abs(r.ratio_exact_vs_houfek - 1.0)
        assert dev < GATED_RTOL, (
            f"anchor (E={r.energy_ha}, v'={r.channel}): exact/houfek ratio "
            f"{r.ratio_exact_vs_houfek:.6f} deviates {dev:.3e}, outside "
            f"GATED_RTOL={GATED_RTOL}"
        )


def test_exact_at_least_as_close_to_houfek_as_lcp_at_gated_anchors() -> None:
    """The exact solver is the ORACLE; it must never agree with Houfek
    *worse* than the LCP approximation it is meant to test does. Measured:
    LCP deviations at the 4 gated anchors are 0.555, 0.226, 0.173, 0.010
    (E=0.2 v'=1,2,3 and E=0.1 v'=1) vs. the exact solver's 1.1e-06, 8.4e-06,
    1.7e-05, 2.7e-04 -- the exact solver is 3-5 orders of magnitude closer
    at every gated anchor, so this holds comfortably. If this ever regressed
    (exact worse than LCP at some gated anchor), the fix is NOT to widen
    this assertion -- per the task brief, that would need an `xfail` with a
    written investigation, not a silently loosened bound.
    """
    for r in compute_exact2d_results():
        if not r.gated:
            continue
        dev_exact = abs(r.ratio_exact_vs_houfek - 1.0)
        dev_lcp = abs(r.ratio_lcp_vs_houfek - 1.0)
        assert dev_exact <= dev_lcp, (
            f"anchor (E={r.energy_ha}, v'={r.channel}): exact deviates "
            f"{dev_exact:.3e} from Houfek, WORSE than LCP's {dev_lcp:.3e}"
        )


def test_v6_exact_model_closes_the_two_lcp_notes() -> None:
    """V6 (a measurement): at the LCP's own two documented structural
    failures -- elastic (E=0.2, v'=0, missing background scattering) and
    near-threshold (E=0.02, v'=1, wrong non-Wigner threshold law) -- does
    the exact 2-D model (real background scattering, an emergent
    energy-dependent effective width) actually close the gap to Houfek?

    Measured: elastic ratio_exact_vs_houfek=1.0000 vs. the LCP's
    ratio_lcp_vs_houfek=0.0401 (off by ~25x); near-threshold
    ratio_exact_vs_houfek=0.9998 vs. the LCP's ratio_lcp_vs_houfek=8133
    (off by ~4 orders of magnitude). Both close dramatically -- this is the
    predicted outcome, confirmed, not assumed; the assertions below record
    that finding rather than merely asserting sigma is finite.
    """
    ungated = [r for r in compute_exact2d_results() if not r.gated]
    elastic = next(r for r in ungated if r.channel == 0)
    near_threshold = next(r for r in ungated if r.channel != 0)

    for r in (elastic, near_threshold):
        dev_exact = abs(r.ratio_exact_vs_houfek - 1.0)
        dev_lcp = abs(r.ratio_lcp_vs_houfek - 1.0)
        assert dev_exact < dev_lcp, (
            f"NOTE anchor (E={r.energy_ha}, v'={r.channel}): expected the exact "
            f"solver to close the LCP's documented gap to Houfek, but exact "
            f"deviation {dev_exact:.3e} >= LCP deviation {dev_lcp:.3e}"
        )
