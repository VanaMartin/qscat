"""The C5 Houfek anchor gate, as a pytest test.

Replaces `projects/n2_ti_cross_section/test_cross_section.py`'s
`test_houfek_anchor_agreement`, which hardcoded the two DOCUMENTED-LIMITED
coordinates; here the GATED / DOCUMENTED-LIMITED split comes from
`validation.n2.cross_section.classify`'s general `(energy, channel)` rule,
so the gate and the `python -m validation.n2.experiment` harness cannot
disagree. `_build_system` is lru_cached, so the ~7 s vres_on_grid walk is
paid once per process however many tests read it.
"""

from __future__ import annotations

from validation.n2 import reference
from validation.n2.cross_section import compute_anchor_results


def test_gated_anchors_within_factor_band() -> None:
    results = compute_anchor_results()
    assert len(results) == 6
    gated = [r for r in results if r.gated]
    assert len(gated) == 4  # the general rule must still gate exactly 4 of 6
    for r in gated:
        assert 1.0 / reference.ANCHOR_FACTOR <= r.ratio <= reference.ANCHOR_FACTOR, (
            f"anchor (E={r.energy_ha}, v'={r.channel}) ratio {r.ratio:.3f} "
            f"outside factor-of-{reference.ANCHOR_FACTOR} band"
        )


def test_documented_limited_anchors_carry_a_mechanism() -> None:
    for r in compute_anchor_results():
        if not r.gated:
            assert r.mechanism  # never silently excluded
        print(
            f"E={r.energy_ha:.4f} Ha, v'={r.channel}: computed={r.sigma_computed:.4e} "
            f"houfek={r.sigma_houfek:.4e} ratio={r.ratio:.3f} "
            f"{'GATED' if r.gated else '[' + r.mechanism + ']'}"
        )
