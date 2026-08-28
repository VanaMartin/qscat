"""Round-trip guard for `qscat.core._archive`, the shared `.npz` mechanism
behind `ResonanceLevels.save`/`.load` (`qscat.core.lcp.levels`) and
`ExactResonanceStates.save`/`.load` (`qscat.core.resonance`).

The archives are a CACHE format: a nuclear eigenproblem or a 2-D pole search
is minutes to tens of minutes of sparse factorizations, so a changed key set
silently invalidates saved work rather than raising. These tests build one
small instance of each holder directly (no solver run needed) and assert two
things neither holder's own existing tests makes explicit together: the
`.npz` key set equals the dataclass's own field names, AND every field
survives a save/load round trip unchanged.
"""

from __future__ import annotations

from dataclasses import fields

import numpy as np
from qscat.core.lcp.levels import ResonanceLevels
from qscat.core.resonance import ExactResonanceStates


def _resonance_levels() -> ResonanceLevels:
    n = 3
    return ResonanceLevels(
        energies=np.array([-0.5 + 0.0j, -0.4 - 1e-3j, -0.3 - 2e-3j], dtype=np.complex128),
        widths=np.array([0.0, 2e-3, 4e-3], dtype=np.float64),
        states=np.arange(n * 5, dtype=np.complex128).reshape(n, 5),
        residuals=np.array([1e-10, 2e-10, 3e-10], dtype=np.float64),
        real_weight=np.array([0.99, 0.98, 0.97], dtype=np.float64),
        golden_rule=np.array([-0.5 + 0.0j, -0.4 - 1e-3j, -0.3 - 2e-3j], dtype=np.complex128),
    )


def _exact_resonance_states() -> ExactResonanceStates:
    n = 2
    return ExactResonanceStates(
        energies=np.array([-0.66 - 4e-3j, -0.60 - 1e-2j], dtype=np.complex128),
        widths=np.array([8e-3, 2e-2], dtype=np.float64),
        states=np.arange(n * 7, dtype=np.complex128).reshape(n, 7),
        residual_electronic=np.array([1e-9, 2e-9], dtype=np.float64),
        residual_nuclear=np.array([3e-9, 4e-9], dtype=np.float64),
    )


def test_resonance_levels_key_set_matches_fields_and_round_trips(tmp_path) -> None:
    out = _resonance_levels()
    path = tmp_path / "levels.npz"
    out.save(path)

    with np.load(path) as z:
        assert set(z.files) == {f.name for f in fields(ResonanceLevels)}

    back = ResonanceLevels.load(path)
    for f in fields(ResonanceLevels):
        np.testing.assert_array_equal(getattr(back, f.name), getattr(out, f.name))


def test_exact_resonance_states_key_set_matches_fields_and_round_trips(tmp_path) -> None:
    out = _exact_resonance_states()
    path = tmp_path / "states.npz"
    out.save(path)

    with np.load(path) as z:
        assert set(z.files) == {f.name for f in fields(ExactResonanceStates)}

    back = ExactResonanceStates.load(path)
    for f in fields(ExactResonanceStates):
        np.testing.assert_array_equal(getattr(back, f.name), getattr(out, f.name))
