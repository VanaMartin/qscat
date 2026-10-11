"""Guard the model-dependent DENPROP ordering before a costly native solve."""

import pytest

from projects.ukrmol_co.sparse_pilot import target_map


def test_contraction_map_follows_denprop_ids_and_checks_model_energies():
    # Interleaved native energy order: it must not be replaced by sector order.
    labels = [
        ("triplet", "B2"),
        ("singlet", "A1"),
        ("triplet", "A1"),
        ("singlet", "B1"),
        ("singlet", "B2"),
        ("triplet", "B1"),
        ("singlet", "A2"),
        ("triplet", "A2"),
    ]
    states = [
        {"spin": spin, "irrep": irrep, "energy_hartree": -113 + i * 0.01}
        for i, (spin, irrep) in enumerate(labels)
    ]
    prop = "\n".join(
        f"5 {i + 1} 0 0 0 1 0 0 {s['energy_hartree']:.12f} "
        f"State No. {i + 1} {s['spin'].title()} {s['irrep']}"
        for i, s in enumerate(states)
    )
    assert target_map(prop, states, roots=1) == [2, 3, 4, 6, 5, 1, 7, 8]
    with pytest.raises(ValueError, match="Incomplete or reordered"):
        target_map(prop.replace("5 1 ", "5 2 "), states, roots=1)
    with pytest.raises(ValueError, match="Incomplete retained target"):
        target_map(prop.replace("Triplet B2", "Triplet A2"), states, roots=1)
    changed = [dict(s) for s in states]
    changed[0]["energy_hartree"] += 1e-5
    with pytest.raises(AssertionError):
        target_map(prop, changed, roots=1)
