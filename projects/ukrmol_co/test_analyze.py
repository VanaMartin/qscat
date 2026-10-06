"""Output-format regression checks; numerical engine validation uses saved runs."""

import numpy as np

from projects.ukrmol_co.analyze import ground_cross_sections


def test_ground_cross_sections_joins_columns_and_excludes_excited_initial_states(tmp_path):
    """CC output splits final-state columns across differently sized blocks."""
    path = tmp_path / "xsec"
    path.write_text(
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 1\n"
        "# I E TOTAL 1 2\n"
        "1 0.1 2.0 2.0 0.0\n2 0.2 3.0 3.0 0.0\n\n"
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 1\n"
        "# I E TOTAL 7\n"
        "1 0.1 0.5\n2 0.2 0.7\n"
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 2\n"
        "# I E TOTAL 1\n"
        "1 0.1 9.0 9.0\n2 0.2 8.0 8.0\n"
    )
    np.testing.assert_allclose(
        ground_cross_sections(path),
        [[1, 0.1, 2, 2, 0, 0.5], [2, 0.2, 3, 3, 0, 0.7]],
        atol=0,
        rtol=0,
    )
