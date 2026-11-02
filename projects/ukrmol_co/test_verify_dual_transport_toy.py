"""Independent scalar-oracle controls; no integration or campaign data reads."""

import pytest

from projects.ukrmol_co.verify_dual_transport_toy import measurements, scalar_oracle


def test_scalar_critical_angle_and_negative_control():
    row = scalar_oracle(45, 20, 1e-6)
    assert abs(row[1]) == pytest.approx(1e-6, abs=1e-20, rel=1e-12)
    result = measurements(row, row, 1e-12)
    assert result["scaled_error"] == pytest.approx(0, abs=1e-15, rel=0)
    assert result["s_identity_error"] < 1e-6 < result["flux_defect"] * 10


def test_independent_scalar_flux_detects_open_normalization_error():
    result = measurements([complex(2), complex(0)], [complex(1), complex(0)], 1e-12)
    assert result["flux_defect"] == pytest.approx(3, abs=1e-15, rel=0)
    assert result["s_identity_error"] == pytest.approx(0, abs=1e-15, rel=0)
