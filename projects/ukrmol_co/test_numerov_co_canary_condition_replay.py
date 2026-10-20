"""A rounding observation cannot hide an operand or residual discrepancy."""

from projects.ukrmol_co.numerov_co_canary_condition_replay import condition_drift, scientific_record


def record(rcond=0.02, matrix_hash="abc"):
    return dict(
        samples=[dict(reciprocal_condition_1=rcond, a_sha256=matrix_hash)], calls=[], checks={}
    )


def test_condition_roundoff_is_separate_from_exact_scientific_records():
    old, fresh = record(), record(0.020000000000000004)
    assert scientific_record(old) == scientific_record(fresh)
    assert 0 < condition_drift(old, fresh) < 32 * 2.0**-52
    assert scientific_record(old) != scientific_record(record(matrix_hash="changed"))
    assert old["samples"][0]["reciprocal_condition_1"] == 0.02


def test_large_condition_change_exceeds_frozen_rounding_budget():
    assert condition_drift(record(), record(0.0200001)) > 32 * 2.0**-52
