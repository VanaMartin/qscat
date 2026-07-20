import experiment


def test_green_groups_pass_pending_never_fail():
    checks = experiment.run_checks()
    statuses = {c[2] for c in checks}
    assert statuses <= {"PASS", "PENDING", "FAIL"}
    # Every Group A and Group C-integrity check must PASS; nothing may FAIL.
    assert not any(c[2] == "FAIL" for c in checks), [c for c in checks if c[2] == "FAIL"]
    assert any(c[2] == "PASS" for c in checks)
    assert any(c[2] == "PENDING" for c in checks)  # resonance / cross-section anchors


def test_main_exits_zero():
    assert experiment.main() == 0
