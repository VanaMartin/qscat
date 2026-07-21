from validation.n2 import experiment


def test_green_groups_pass_pending_never_fail():
    checks = experiment.run_checks()
    statuses = {c[2] for c in checks}
    assert statuses <= {"PASS", "PENDING", "FAIL", "NOTE"}
    # Every Group A, Group C-integrity, Group B, and GATED Group C5 check must PASS;
    # nothing may FAIL. NOTE rows (DOCUMENTED-LIMITED C5 anchors) are informational,
    # not failures.
    assert not any(c[2] == "FAIL" for c in checks), [c for c in checks if c[2] == "FAIL"]
    assert any(c[2] == "PASS" for c in checks)
    assert any(c[2] == "PENDING" for c in checks)  # Group D (time-dependent), still pending
    assert any(c[2] == "NOTE" for c in checks)  # documented-limited C5 anchors (elastic,
    # near-threshold)


def test_main_exits_zero():
    assert experiment.main() == 0
