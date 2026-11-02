"""Missing versus changed prerequisite files are distinct; no model or ODE."""

import pytest

from projects.ukrmol_co.verify_stretched_real_completion_stop import sha, source_prerequisites


def test_missing_pointer_is_named_instead_of_qualifying_the_capsule():
    assert source_prerequisites(
        {"present": b"ok"}, {"present": sha(b"ok"), "missing": sha(b"pointer")}
    ) == ["missing"]


def test_changed_present_source_is_not_misreported_as_missing():
    with pytest.raises(ValueError, match="Changed present prerequisite"):
        source_prerequisites({"present": b"changed"}, {"present": sha(b"ok")})
