"""Original launch/service/native exits must not be silently rewritten by qualification."""

import pytest

from projects.ukrmol_co.verify_stretched_augtz_preparation import original_exits


def test_all_original_exit_sources_are_required():
    owner = dict(original_command_exit=0, original_container_exit=0, original_owner_exit=0)
    result = dict(scientific_exit=0, stages=[dict(original_exit=0)])
    resources = dict(exit_code=0)
    service = dict(original_service_exit=0)
    process = dict(original_exit=0)
    assert original_exits(owner, result, resources, service, process) == [0] * 8
    for key in owner:
        with pytest.raises(ValueError):
            original_exits(owner | {key: 2}, result, resources, service, process)
    with pytest.raises(ValueError):
        original_exits(owner, result, resources, service, dict(original_exit=2))
    with pytest.raises(ValueError):
        original_exits(owner, result, resources, dict(original_service_exit=1), process)
