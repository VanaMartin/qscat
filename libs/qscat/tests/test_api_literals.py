"""lib-M14 (2026-08-25 API surface pass): the free-form str parameters carry
Literal types, and each Literal is the single source of its legal values --
the runtime `raise` paths (which stay, for callers holding a plain str)
validate against the SAME tuple the type checker sees."""

from __future__ import annotations

from typing import get_args


def test_method_literal_names_the_three_extraction_methods() -> None:
    from qscat.core.time_dependent import Method

    assert get_args(Method) == ("tw", "delta", "flow")


def test_axis_literal_is_the_single_source_of_the_axes_tuple() -> None:
    from qscat.core.td_extractors import _AXES, Axis

    assert get_args(Axis) == ("electronic", "nuclear")
    assert _AXES == get_args(Axis)


def test_channel_literal_names_the_two_mesh_channels() -> None:
    from qscat.tuning.propose import Channel

    assert get_args(Channel) == ("ve", "dissociation")


def test_refinement_coordinate_literal_names_the_two_grids() -> None:
    # Public through the package, since it is the vocabulary a recorded
    # refinement step speaks (`Refine2dStep["coordinate"]`).
    from qscat.tuning import RefinementCoordinate
    from qscat.tuning.refine2d import RefinementCoordinate as module_coordinate

    assert RefinementCoordinate is module_coordinate
    assert get_args(RefinementCoordinate) == ("electronic", "nuclear")


def test_verdict_literal_is_public_and_names_the_seven_verdicts() -> None:
    # Public through both the module and the qscat.core re-export (it is the
    # vocabulary `OverlapPair.verdict` speaks; typed user code needs it).
    from qscat.core import Verdict as core_verdict
    from qscat.core.assignment import OverlapPair, Verdict

    assert core_verdict is Verdict
    assert get_args(Verdict) == (
        "ok",
        "spurious",
        "basis-limited",
        "box-limited",
        "weak",
        "mixed",
        "distant",
    )
    assert OverlapPair.__dataclass_fields__["verdict"].type == "Verdict"
