"""Keep the CI trial-space refinement separate from orbital-ensemble changes."""

import sys

import pytest

from projects.ukrmol_co import run
from projects.ukrmol_co.target import ci_trial_space


@pytest.mark.parametrize(("roots", "expected"), [(1, 40), (5, 40), (10, 80)])
def test_default_preserves_the_established_trial_space(roots, expected):
    assert ci_trial_space(roots, None) == expected


@pytest.mark.parametrize("space", [40, 80, 160])
def test_override_controls_trial_vectors_without_changing_roots(space):
    assert ci_trial_space(5, space) == space


@pytest.mark.parametrize("space", [-1, 0, 5])
def test_trial_space_must_exceed_requested_ensemble_roots(space):
    with pytest.raises(ValueError, match="exceed"):
        ci_trial_space(5, space)


@pytest.mark.parametrize("space", [0, 5, 10])
def test_cli_rejects_insufficient_space_for_nonuniform_ensemble(monkeypatch, tmp_path, space):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run",
            "--workdir",
            str(tmp_path / "run"),
            "--model",
            "CAS-A",
            "--orbitals",
            "state-averaged",
            "--qc-only",
            "--target-roots",
            "10",
            "--sa-singlet-roots",
            "10",
            "5",
            "5",
            "5",
            "--sa-triplet-roots",
            "10",
            "5",
            "5",
            "5",
            "--target-ci-max-space",
            str(space),
        ],
    )
    with pytest.raises(SystemExit) as error:
        run.main()
    assert error.value.code == 2
    assert not (tmp_path / "run").exists()
