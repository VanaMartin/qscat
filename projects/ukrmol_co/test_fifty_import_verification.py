"""An aggregate coverage pass must not hide missing or inaccurate high roots."""

import copy

import pytest

from projects.ukrmol_co.fifty_import_verification import IRREPS, covered_sectors


def coverage():
    probes = []
    for spin in (0, 2):
        for irrep in IRREPS:
            for roots in (10, 13) if irrep == "A1" else (5, 8):
                for space in (40, 80, 160):
                    probes.append(
                        dict(
                            spin=spin,
                            irrep=irrep,
                            roots=roots,
                            space=space,
                            cycles=600,
                            checkpoint_sha256="fixed",
                            all_roots_pass=True,
                            ci_converged=[True] * roots,
                            states=[
                                dict(
                                    energy_hartree=-100 + i / 10,
                                    spin_square=spin,
                                    physical_residual_norm_hartree=1e-11,
                                    penalized_residual_norm_hartree=1e-11,
                                    spin_penalty_action_norm_hartree=1e-11,
                                )
                                for i in range(roots)
                            ],
                        )
                    )
    return dict(all_roots_pass=True, probes=probes)


def test_complete_nonuniform_coverage_returns_all74_roots():
    selected = covered_sectors(coverage(), "fixed")
    assert sum(len(v) for v in selected.values()) == 74


@pytest.mark.parametrize(
    "mutation", ["duplicate", "checkpoint", "residual", "spin", "count", "high_energy"]
)
def test_corrupt_coverage_cannot_be_hidden_by_passing_flags(mutation):
    scan = coverage()
    if mutation == "duplicate":
        scan["probes"][-1] = copy.deepcopy(scan["probes"][-2])
    elif mutation == "checkpoint":
        scan["probes"][0]["checkpoint_sha256"] = "other"
    elif mutation == "residual":
        scan["probes"][0]["states"][-1]["physical_residual_norm_hartree"] = 2e-9
    elif mutation == "spin":
        scan["probes"][0]["states"][-1]["spin_square"] = float("nan")
    elif mutation == "count":
        scan["probes"][0]["ci_converged"].pop()
    else:
        scan["probes"][3]["states"][-1]["energy_hartree"] += 1e-5
    with pytest.raises((ValueError, AssertionError)):
        covered_sectors(scan, "fixed")
