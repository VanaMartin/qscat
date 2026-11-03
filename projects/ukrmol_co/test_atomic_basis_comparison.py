"""Known offset cancellations and diagnostic rejection are not absolute convergence."""

import pytest
from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.atomic_basis_comparison import compare


def fixture(shift=0.0):
    qz = dict(
        finite_method_endpoint_reconstruction_verdict=True,
        fragment_sum_hartree={"C+O": -112.0, "C+O-": -112.05, "C-+O": -112.04},
        electron_affinity_ev={"O": 0.05 * HARTREE_TO_EV, "C": 0.04 * HARTREE_TO_EV},
    )
    f = qz | dict(
        fragment_sum_hartree={k: v - 0.1 + shift for k, v in qz["fragment_sum_hartree"].items()}
    )
    neutral = dict(
        runs={
            b: dict(
                result=dict(
                    basis=b,
                    bond_length_bohr=2.1323,
                    frozen_spatial_orbitals=2,
                    rhf_internally_stable=True,
                    rhf_externally_stable=True,
                    ccsd_converged=True,
                    lambda_converged=True,
                    neutral_energy_hartree=e,
                )
            )
            for b, e in [("aug-cc-pVQZ", -113.0), ("aug-cc-pV5Z", -113.1)]
        }
    )
    return (
        qz,
        f,
        neutral,
        dict(
            ea_basis_change_ev=0.01, atomic_to_fixed_molecular_reference_offset_basis_change_ev=0.02
        ),
    )


def test_matching_atomic_and_molecular_basis_shifts_cancel_without_a_de_claim():
    r = compare(*fixture())
    assert r["finite_basis_diagnostic_passed"] is True
    assert r["atomic_to_fixed_reference_offset_basis_change_ev"]["C+O"] == pytest.approx(
        0.0, abs=1e-12, rel=0
    )
    assert (
        r["reference_offset_is_dissociation_depth"]
        is r["complete_basis_convergence_qualified"]
        is False
    )


def test_large_reference_shift_rejects_without_basis_ladder():
    r = compare(*fixture(shift=0.002))
    assert (
        r["reference_offset_basis_diagnostic_passed"]
        is r["finite_basis_diagnostic_passed"]
        is False
    )
    assert r["no_automatic_larger_basis"] is True


def test_unqualified_or_duplicate_reference_is_not_used():
    q, f, n, b = fixture()
    n["runs"]["duplicate"] = next(iter(n["runs"].values()))
    with pytest.raises(ValueError, match="unique"):
        compare(q, f, n, b)
    q, f, n, b = fixture()
    n["runs"]["aug-cc-pV5Z"]["result"]["rhf_externally_stable"] = False
    with pytest.raises(ValueError, match="Unqualified"):
        compare(q, f, n, b)
