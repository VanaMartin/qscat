"""Exercise PySCF checkpoint projection across larger and smaller AO bases."""

import numpy as np
import pytest

from projects.ukrmol_co.target import project_initial_orbitals

gto = pytest.importorskip("pyscf.gto")
lib = pytest.importorskip("pyscf.lib")
mcscf = pytest.importorskip("pyscf.mcscf")
scf = pytest.importorskip("pyscf.scf")
symm = pytest.importorskip("pyscf.symm")


@pytest.fixture(scope="module")
def orbitals():
    lib.num_threads(1)
    result = {}
    for basis, bond in (("cc-pvdz", 2.1323), ("sto-3g", 2.1323), ("sto-3g", 2.5)):
        mol = gto.M(
            atom=f"C 0 0 0; O 0 0 {bond}", unit="Bohr", basis=basis, symmetry="C2v", verbose=0
        )
        hf = scf.RHF(mol).run(conv_tol=1e-11)
        assert hf.converged
        mc = mcscf.CASSCF(hf, 6, 10, ncore=2)
        mo = mcscf.sort_mo_by_irrep(mc, hf.mo_coeff, {"A1": 4, "B1": 1, "B2": 1}, {"A1": 2})
        result[basis, bond] = mc, mo
    return result


def test_downprojection_preserves_core_and_active_orbitals(orbitals):
    source, mo = orbitals["cc-pvdz", 2.1323]
    dest, _ = orbitals["sto-3g", 2.1323]
    with pytest.raises(RuntimeError, match="Too many orbitals"):
        mcscf.project_init_guess(dest, mo, prev_mol=source.mol)
    actual = project_initial_orbitals(dest, mo, source.mol)
    expected = mcscf.project_init_guess(dest, mo[:, :8], prev_mol=source.mol)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)
    metric = dest._scf.get_ovlp()
    assert actual.shape == dest._scf.mo_coeff.shape
    np.testing.assert_allclose(
        actual.T @ metric @ actual, np.eye(actual.shape[1]), atol=1e-10, rtol=0
    )
    labels = symm.label_orb_symm(dest.mol, dest.mol.irrep_name, dest.mol.symm_orb, actual)
    assert list(labels[:2]) == ["A1", "A1"]
    assert {name: np.count_nonzero(labels[2:8] == name) for name in ("A1", "B1", "B2", "A2")} == {
        "A1": 4,
        "B1": 1,
        "B2": 1,
        "A2": 0,
    }


@pytest.mark.parametrize("source_basis", ["cc-pvdz", "sto-3g"])
def test_other_projection_routes_match_the_original_api(orbitals, source_basis):
    source, mo = orbitals[source_basis, 2.1323]
    dest, _ = orbitals["cc-pvdz", 2.1323]
    expected = mcscf.project_init_guess(dest, mo, prev_mol=source.mol)
    actual = project_initial_orbitals(dest, mo, source.mol)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)


def test_simultaneous_geometry_and_basis_change_is_still_rejected(orbitals):
    source, mo = orbitals["cc-pvdz", 2.1323]
    dest, _ = orbitals["sto-3g", 2.5]
    with pytest.raises(NotImplementedError, match="different system"):
        project_initial_orbitals(dest, mo, source.mol)
