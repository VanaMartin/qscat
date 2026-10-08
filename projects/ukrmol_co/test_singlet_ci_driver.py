"""Physical controls for the experimental driver and mixed-spin CASSCF path."""

import numpy as np
import pytest

from projects.ukrmol_co.target import _fresh_ci_kernel, _sector_ci_driver

pyscf = pytest.importorskip("pyscf")
from pyscf import fci, gto, lib, mcscf, scf  # noqa: E402


@pytest.mark.parametrize("choice", ["spin1", "spin0"])
def test_hubbard_singlet_matches_analytic_energy_and_full_hamiltonian(choice):
    lib.num_threads(1)
    mol = gto.M(atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="sto-3g", symmetry="C2v", verbose=0)
    driver = _sector_ci_driver(fci, 0, "A1", {"target_singlet_a1_driver": choice})
    solver = driver(mol)
    solver.orbsym, solver.wfnsym = np.array([0, 0]), "A1"
    solver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
    solver.conv_tol, solver.conv_tol_residual = 1e-14, 1e-12
    h1 = np.array([[0.0, -1.0], [-1.0, 0.0]])
    h2 = np.zeros((2, 2, 2, 2))
    h2[0, 0, 0, 0] = h2[1, 1, 1, 1] = 4.0
    energy, ci = solver.kernel(h1, h2, 2, (1, 1), ecore=3.0)
    np.testing.assert_allclose(energy, 3 + (4 - np.sqrt(32)) / 2, atol=1e-12, rtol=0)
    effective = solver.absorb_h1e(h1, h2, 2, (1, 1), 0.5)
    action = fci.direct_spin1.contract_2e(effective, ci, 2, (1, 1))
    assert np.linalg.norm(action - (energy - 3) * ci) / np.linalg.norm(ci) <= 1e-12
    np.testing.assert_allclose(solver.spin_square(ci, 2, (1, 1))[0], 0, atol=1e-10, rtol=0)


@pytest.mark.parametrize("choice", ["spin1", "spin0"])
def test_mixed_spin_optimizer_preserves_full_ci_energies_and_spins(choice):
    lib.num_threads(1)
    mol = gto.M(atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="6-31g", symmetry="C2v", verbose=0)
    hf = scf.RHF(mol).run(conv_tol=1e-12)
    energies, solvers = [], []
    config = {"target_singlet_a1_driver": choice}
    for spin in (0, 2):
        driver = _sector_ci_driver(fci, spin, "A1", config)
        solver = driver(mol)
        solver.spin, solver.wfnsym = spin, "A1"
        solver.conv_tol, solver.conv_tol_residual = 1e-12, 1e-10
        solver = fci.addons.fix_spin_(solver, shift=1.0, ss=spin * (spin + 2) / 4)
        solvers.append(solver)
    mc = mcscf.state_average_mix_(mcscf.CASSCF(hf, 2, 2), solvers, [0.5, 0.5])
    for solver in solvers:
        solver.kernel = _fresh_ci_kernel(solver.kernel)
    mc.conv_tol, mc.conv_tol_grad = 1e-11, 1e-8
    mc.kernel()
    assert mc.converged and all(np.all(s.converged) for s in solvers)
    overlap = hf.mo_coeff[:, :2].T @ hf.get_ovlp() @ mc.mo_coeff[:, :2]
    assert min(np.linalg.svd(overlap, compute_uv=False)) < 1 - 1e-8
    reference = mcscf.CASCI(hf, 2, 2)
    h1, ecore = reference.get_h1eff(mc.mo_coeff)
    h2 = reference.get_h2eff(mc.mo_coeff)
    for spin in (0, 2):
        electrons = ((2 + spin) // 2, (2 - spin) // 2)
        oracle = fci.direct_spin1.FCI(mol)
        energy, ci = oracle.kernel(h1, h2, 2, electrons, ecore=ecore)
        np.testing.assert_allclose(
            oracle.spin_square(ci, 2, electrons)[0], spin * (spin + 2) / 4, atol=1e-9, rtol=0
        )
        energies.append(energy)
    np.testing.assert_allclose(mc.e_states, energies, atol=1e-9, rtol=0)
    np.testing.assert_allclose(mc.e_tot, np.mean(energies), atol=1e-9, rtol=0)
    for spin, solver, ci in zip((0, 2), solvers, mc.ci, strict=True):
        electrons = ((2 + spin) // 2, (2 - spin) // 2)
        np.testing.assert_allclose(
            solver.spin_square(ci, 2, electrons)[0], spin * (spin + 2) / 4, atol=1e-9, rtol=0
        )
