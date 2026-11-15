"""Singlet-driver selection and nontrivial-irrep full-Hamiltonian controls."""

import importlib.util
import unittest
from types import SimpleNamespace

import numpy as np

from projects.ukrmol_co.target import _fresh_ci_kernel, _sector_ci_driver


class DriverSelectionTests(unittest.TestCase):
    def test_only_requested_singlet_irreps_change(self):
        ordinary, restricted = object(), object()
        fci = SimpleNamespace(
            direct_spin1_symm=SimpleNamespace(FCI=ordinary),
            direct_spin0_symm=SimpleNamespace(FCI=restricted),
        )
        config = {"target_singlet_irrep_drivers": {"B1": "spin0", "B2": "spin0"}}
        for spin in (0, 2):
            for irrep in ("A1", "B1", "B2", "A2"):
                expected = restricted if spin == 0 and irrep in ("B1", "B2") else ordinary
                self.assertIs(_sector_ci_driver(fci, spin, irrep, config), expected)
                self.assertIs(_sector_ci_driver(fci, spin, irrep, {}), ordinary)
        legacy = {"target_singlet_a1_driver": "spin0"}
        self.assertIs(_sector_ci_driver(fci, 0, "A1", legacy), restricted)
        self.assertIs(_sector_ci_driver(fci, 0, "B1", legacy), ordinary)

    def test_invalid_driver_maps_reject(self):
        for overrides in (None, [], {"B3": "spin0"}, {"B1": "unknown"}):
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                _sector_ci_driver(None, 0, "B1", {"target_singlet_irrep_drivers": overrides})


@unittest.skipIf(importlib.util.find_spec("pyscf") is None, "PySCF is optional")
class NontrivialIrrepControls(unittest.TestCase):
    def test_five_singlets_match_analytic_spectrum_and_full_action(self):
        from pyscf import fci, gto, lib, symm

        lib.num_threads(1)
        mol = gto.M(
            atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="sto-3g", symmetry="C2v", verbose=0
        )
        diagonal = np.array([-1.0, -0.63, -0.21, 0.12, 0.43, 0.89])
        h1 = np.diag(diagonal)
        h2 = np.zeros((6, 6, 6, 6))
        for i in range(6):
            h2[i, i, i, i] = 4.0
        for i in range(3):
            for j in range(3, 6):
                h2[i, i, j, j] = h2[j, j, i, i] = 0.3
                for indices in ((i, j, j, i), (j, i, i, j), (i, j, i, j), (j, i, j, i)):
                    h2[indices] = 0.02
        expected = sorted(diagonal[i] + diagonal[j] + 0.32 for i in range(3) for j in range(3, 6))[
            :5
        ]
        effective = fci.direct_spin1.absorb_h1e(h1, h2, 6, (1, 1), 0.5)
        for irrep in ("B1", "B2"):
            with self.subTest(irrep=irrep):
                driver = _sector_ci_driver(
                    fci, 0, irrep, {"target_singlet_irrep_drivers": {irrep: "spin0"}}
                )
                solver = driver(mol)
                solver.orbsym = np.array([0] * 3 + [symm.irrep_name2id("C2v", irrep)] * 3)
                solver.wfnsym, solver.spin, solver.nroots = irrep, 0, 5
                solver.conv_tol, solver.conv_tol_residual = 1e-13, 1e-11
                solver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
                energies, vectors = solver.kernel(h1, h2, 6, (1, 1))
                self.assertTrue(np.all(solver.converged))
                np.testing.assert_allclose(energies, expected, atol=1e-11, rtol=0)
                for energy, vector in zip(energies, vectors, strict=True):
                    action = fci.direct_spin1.contract_2e(effective, vector, 6, (1, 1))
                    self.assertLessEqual(
                        np.linalg.norm(action - energy * vector) / np.linalg.norm(vector), 1e-10
                    )
                    np.testing.assert_allclose(
                        fci.spin_op.spin_square(vector, 6, (1, 1))[0], 0, atol=1e-10, rtol=0
                    )

    def test_moving_mixed_spin_orbitals_match_dense_physical_oracle(self):
        from pyscf import fci, gto, lib, mcscf, scf, symm

        lib.num_threads(1)
        mol = gto.M(
            atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="cc-pVDZ", symmetry="C2v", verbose=0
        )
        hf = scf.RHF(mol).run(conv_tol=1e-12)
        for irrep in ("B1", "B2"):
            with self.subTest(irrep=irrep):
                mc = mcscf.CASSCF(hf, 2, 2)
                initial = mcscf.sort_mo_by_irrep(mc, hf.mo_coeff, {"A1": 1, irrep: 1})
                labels = [(0, "A1"), (0, irrep), (2, irrep)]
                config = {"target_singlet_irrep_drivers": {irrep: "spin0"}}
                solvers = []
                for spin, sector in labels:
                    solver = _sector_ci_driver(fci, spin, sector, config)(mol)
                    solver.spin, solver.wfnsym = spin, sector
                    solver.conv_tol, solver.conv_tol_residual = 1e-12, 1e-10
                    solvers.append(
                        fci.addons.fix_spin_(solver, shift=1.0, ss=spin * (spin + 2) / 4)
                    )
                mc = mcscf.state_average_mix_(mc, solvers, np.full(3, 1 / 3))
                for solver in solvers:
                    solver.kernel = _fresh_ci_kernel(solver.kernel)
                mc.conv_tol, mc.conv_tol_grad = 1e-11, 1e-8
                # Match the target path's AH precision, rather than its looser defaults.
                mc.max_cycle_macro = 150
                mc.ah_lindep, mc.ah_conv_tol, mc.ah_start_tol = 1e-24, 1e-16, 1e-8
                mc.kernel(initial)
                self.assertTrue(
                    mc.converged and all(np.all(s.converged) for s in solvers),
                    f"orbital={mc.converged}, CI={[s.converged for s in solvers]}, "
                    f"energies={mc.e_states}",
                )
                overlap = initial[:, :2].T @ hf.get_ovlp() @ mc.mo_coeff[:, :2]
                self.assertLess(min(np.linalg.svd(overlap, compute_uv=False)), 1 - 1e-8)
                h1, ecore = mc.get_h1eff()
                h2 = mc.get_h2eff()
                effective = fci.direct_spin1.absorb_h1e(h1, h2, 2, (1, 1), 0.5)
                basis = np.eye(4).reshape(4, 2, 2)
                physical = np.column_stack(
                    [fci.direct_spin1.contract_2e(effective, c, 2, (1, 1)).ravel() for c in basis]
                )
                np.testing.assert_allclose(physical, physical.T, atol=1e-12, rtol=0)
                eigenvalues, eigenvectors = np.linalg.eigh(physical)
                ids = [symm.irrep_name2id("C2v", name) for name in ("A1", irrep)]
                dense_states = []
                for energy, vector in zip(eigenvalues, eigenvectors.T, strict=True):
                    ci = vector.reshape(2, 2)
                    spin_square = fci.spin_op.spin_square(ci, 2, (1, 1))[0]
                    for sector in ("A1", irrep):
                        mask = np.array(
                            [
                                [
                                    ids[a] ^ ids[b] == symm.irrep_name2id("C2v", sector)
                                    for b in range(2)
                                ]
                                for a in range(2)
                            ]
                        )
                        if np.sum(abs(ci[mask]) ** 2) >= 1 - 1e-9:
                            dense_states.append((energy + ecore, spin_square, sector))
                expected = [
                    min(
                        e
                        for e, ss, sector in dense_states
                        if sector == label and abs(ss - spin * (spin + 2) / 4) <= 1e-8
                    )
                    for spin, label in labels
                ]
                np.testing.assert_allclose(mc.e_states, expected, atol=1e-9, rtol=0)
                np.testing.assert_allclose(mc.e_tot, np.mean(expected), atol=1e-9, rtol=0)


if __name__ == "__main__":
    unittest.main()
