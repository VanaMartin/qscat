"""Previously validated fixed-orbital CI residual controls for frozen diagnostics."""


def state_residual(solver, effective, ci, energy, ecore, ncas, electrons):
    import numpy as np

    norm = float(np.linalg.norm(ci))
    hci = solver.contract_2e(effective, ci, ncas, electrons)
    physical_hci = solver.undo_fix_spin().contract_2e(effective, ci, ncas, electrons)
    physical_rayleigh = float(np.vdot(ci, physical_hci).real / norm**2) + ecore
    return {
        "energy_hartree": float(energy),
        "spin_square": float(solver.spin_square(ci, ncas, electrons)[0]),
        "penalized_residual_norm_hartree": float(
            np.linalg.norm(hci - (energy - ecore) * ci) / norm
        ),
        "physical_rayleigh_energy_hartree": physical_rayleigh,
        "physical_residual_norm_hartree": float(
            np.linalg.norm(physical_hci - (physical_rayleigh - ecore) * ci) / norm
        ),
        "spin_penalty_action_norm_hartree": float(np.linalg.norm(hci - physical_hci) / norm),
    }


def benchmark():
    """Two-site/two-electron Hubbard singlet has a closed-form ground energy."""
    import numpy as np
    from pyscf import fci, lib

    lib.num_threads(4)
    h1 = np.array([[0.0, -1.0], [-1.0, 0.0]])
    h2 = np.zeros((2, 2, 2, 2))
    h2[0, 0, 0, 0] = h2[1, 1, 1, 1] = 4.0
    solver = fci.addons.fix_spin_(fci.direct_spin1.FCI(), shift=1.0, ss=0.0)
    solver.conv_tol, solver.conv_tol_residual = 1e-14, 1e-12
    offset = 3.0
    energy, ci = solver.kernel(h1, h2, 2, (1, 1), ecore=offset)
    expected = offset + (4 - np.sqrt(4**2 + 16)) / 2
    np.testing.assert_allclose(energy, expected, atol=1e-12, rtol=0)
    effective = solver.absorb_h1e(h1, h2, 2, (1, 1), 0.5)
    row = state_residual(solver, effective, ci, energy, offset, 2, (1, 1))
    assert row["physical_residual_norm_hartree"] <= 1e-12
    assert row["penalized_residual_norm_hartree"] <= 1e-12
    wrong = np.array(ci, copy=True)
    wrong[0, 0] += 0.01
    perturbed = state_residual(solver, effective, wrong, energy, offset, 2, (1, 1))
    assert perturbed["physical_residual_norm_hartree"] >= 1e-3
    return {
        "model": "Hubbard dimer, U=4, t=1, two electrons; energy offset 3",
        "analytic_ground_energy": float(expected),
        "exact_eigenvector": row,
        "perturbed_vector": perturbed,
    }


def solve(checkpoint, spin, irrep, roots, space, cycles):
    import numpy as np
    from pyscf import fci, lib, mcscf, scf

    lib.num_threads(4)
    mol = lib.chkfile.load_mol(str(checkpoint))
    mol.verbose = 4
    coefficients = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
    ncas = int(lib.chkfile.load(str(checkpoint), "mcscf/ncas"))
    ncore = int(lib.chkfile.load(str(checkpoint), "mcscf/ncore"))
    assert ncore == 2 and ncas in (11, 12)
    electrons = ((10 + spin) // 2, (10 - spin) // 2)
    cas = mcscf.CASCI(scf.RHF(mol), ncas, electrons, ncore=ncore)
    cas.canonicalization = False
    solver = fci.direct_spin1_symm.FCI(mol)
    solver.spin, solver.wfnsym = spin, irrep
    solver.nroots, solver.max_space, solver.max_cycle = roots, space, cycles
    solver.conv_tol, solver.conv_tol_residual, solver.lindep = 1e-12, 1e-9, 1e-22
    expected_spin = spin * (spin + 2) / 4
    cas.fcisolver = fci.addons.fix_spin_(solver, shift=1.0, ss=expected_spin)
    cas.kernel(coefficients)
    energies = np.atleast_1d(cas.e_tot)
    h1, ecore = cas.get_h1eff(coefficients)
    h2 = cas.get_h2eff(coefficients)
    effective = cas.fcisolver.absorb_h1e(h1, h2, ncas, electrons, 0.5)
    rows, vectors = [], []
    for energy, ci in zip(energies, cas.ci, strict=True):
        norm = float(np.linalg.norm(ci))
        rows.append(state_residual(cas.fcisolver, effective, ci, energy, ecore, ncas, electrons))
        vectors.append(ci.ravel() / norm)
    overlap = np.array(vectors) @ np.array(vectors).T
    orthogonality = float(np.max(np.abs(overlap - np.eye(roots))))
    flags = np.atleast_1d(cas.fcisolver.converged).tolist()
    return {
        "spin": spin,
        "irrep": irrep,
        "roots": roots,
        "space": space,
        "cycles": cycles,
        "ci_converged": flags,
        "states": rows,
        "ci_vector_orthogonality_error": orthogonality,
        "all_roots_pass": bool(
            all(flags)
            and orthogonality <= 1e-9
            and all(
                np.isfinite(list(row.values())).all()
                and abs(row["spin_square"] - expected_spin) <= 1e-6
                and row["penalized_residual_norm_hartree"] <= 1e-9
                and row["physical_residual_norm_hartree"] <= 1e-9
                for row in rows
            )
        ),
    }
