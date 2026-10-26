"""Small independent two-electron oracle and bounded A2 driver controls."""

from __future__ import annotations

import math
from itertools import product

import numpy as np


def two_electron_matrix(h1, eri):
    """Second-quantized M_s=0 matrix, with alpha modes before beta modes."""
    n = len(h1)
    if h1.shape != (n, n) or eri.shape != (n, n, n, n):
        raise ValueError("Square one-/two-electron tensors required")
    basis = [(1 << a) | (1 << (n + b)) for a, b in product(range(n), repeat=2)]
    addresses = {state: i for i, state in enumerate(basis)}

    def apply(state, operators):
        sign = 1
        for mode, create in operators:
            mode = int(mode)
            occupied = (state >> mode) & 1
            if occupied == create:
                return None, 0
            sign *= -1 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1
            state ^= 1 << mode
        return state, sign

    matrix = np.zeros((n * n, n * n))
    for column, state in enumerate(basis):
        for p, q, spin in product(range(n), range(n), range(2)):
            final, sign = apply(state, [(q + spin * n, False), (p + spin * n, True)])
            if sign:
                matrix[addresses[final], column] += sign * h1[p, q]
        for p, q, r, s in zip(*np.nonzero(eri), strict=True):
            for sigma, tau in product(range(2), repeat=2):
                final, sign = apply(
                    state,
                    [
                        (q + sigma * n, False),
                        (s + tau * n, False),
                        (r + tau * n, True),
                        (p + sigma * n, True),
                    ],
                )
                if sign:
                    matrix[addresses[final], column] += 0.5 * sign * eri[p, q, r, s]
    return matrix


def spin_blocks(matrix, orbsym, wfnsym=1):
    """Project symmetric singlets/antisymmetric triplets in one spatial block."""
    n = len(orbsym)
    vectors = [[], []]
    for a in range(n):
        for b in range(a, n):
            if orbsym[a] ^ orbsym[b] != wfnsym:
                continue
            for spin in (0, 1):
                if spin and a == b:
                    continue
                v = np.zeros((n, n))
                v[a, b] = 1 if a == b else 1 / math.sqrt(2)
                if a != b:
                    v[b, a] = (-1 if spin else 1) / math.sqrt(2)
                vectors[spin].append(v.ravel())
    return [np.array(rows) @ matrix @ np.array(rows).T for rows in vectors]


def synthetic_integrals():
    """C2v-covariant positive pair metric with known non-A1 singlet/triplet blocks."""
    h1 = np.diag([-1.2, -0.6, 0.1, 0.7])
    eri = np.zeros((4, 4, 4, 4))

    def put(p, q, r, s, value):
        for a, b, c, d in [
            (p, q, r, s),
            (q, p, r, s),
            (p, q, s, r),
            (q, p, s, r),
            (r, s, p, q),
            (s, r, p, q),
            (r, s, q, p),
            (s, r, q, p),
        ]:
            eri[a, b, c, d] = value

    for i, value in enumerate([2.7, 2.2, 1.9, 1.4]):
        put(i, i, i, i, value)
    for p, q, j, k in [(0, 3, 0.8, 0.12), (1, 2, 0.65, 0.08)]:
        put(p, p, q, q, j)
        put(p, q, p, q, k)
    for p, q in [(0, 1), (2, 3)]:
        put(p, q, p, q, 0.2)
    for p, q in [(0, 2), (1, 3)]:
        put(p, q, p, q, 0.15)
    put(0, 1, 3, 2, 0.09)
    put(0, 2, 3, 1, 0.04)
    return h1, eri, np.array([0, 2, 3, 1])


def exact_controls():
    """Validate both external drivers against analytic and independent Fock matrices."""
    from pyscf import fci, lib

    lib.num_threads(1)
    h1, eri, labels = synthetic_integrals()
    matrix = two_electron_matrix(h1, eri)
    singlet, triplet = spin_blocks(matrix, labels)
    np.testing.assert_allclose(singlet, [[0.42, 0.13], [0.13, 0.23]], atol=1e-13, rtol=0)
    np.testing.assert_allclose(triplet, [[0.18, 0.05], [0.05, 0.07]], atol=1e-13, rtol=0)
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-13, rtol=0)
    swap = np.eye(16).reshape(16, 4, 4).transpose(0, 2, 1).reshape(16, 16)
    spin = np.eye(16) - swap
    np.testing.assert_allclose(matrix @ spin, spin @ matrix, atol=1e-13, rtol=0)
    sectors = (labels[:, None] ^ labels).ravel()
    np.testing.assert_allclose(matrix[sectors[:, None] != sectors], 0, atol=0, rtol=0)
    expected = np.array(
        [0.325 - math.sqrt(0.095**2 + 0.13**2), 0.325 + math.sqrt(0.095**2 + 0.13**2)]
    )
    np.testing.assert_allclose(np.linalg.eigvalsh(singlet), expected, atol=1e-13, rtol=0)
    triplet_energy, triplet_vectors = np.linalg.eigh(triplet)
    negative = np.zeros((4, 4))
    for value, (i, j) in zip(triplet_vectors[:, 0], [(0, 3), (1, 2)], strict=True):
        negative[i, j] = value / math.sqrt(2)
        negative[j, i] = -value / math.sqrt(2)
    negative_residual = float(
        np.linalg.norm(matrix @ negative.ravel() - triplet_energy[0] * negative.ravel())
    )
    negative_spin = float(np.vdot(negative.ravel(), spin @ negative.ravel()).real)
    assert negative_residual <= 1e-11 and abs(negative_spin - 2) <= 1e-11
    wrong = negative.copy()
    wrong[0, 3] += 0.01
    perturbed_residual = float(
        np.linalg.norm(matrix @ wrong.ravel() - triplet_energy[0] * wrong.ravel())
    )
    assert perturbed_residual > 1e-3
    rows = []
    for name, driver in [
        ("spin1", fci.direct_spin1_symm.FCI),
        ("spin0", fci.direct_spin0_symm.FCI),
    ]:
        solver = driver()
        solver.orbsym, solver.wfnsym, solver.spin = labels, 1, 0
        solver.nroots, solver.max_space, solver.max_cycle = 2, 160, 200
        solver.conv_tol, solver.conv_tol_residual, solver.lindep = 1e-12, 1e-9, 1e-22
        solver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
        energies, vectors = solver.kernel(h1, eri, 4, (1, 1))
        np.testing.assert_allclose(energies, expected, atol=1e-11, rtol=0)
        assert np.all(solver.converged)
        for energy, ci in zip(energies, vectors, strict=True):
            assert np.linalg.norm(matrix @ ci.ravel() - energy * ci.ravel()) <= 1e-11
            assert np.linalg.norm(ci - ci.T) <= 1e-11
            np.testing.assert_allclose(solver.spin_square(ci, 4, (1, 1))[0], 0, atol=1e-10, rtol=0)
        rows.append(
            dict(
                driver=name,
                energy_hartree=np.asarray(energies).tolist(),
                ci_vectors=[np.asarray(ci).tolist() for ci in vectors],
                ci_converged=np.atleast_1d(solver.converged).tolist(),
                passed=True,
            )
        )
    return dict(
        drivers=rows,
        analytic_singlet_energy_hartree=expected.tolist(),
        singlet_block=singlet.tolist(),
        triplet_block=triplet.tolist(),
        h1=h1.tolist(),
        eri=eri.tolist(),
        orbsym=labels.tolist(),
        negative_controls=dict(
            exact_triplet_residual=negative_residual,
            exact_triplet_spin_square=negative_spin,
            perturbed_vector_residual=perturbed_residual,
        ),
    )


def mixed_optimizer_control():
    """Nontrivial A2 mixed-spin orbital movement, checked against the final exact block."""
    from pyscf import ao2mo, fci, gto, lib, mcscf, scf, symm

    from projects.ukrmol_co.target import _fresh_ci_kernel

    lib.num_threads(1)
    mol = gto.M(
        atom="H 1 0 0; H -1 0 0; H 0 1.6 0; H 0 -1.6 0",
        unit="Bohr",
        basis="sto-3g",
        symmetry="C2v",
        verbose=0,
    )
    hf = scf.RHF(mol).run(conv_tol=1e-12)
    assert hf.converged
    solvers = []
    for spin, driver in [(0, fci.direct_spin0_symm.FCI), (2, fci.direct_spin1_symm.FCI)]:
        solver = driver(mol)
        solver.spin, solver.wfnsym = spin, "A2"
        solver.max_cycle, solver.max_space, solver.lindep = 200, 160, 1e-22
        solver.conv_tol, solver.conv_tol_residual = 1e-12, 1e-10
        solvers.append(fci.addons.fix_spin_(solver, shift=1.0, ss=spin * (spin + 2) / 4))
    mc = mcscf.CASSCF(hf, 2, 2, ncore=1)
    initial = mcscf.sort_mo_by_irrep(mc, hf.mo_coeff, {"B1": 1, "B2": 1}, {"A1": 1})
    mc = mcscf.state_average_mix_(mc, solvers, [0.5, 0.5])
    for solver in solvers:
        solver.kernel = _fresh_ci_kernel(solver.kernel)
    mc.conv_tol, mc.conv_tol_grad, mc.max_cycle_macro = 1e-11, 1e-8, 50
    mc.kernel(initial)
    assert mc.converged and all(np.all(s.converged) for s in solvers)
    overlap = initial[:, :1].T @ hf.get_ovlp() @ mc.mo_coeff[:, :1]
    movement = float(np.linalg.svd(overlap, compute_uv=False).min())
    assert movement < 1 - 1e-8
    h1, ecore = mc.get_h1eff(mc.mo_coeff)
    eri = ao2mo.restore(1, mc.get_h2eff(mc.mo_coeff), 2)
    labels = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mc.mo_coeff[:, 1:3])
    ids = np.array([{"A1": 0, "A2": 1, "B1": 2, "B2": 3}[str(label)] for label in labels])
    matrix = two_electron_matrix(h1, eri)
    blocks = spin_blocks(matrix, ids)
    reference = np.array([float(block[0, 0]) + ecore for block in blocks])
    np.testing.assert_allclose(mc.e_states, reference, atol=1e-9, rtol=0)
    np.testing.assert_allclose(mc.e_tot, reference.mean(), atol=1e-9, rtol=0)
    actions = []
    for spin, solver, ci in zip((0, 2), solvers, mc.ci, strict=True):
        electrons = ((2 + spin) // 2, (2 - spin) // 2)
        np.testing.assert_allclose(
            solver.spin_square(ci, 2, electrons)[0], spin * (spin + 2) / 4, atol=1e-9, rtol=0
        )
        effective = fci.direct_spin1.absorb_h1e(h1, eri, 2, electrons, 0.5)
        action = fci.direct_spin1.contract_2e(effective, ci, 2, electrons)
        assert np.linalg.norm(action - (mc.e_states[spin // 2] - ecore) * ci) <= 1e-9
        actions.append(np.asarray(action).tolist())
    return dict(
        passed=True,
        core_overlap_singular_value=movement,
        energy_hartree=np.asarray(mc.e_states).tolist(),
        independent_energy_hartree=reference.tolist(),
        ensemble_energy_hartree=float(mc.e_tot),
        h1=h1.tolist(),
        eri=eri.tolist(),
        orbsym=ids.tolist(),
        ecore=float(ecore),
        initial_mo_coefficients=initial.tolist(),
        final_mo_coefficients=mc.mo_coeff.tolist(),
        ao_overlap=hf.get_ovlp().tolist(),
        ci_vectors=[np.asarray(ci).tolist() for ci in mc.ci],
        full_actions=actions,
        optimizer_converged=bool(mc.converged),
        ci_converged=[bool(np.all(s.converged)) for s in solvers],
    )
