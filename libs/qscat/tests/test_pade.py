"""Order-N Pade propagator: accuracy vs `expm`, CN equivalence, convergence order.

`make_pade_stepper(H, dt, order)` approximates `exp(-i H dt)` to
`O(dt^(2*order+1))` per step. These tests gate: (1) order 1 == the CN stepper;
(2) the per-step error against dense `scipy.linalg.expm(-1j*H*dt)` scales as the
claimed order (order 3 far tighter than order 1); (3) the roots really encode a
valid `exp` approximant (scalar check). See `qscat/evolution/pade.py`.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from qscat.evolution import make_pade_stepper, make_sparse_cn_stepper, pade_roots
from scipy.linalg import expm


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A small complex-symmetric sparse H (the ECS case), well-conditioned."""
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    a = 0.5 * (a + a.T)  # symmetric (NOT Hermitian) -- like an ECS Hamiltonian
    a += np.diag(rng.standard_normal(n) + 5.0)  # dominant real diagonal
    return sp.csc_matrix(a)


def test_pade_roots_encode_exp_approximant() -> None:
    """Scalar check: prod_i (1 - z/r_i)/(1 + z/r_i) with z = i*h*dt reproduces
    exp(-i*h*dt) to the claimed order for each Pade order."""
    h, dt = 0.7, 0.3
    exact = np.exp(-1j * h * dt)
    prev_err = 1.0
    for order in (1, 2, 3):
        roots = pade_roots(order)
        z = 1j * h * dt
        approx = np.prod((1.0 - z / roots) / (1.0 + z / roots))
        err = abs(approx - exact)
        assert err < prev_err  # higher order is strictly more accurate
        prev_err = err
    assert prev_err < 1e-9  # order 3 essentially exact at this dt


def test_pade_order1_equals_crank_nicolson() -> None:
    """order=1 (root=2) must reproduce make_sparse_cn_stepper bit-for-bit."""
    n = 40
    H = _complex_symmetric(n, seed=1)
    dt = 0.5
    rng = np.random.default_rng(2)
    psi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_pade = make_pade_stepper(H, dt, order=1)(psi)
    x_cn = make_sparse_cn_stepper(H, dt)(psi)
    assert np.linalg.norm(x_pade - x_cn) <= 1e-12 * np.linalg.norm(x_cn)


def test_pade_per_step_accuracy_scales_with_order() -> None:
    """Per-step error vs dense expm(-i H dt): order 3 far tighter than order 1,
    each order improving on the previous. dt chosen so |H*dt|~1 (spectral
    radius ~10), the asymptotic regime where the order hierarchy is clean."""
    n = 30
    H = _complex_symmetric(n, seed=3)
    H_dense = H.toarray()
    dt = 0.1
    exact = expm(-1j * H_dense * dt)
    rng = np.random.default_rng(4)
    psi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    ref = exact @ psi
    errs = {}
    for order in (1, 2, 3):
        x = make_pade_stepper(H, dt, order=order)(psi)
        errs[order] = float(np.linalg.norm(x - ref) / np.linalg.norm(ref))
    assert errs[1] > errs[2] > errs[3]
    assert errs[1] > 1e-2  # CN visibly inexact
    assert errs[3] < 1e-4  # order 3 orders-of-magnitude tighter (~5e-6 measured)


def test_pade_order3_convergence_rate() -> None:
    """Halving dt cuts the order-3 per-step error by ~2^7 (O(dt^7)); the
    measured log-slope must exceed 6 (well above CN's 3)."""
    n = 24
    H = _complex_symmetric(n, seed=5)
    H_dense = H.toarray()
    rng = np.random.default_rng(6)
    psi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    errs = []
    dts = [0.2, 0.1]
    for dt in dts:
        ref = expm(-1j * H_dense * dt) @ psi
        x = make_pade_stepper(H, dt, order=3)(psi)
        errs.append(float(np.linalg.norm(x - ref) / np.linalg.norm(ref)))
    slope = np.log(errs[0] / errs[1]) / np.log(dts[0] / dts[1])
    assert slope > 6.0
