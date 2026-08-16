"""V6: the library survives a production-scale assembly.

Grid dimensions are the real eMoScat N2 2-D deck
(`reference/eMoScat/input/experimental/N2-model.json`), but the potential here
is a generic analytic function: `libs/qscat` must not depend on `validation/`
or `projects/`, and the N2-specific assembly belongs to sub-project #6.

Deliberately does NOT factorize. A measured spike put that at 128 s and
13.6 GB peak RSS -- real, and no business in a routine test run.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, TensorGrid, hamiltonian_nd

MU = 12766.36


def _ecs_tail(base: float, n: int, *, skip: int = 2, alpha: float = 0.2) -> list[float]:
    """eMoScat `uniform_increment`/`exp`: `skip` elements at `base`, then growing."""
    return [base if i < skip else base * float(np.exp(alpha * (i - skip + 1))) for i in range(n)]


def _build(
    segments: list[tuple[float, float]],
    order: int,
    tail_base: float,
    n_tail: int,
    angle: float = 35.0,
) -> FemDvrEcsGrid:
    els: list[ElementSpec] = []
    start = 0.0
    for end, length in segments:
        k = round((end - start) / length)
        els += [ElementSpec((end - start) / k) for _ in range(k)]
        start = end
    els += [ElementSpec(h, angle) for h in _ecs_tail(tail_base, n_tail)]
    return FemDvrEcsGrid(GridSpec(quadrature=order, elements=els, x_min=0.0))


@pytest.mark.slow
def test_production_scale_2d_assembly() -> None:
    g_el = _build([(1.0, 0.2), (5.0, 1.0), (7.0, 2.0), (10.0, 3.0), (98.0, 4.0)], 8, 4.0, 15)
    g_nu = _build([(1.5, 0.5), (3.0, 0.15), (4.0, 0.5), (12.0, 1.0)], 14, 1.0, 10)

    assert (g_el.n, g_nu.n) == (335, 428)

    tg = TensorGrid([g_el, g_nu])
    assert tg.size == 143_380

    def V(r: npt.NDArray[np.complex128], R: npt.NDArray[np.complex128]) -> npt.ArrayLike:
        return 1.0 / (1.0 + r**2) + 1.0 / (1.0 + R**2)

    H = hamiltonian_nd(tg, [1.0, MU], V)
    assert H.shape == (143_380, 143_380)
    # matches eMoScat's own nnz formula, independently re-derived
    assert H.nnz == 3_276_450
    assert abs(H - H.T).max() < 1e-10  # complex symmetric, never Hermitian
