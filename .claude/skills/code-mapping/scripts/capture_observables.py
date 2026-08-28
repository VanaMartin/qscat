#!/usr/bin/env python3
"""Capture numeric observables so a refactor can be proven to change nothing.

This is the lane-B gate of the structure audit. Run it BEFORE a change that touches
a numerics path, make the change, run it again, and compare bitwise. A
behaviour-preserving refactor reorders no floating-point operation, so identical
bits are the correct bar; a mismatch has found a change in behaviour, and the
answer is to explain or revert it, never to loosen the comparison.

Cases are deliberately toy-scale — this detects change, it does not validate physics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Callable
from pathlib import Path

import numpy as np


def _hash(values: np.ndarray) -> str:
    """Bitwise digest of an array, independent of how it is printed."""
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def case_dvr_eigenvalues() -> tuple[list[str], np.ndarray]:
    """Lowest oscillator eigenvalues on a small FEM-DVR-ECS grid."""
    from qscat.dvr import FemDvrEcsGrid, eigen, hamiltonian
    from qscat.dvr.spec import ElementSpec, GridSpec

    spec = GridSpec(
        quadrature=6,
        elements=[ElementSpec(2.0), ElementSpec(2.0), ElementSpec(2.0, 20.0)],
    )
    grid = FemDvrEcsGrid(spec)
    energies, _ = eigen(hamiltonian(grid, 0.5 * grid.points**2, mass=1.0))
    return ["qscat.dvr", "qscat.ecs"], np.asarray(energies[:6])


def case_sparse_lu_solve() -> tuple[list[str], np.ndarray]:
    """Solution of a small complex-symmetric system through SparseLU."""
    import scipy.sparse as sp
    from qscat.linalg import SparseLU

    n = 40
    diag = np.arange(1, n + 1, dtype=complex) + 0.5j
    matrix = sp.diags([diag, -np.ones(n - 1), -np.ones(n - 1)], [0, 1, -1], format="csc")
    rhs = np.ones(n, dtype=complex)
    return ["qscat.linalg"], SparseLU(matrix, backend="scipy").solve(rhs)


CASES: dict[str, Callable[[], tuple[list[str], np.ndarray]]] = {
    "dvr_eigenvalues": case_dvr_eigenvalues,
    "sparse_lu_solve": case_sparse_lu_solve,
}


def main() -> None:
    """Run the selected cases and write their digests and values."""
    parser = argparse.ArgumentParser(description="Capture observables for the lane-B gate.")
    parser.add_argument("--out", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--case",
        action="append",
        choices=sorted(CASES),
        help="Case to run; repeatable. Default: all.",
    )
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for name in sorted(args.case or CASES):
        modules, values = CASES[name]()
        array = np.asarray(values)
        (args.out / f"{name}.json").write_text(
            json.dumps(
                {"real": array.real.tolist(), "imag": array.imag.tolist()},
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        manifest.append({"case": name, "modules": sorted(modules), "sha256": _hash(array)})
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
