"""Inspect saved CO channel sparsity and early Numerov fill, without truncation."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import csc_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import digest, load_model, verify_inputs
from projects.ukrmol_co.numerov_typed_coupled import potential_derivatives, taylor_change


def support(matrix):
    """Exact stored nonzeros; small coefficients are retained."""
    sparse = csc_matrix(matrix)
    return dict(
        shape=list(matrix.shape),
        nonzeros=sparse.nnz,
        density=sparse.nnz / matrix.size,
        dense_bytes=matrix.nbytes,
        csc_bytes=sum(x.nbytes for x in [sparse.data, sparse.indices, sparse.indptr]),
        binary64_sha256=hashlib.sha256(np.asarray(matrix, dtype="<f8").tobytes()).hexdigest(),
    )


def factor_audit(a, b):
    """SuperLU fill and same-input binary64 differential solve diagnostics."""
    sparse = csc_matrix(a)
    factor = splu(sparse, permc_spec="COLAMD")
    dense_x = lu_solve(lu_factor(a), b)
    sparse_x = factor.solve(b)
    norm = lambda x: float(np.linalg.norm(x, ord=np.inf))  # noqa: E731
    denominator = norm(a) * norm(sparse_x) + norm(b)
    return dict(
        a=support(a),
        b=support(b),
        solution=support(dense_x),
        lu_nonzeros_including_both_diagonals=factor.L.nnz + factor.U.nnz,
        lu_fill_ratio=(factor.L.nnz + factor.U.nnz) / sparse.nnz,
        lu_csc_bytes=sum(
            x.nbytes
            for matrix in [factor.L, factor.U]
            for x in [matrix.data, matrix.indices, matrix.indptr]
        ),
        sparse_backward_error=norm(b - a @ sparse_x) / denominator if denominator else 0.0,
        relative_solution_difference=norm(dense_x - sparse_x) / max(norm(dense_x), 1.0),
    )


def inspect_model(model, index):
    """Four binary64 intervals from the retained typed-coupled recurrence algebra.

    This is a structure probe, not an exterior K calculation or an extended-
    precision solver qualification. No magnitude threshold removes entries.
    """
    energy = model["energies"][index]
    q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
    c2, c3 = model["c2"], model["c3"]
    identity = np.eye(len(q))
    diagonal = np.diag(model["thresholds"] - energy)
    radius, h = 18.0, 0.0015625

    def w(r):
        return diagonal + c2 / r**2 + c3 / r**3

    def b(r):
        return identity - h * h * w(r) / 12

    graph = (c2 != 0) | (c3 != 0)
    _, labels = connected_components(csc_matrix(graph), directed=False)
    across = labels[:, None] != labels[None, :]
    change = taylor_change(q, identity, potential_derivatives(diagonal, c2, c3, radius), h)
    b1 = b(radius + h)
    difference = -h * h * (w(radius + h) - w(radius)) / 12
    delta = lu_solve(lu_factor((b1 @ (q + change)).T), (difference @ q + b1 @ change).T).T
    initial = support(delta)
    rows = []
    for i in range(1, 5):
        current_w = w(radius + i * h)
        current_b = identity - h * h * current_w / 12
        increment = h * h * lu_solve(lu_factor(current_b), current_w)
        x = delta + increment
        a = identity + x
        row = dict(interval=i, potential_a=support(current_b), ratio_a=support(a))
        if i == 1:
            row["potential_factor"] = factor_audit(current_b, current_w)
            row["ratio_factor"] = factor_audit(a, x)
        delta = lu_solve(lu_factor(a), x)
        row["delta"] = support(delta)
        rows.append(row)
    return dict(
        energy_index=index,
        energy_rydberg_hex=float(energy).hex(),
        q=support(q),
        potential_component_sizes=sorted(np.bincount(labels).tolist(), reverse=True),
        q_cross_component_maximum_relative=float(np.max(abs(q[across])) / np.max(abs(q))),
        initial_delta=initial,
        intervals=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path, help="Checksum-pinned six-case corner raw directory")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    contract_path = Path("projects/ukrmol_co/numerov-corner-contract.json")
    contract = json.loads(contract_path.read_text())
    verify_inputs(args.raw, contract)
    rows = []
    for case in sorted({item["case"] for item in contract["identities"]}):
        model = load_model(args.raw / case)
        rows.append(
            dict(
                case=case,
                c2=support(model["c2"]),
                c3=support(model["c3"]),
                energies=[inspect_model(model, i) for i in [0, 48, 98]],
            )
        )
    verify_inputs(args.raw, contract)
    report = dict(
        scope="Exact-zero structure and binary64 early-interval LU only; no CO K qualification",
        platform=platform.platform(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        contract_sha256=digest(contract_path),
        input_sha256=contract["input_sha256"],
        source_sha256={
            name: digest(Path(name))
            for name in [
                "benchmarks/co_exterior_sparsity.py",
                "projects/ukrmol_co/native_full_reconstruction.py",
                "projects/ukrmol_co/native_boundary.py",
                "projects/ukrmol_co/numerov_typed_coupled.py",
            ]
        },
        step_bohr=0.0015625,
        radius_bohr=18.0,
        rows=rows,
    )
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
