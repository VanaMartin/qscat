"""Experimental native-residual mixed solves around the retained coupled oracle."""

from __future__ import annotations

import warnings
from contextlib import contextmanager

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve

import projects.ukrmol_co.numerov_typed_coupled as reference
from projects.ukrmol_co.numerov_native_residual import native_residual

REFERENCE_SOLVE = reference.refined_solve


def native_refined_solve(a, b, stats, *, observer=None):
    """Mirror the frozen refinement algorithm, substituting only its residual."""
    if len(a) <= 3:
        return REFERENCE_SOLVE(a, b, stats)
    dtype = a.dtype.type
    if b.dtype != a.dtype or a.dtype.kind != "f":
        raise ValueError("Compatible real arithmetic required")

    def norm(x):
        return np.max(np.sum(abs(x), axis=1))

    with warnings.catch_warnings():
        warnings.simplefilter("error", LinAlgWarning)
        factors = lu_factor(np.asarray(a, float), check_finite=True)
        x = np.asarray(lu_solve(factors, np.asarray(b, float)), dtype=dtype)
        for iteration in range(4):
            residual = native_residual(a, b, x)
            denominator = norm(a) * norm(x) + norm(b)
            error = norm(residual) / denominator if denominator else dtype(0)
            if observer is not None:
                observer(a, b, x, residual, error, iteration)
            if not np.isfinite(error):
                raise ValueError("Nonfinite mixed-precision residual")
            if error <= 32 * np.finfo(dtype).eps:
                stats["maximum_solve_backward_error"] = max(
                    stats["maximum_solve_backward_error"], float(error)
                )
                stats["maximum_refinement_corrections"] = max(
                    stats["maximum_refinement_corrections"], iteration
                )
                stats["large_system_solves"] += 1
                return x
            if iteration < 3:
                x += np.asarray(lu_solve(factors, np.asarray(residual, float)), dtype=dtype)
    raise np.linalg.LinAlgError("Mixed-precision refinement exhausted three corrections")


@contextmanager
def installed_solver(solver):
    """Single-worker experiment hook; restore the oracle binding on every exit."""
    if reference.refined_solve is not REFERENCE_SOLVE:
        raise RuntimeError("Coupled solver is already instrumented")
    reference.refined_solve = solver
    try:
        yield
    finally:
        reference.refined_solve = REFERENCE_SOLVE


def native_coupled_k(*args, **kwargs):
    """Same finite exterior equation; native large solves and original small path."""
    kwargs.setdefault("dtype", np.longdouble)
    with installed_solver(native_refined_solve):
        return reference.typed_coupled_k(*args, **kwargs)


class ResidualAudit:
    """Independently reevaluate every observed residual and accept/reject decision."""

    def __init__(self):
        self.evaluations = 0
        self.accepted_solves = 0
        self.exact_residuals = 0
        self.gate_disagreements = 0
        self.maximum_scaled_difference = 0.0
        self.maximum_accepted_numpy_error = 0.0
        self.maximum_iteration = 0
        self.precision_preserved = True
        self.last = None

    def __call__(self, a, b, x, residual, error, iteration):
        expected = b - a @ x

        def norm(value):
            return np.max(np.sum(abs(value), axis=1))

        denominator = norm(a) * norm(x) + norm(b)
        numpy_error = norm(expected) / denominator if denominator else a.dtype.type(0)
        difference = (
            norm(expected - residual) / denominator if denominator else norm(expected - residual)
        )
        bound = 32 * np.finfo(a.dtype).eps
        accepted = error <= bound
        self.accepted_solves += int(accepted)
        self.evaluations += 1
        self.exact_residuals += int(np.array_equal(expected, residual))
        self.gate_disagreements += int(accepted != (numpy_error <= bound))
        self.maximum_scaled_difference = max(self.maximum_scaled_difference, float(difference))
        if accepted:
            self.maximum_accepted_numpy_error = max(
                self.maximum_accepted_numpy_error, float(numpy_error)
            )
        self.maximum_iteration = max(self.maximum_iteration, iteration)
        self.precision_preserved &= (
            a.dtype == b.dtype == x.dtype == residual.dtype
            and np.finfo(residual.dtype).nmant + 1 == 64
        )
        self.last = dict(
            a_shape=list(a.shape),
            b_shape=list(b.shape),
            dtype=a.dtype.name,
            native_error=float(error),
            numpy_error=float(numpy_error),
            iteration=iteration,
        )

    def summary(self):
        return dict(
            evaluations=self.evaluations,
            accepted_solves=self.accepted_solves,
            exact_residuals=self.exact_residuals,
            gate_disagreements=self.gate_disagreements,
            maximum_scaled_difference=self.maximum_scaled_difference,
            maximum_accepted_numpy_error=self.maximum_accepted_numpy_error,
            maximum_iteration=self.maximum_iteration,
            precision_preserved=bool(self.precision_preserved),
            last=self.last,
        )
