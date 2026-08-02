"""Typed exception hierarchy for qscat.

Every qscat-specific error derives from `QscatError`, so a caller can catch
*any* qscat failure with a single `except QscatError`. Each concrete class also
subclasses the built-in it semantically replaces (`GridError` is a `ValueError`,
`BackendError` is a `RuntimeError`, …), so existing `except ValueError` /
`pytest.raises(ValueError)` code keeps working unchanged — adopting the hierarchy
is backward compatible.

Scope: these classes name the major *recoverable* categories a user might want
to branch on (a bad grid spec vs. an unavailable solver backend vs. a
non-converged iteration). Generic argument validation that a user cannot
meaningfully recover from may still raise a plain built-in `ValueError`/
`TypeError` — the hierarchy is for the categories worth catching, not a blanket
wrapper around every raise.
"""

from __future__ import annotations

__all__ = [
    "QscatError",
    "GridError",
    "ModelError",
    "BackendError",
    "ConvergenceError",
]


class QscatError(Exception):
    """Base class for all errors raised by qscat."""


class GridError(QscatError, ValueError):
    """Invalid FEM-DVR-ECS grid / discretisation specification.

    Raised when a grid, element, or quadrature specification is malformed or
    out of range (e.g. a non-increasing segment endpoint, a quadrature order
    below 2, an `r_max` inside the fixed inner segments).

    Examples
    --------
    Each qscat error subclasses the built-in it replaces, so existing
    ``except``/``pytest.raises`` code keeps working:

    >>> from qscat.exceptions import GridError, QscatError
    >>> issubclass(GridError, ValueError)
    True
    >>> issubclass(GridError, QscatError)
    True
    """


class ModelError(QscatError, ValueError):
    """Invalid model specification, parameters, or registry lookup."""


class BackendError(QscatError, RuntimeError):
    """A linear-algebra backend is unavailable or failed at runtime.

    Raised when a requested sparse backend (e.g. MUMPS) is not installed, or a
    backend operation fails in a way that is about the backend rather than the
    input (e.g. an optional dependency missing on the chosen code path).
    """


class ConvergenceError(QscatError, RuntimeError):
    """A numerical procedure failed to reach its convergence target."""
