"""The typed exception hierarchy and the lazy top-level namespace (Phase B)."""

from __future__ import annotations

import warnings

import pytest
import qscat
from qscat.exceptions import (
    BackendError,
    ConvergenceError,
    GridError,
    ModelError,
    QscatError,
)


def test_exception_hierarchy_and_builtin_compat() -> None:
    # Everything is a QscatError...
    for exc in (GridError, ModelError, BackendError, ConvergenceError):
        assert issubclass(exc, QscatError)
    # ...and each also IS the builtin it replaces, so `except <builtin>` still
    # catches it (adopting the hierarchy is backward compatible).
    assert issubclass(GridError, ValueError)
    assert issubclass(ModelError, ValueError)
    assert issubclass(BackendError, RuntimeError)
    assert issubclass(ConvergenceError, RuntimeError)


def test_grid_error_raised_and_catchable_both_ways() -> None:
    from qscat.core.grids import electronic_grid

    # r_max below the fixed inner segments is an invalid grid spec.
    with pytest.raises(GridError):
        electronic_grid(r_max=0.1, order=6, n_complex=3)
    with pytest.raises(ValueError):  # still a ValueError
        electronic_grid(r_max=0.1, order=6, n_complex=3)
    with pytest.raises(QscatError):  # and a QscatError
        electronic_grid(r_max=0.1, order=6, n_complex=3)


def test_lazy_namespace_exposes_submodules() -> None:
    # Attribute access triggers the import (PEP 562).
    assert qscat.dvr.__name__ == "qscat.dvr"
    assert qscat.linalg.__name__ == "qscat.linalg"
    for name in ("core", "model", "tuning", "special", "ecs", "evolution", "units"):
        assert name in dir(qscat)
    with pytest.raises(AttributeError):
        _ = qscat.not_a_submodule


def test_lazy_submodule_import_does_not_warn() -> None:
    # `qscat.__getattr__` is a lazy submodule importer (PEP 562), not a
    # deprecation shim, and must resolve submodules silently.
    #
    # Call `qscat.__getattr__` directly (rather than `qscat.dvr`) so the test
    # still exercises the hook even if some other test already imported and
    # cached the submodule as a real module attribute.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        module = qscat.__getattr__("dvr")
    assert module.__name__ == "qscat.dvr"
