"""QSCAT standard library — validated quantum-scattering numerics.

The single source of truth for the package version is ``__version__`` below;
``pyproject.toml`` reads it dynamically (``[tool.hatch.version]``). Bump it here
only, then tag the release ``qscat-v<version>``.
"""

__version__ = "0.1.0.dev0"

from typing import TYPE_CHECKING

from qscat.exceptions import (
    BackendError,
    ConvergenceError,
    GridError,
    ModelError,
    QscatError,
)

# Submodules are exposed lazily (PEP 562): `qscat.dvr` triggers the import on
# first access, so `import qscat` is cheap and pulls in numpy/scipy only for the
# submodules a caller actually touches. Mirrors SciPy's namespace ergonomics.
_SUBMODULES = (
    "units",
    "linalg",
    "special",
    "dvr",
    "ecs",
    "evolution",
    "core",
    "model",
    "tuning",
)

if TYPE_CHECKING:  # let type-checkers/IDEs resolve `qscat.<submodule>` statically
    from qscat import (  # noqa: F401
        core,
        dvr,
        ecs,
        evolution,
        linalg,
        model,
        special,
        tuning,
        units,
    )


def __getattr__(name: str) -> object:
    if name in _SUBMODULES:
        import importlib

        module = importlib.import_module(f"qscat.{name}")
        globals()[name] = module  # cache so subsequent access skips __getattr__
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), *_SUBMODULES})


__all__ = [
    "__version__",
    "QscatError",
    "GridError",
    "ModelError",
    "BackendError",
    "ConvergenceError",
    *_SUBMODULES,
]
