"""Compatibility shim: the ansatz was promoted to `qscat.model.flexible`."""

from qscat.model.flexible import (
    FlexibleDiatomicModel,
    SmoothR,
    TailR,
    from_diatomic,
    pack,
    params,
    unpack,
    with_params,
    y_p,
)

__all__ = [
    "FlexibleDiatomicModel",
    "SmoothR",
    "TailR",
    "from_diatomic",
    "pack",
    "params",
    "unpack",
    "with_params",
    "y_p",
]
