"""Shared `.npz` round-trip for `qscat.core`'s frozen dataclass result holders.

`ResonanceLevels.save` (`qscat.core.lcp.levels`) and `ExactResonanceStates.save`
(`qscat.core.resonance`) had byte-identical bodies: write every dataclass field
to an `.npz` archive under the field's own name. Both are caches for
expensive results -- a 2-D pole search is minutes to tens of minutes of
sparse factorizations -- and a hand-rolled cache once drifted: one call site
stored `res_el`/`res_nuc` where `ExactResonanceStates` calls its fields
`residual_electronic`/`residual_nuclear`, a rename away from silently loading
garbage. Round-tripping through the dataclass's own field names keeps them
the dataclass's business; this module is the one copy of that mechanism.

`load` is shared only as far as reading the archive back into a
field-name -> array dict, with the missing-field check. Each holder keeps its
own validation after that (`ExactResonanceStates.load` additionally rejects
an archive written before its row-per-state layout flip) -- that guard is
deliberately NOT here, so it cannot leak onto a holder that never had it.

Private (leading underscore): this is implementation plumbing, not part of
the public API surface.
"""

from __future__ import annotations

import os
from dataclasses import fields
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from _typeshed import DataclassInstance


def save_dataclass_npz(obj: DataclassInstance, path: str | os.PathLike[str]) -> None:
    """Write every dataclass field of `obj` to an `.npz` archive under its own name."""
    np.savez(path, **{f.name: getattr(obj, f.name) for f in fields(obj)})


def load_dataclass_npz(
    cls: type[DataclassInstance], path: str | os.PathLike[str]
) -> dict[str, npt.NDArray[Any]]:
    """Read a `save_dataclass_npz` archive into a field-name -> array dict.

    Raises `ValueError` if any of `cls`'s fields are missing from the
    archive. Does not construct `cls` -- the caller applies its own
    validation (if any) and does that itself.
    """
    with np.load(path) as z:
        missing = [f.name for f in fields(cls) if f.name not in z]
        if missing:
            article = "an" if cls.__name__[0] in "AEIOU" else "a"
            raise ValueError(f"{path} is not {article} {cls.__name__} archive: missing {missing}")
        return {f.name: z[f.name] for f in fields(cls)}
