"""Bounded selected CONGEN determinant decoding, not a continuum-product adapter.

Only the electron-only D2h int64/float64 six-record format is admitted. See
docs/physics/co-native-continuum-product-map.md for phase and coupling gates.
"""

from __future__ import annotations

import struct
from itertools import islice
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_l2_determinants import _record_locations
from projects.ukrmol_co.native_target_determinants import parity


def selected_csf_determinants(path: Path, indices, *, max_terms: int = 100_000) -> dict:
    """Decode selected zero-based CSFs into canonical alpha-then-beta determinants.

    Return source inventories, target-phase/group metadata and selected sparse
    amplitudes. No target phase correction or CG factor is applied implicitly.
    Large payloads remain read-only memory maps; unselected CSFs are not expanded.
    """
    indices = tuple(islice(iter(indices), 33))
    if (
        not 1 <= len(indices) <= 32
        or len(set(indices)) != len(indices)
        or any(type(i) is not int or i < 0 for i in indices)
        or type(max_terms) is not int
        or not 1 <= max_terms <= 100_000
    ):
        raise ValueError("One to32 unique CSF indices and bounded term cap required")
    locations = _record_locations(path)
    if locations[0][1] != 304:
        raise ValueError("Unsupported native CSF header")
    with path.open("rb") as stream:
        stream.seek(locations[0][0])
        header = stream.read(locations[0][1])
    h = struct.unpack("<qdddd9q6qdqq", header[120:])
    norb, nso, ncsf, nelec, terms, nsym, packed_size = (int(h[k]) for k in [5, 6, 7, 8, 9, 11, 13])
    nphase, ngroups = int(h[21]), int(h[22])
    if (
        not 1 <= norb <= 256
        or nso != 2 * norb
        or not 1 <= nelec <= 16
        or not 1 <= ncsf <= 1_000_000
        or max(indices) >= ncsf
        or not 1 <= nsym <= 8
        or h[12] != 2
        or not np.isfinite(h[1:5]).all()
        or not 0 < nphase <= ncsf
        or not -1 <= ngroups <= 16
    ):
        raise ValueError("Electron-only bounded D2h CSF layout required")

    def mapped(record, dtype):
        offset, length = locations[record]
        if length % 8:
            raise ValueError("Non-eight-byte CSF record")
        return np.memmap(path, mode="r", dtype=dtype, offset=offset, shape=(length // 8,))

    phase_data, meta, offsets, packed, coefficients = (
        mapped(i, dtype)
        for i, dtype in [(1, "<i8"), (2, "<i8"), (3, "<i8"), (4, "<i8"), (5, "<f8")]
    )
    if (
        len(phase_data) != nphase + 4 * max(ngroups, 0)
        or len(meta) != nsym + nelec + ncsf + 1 + 5 * nsym
        or len(offsets) != 2 * (ncsf + 1)
        or len(packed) != packed_size
        or len(coefficients) != terms
        or not np.isin(phase_data[:nphase], [-1, 1]).all()
    ):
        raise ValueError("CSF record inventories disagree")
    inventory = tuple(map(int, meta[:nsym]))
    reference = tuple(map(int, meta[nsym : nsym + nelec]))
    counts = meta[nsym + nelec : nsym + nelec + ncsf]
    tail = meta[nsym + nelec + ncsf :]
    co, po = np.split(offsets, 2)
    if (
        sum(inventory) != norb
        or any(i < 0 for i in inventory)
        or len(set(reference)) != nelec
        or any(not 1 <= i <= nso for i in reference)
        or tail[0] != 0
        or co[0] != 1
        or co[-1] != terms + 1
        or po[0] != 1
        or po[-1] != packed_size + 1
        or np.any(counts <= 0)
        or not np.array_equal(np.diff(co), counts)
        or np.any(np.diff(po) <= 0)
        or sum(int(counts[i]) for i in indices) > max_terms
    ):
        raise ValueError("Invalid occupations/offsets or selected determinant term cap")
    groups = tuple(tuple(map(int, values)) for values in np.split(phase_data[nphase:], 4))
    if ngroups > 0 and (sum(groups[0]) != nphase or any(i <= 0 for i in groups[0])):
        raise ValueError("Target prototype group sizes disagree")
    selected = {}
    for csf in indices:
        cursor, stop = int(po[csf] - 1), int(po[csf + 1] - 1)
        determinant = {}
        for term in range(int(counts[csf])):
            degree = int(packed[cursor])
            end = cursor + 1 + 2 * degree
            if not 0 <= degree <= nelec or end > stop:
                raise ValueError("Selected packed determinant bounds")
            holes = tuple(map(int, packed[cursor + 1 : cursor + 1 + degree]))
            particles = tuple(map(int, packed[cursor + 1 + degree : end]))
            if len(set(holes)) != degree or any(i not in reference for i in holes):
                raise ValueError("Selected determinant has invalid substitution holes")
            occupied = list(reference)
            for hole, particle in zip(holes, particles, strict=True):
                occupied[reference.index(hole)] = particle
            if len(set(occupied)) != nelec or any(not 1 <= i <= nso for i in occupied):
                raise ValueError("Selected determinant has invalid/repeated spin orbitals")
            keys = [((i - 1) % 2) * norb + (i - 1) // 2 for i in occupied]
            if 2 * sum(i < norb for i in keys) - nelec != 2 * h[2]:
                raise ValueError("Selected determinant spin projection disagrees")
            value = float(coefficients[int(co[csf] - 1) + term]) * parity(keys)
            key = tuple(sorted(keys))
            if key in determinant or not np.isfinite(value):
                raise ValueError("Duplicate determinant or nonfinite CSF coefficient")
            determinant[key] = value
            cursor = end
        if cursor != stop or abs(sum(v * v for v in determinant.values()) - 1) > 1e-10:
            raise ValueError("Unconsumed or nonunit selected CSF")
        selected[csf] = determinant
    return dict(
        norb=norb,
        electrons=nelec,
        spin=float(h[1]),
        projection=float(h[2]),
        csfs=ncsf,
        inventory=inventory,
        phase_labels=phase_data[:nphase].copy(),
        target_groups=groups,
        selected=selected,
    )


def corrected_target_coefficients(coefficients, target_phases, scattering_phases) -> np.ndarray:
    """Apply the explicit native row-phase product; do not normalize or infer signs."""
    c, p, q = map(np.asarray, (coefficients, target_phases, scattering_phases))
    if (
        c.ndim != 2
        or c.size > 1_000_000
        or p.shape != (len(c),)
        or q.shape != p.shape
        or not np.isfinite(c).all()
        or not np.isin(p, [-1, 1]).all()
        or not np.isin(q, [-1, 1]).all()
    ):
        raise ValueError(
            "Finite configuration-by-state coefficients and explicit signed phases required"
        )
    return (p * q)[:, None] * c
