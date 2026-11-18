"""Decode only the uncontracted bound-orbital tail of native scattering CSFs.

The continuum contraction is deliberately unsupported. Its omitted norm belongs
in a projection bound, not in an invented determinant expansion.
"""

from __future__ import annotations

import struct
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from projects.ukrmol_co.native_target_determinants import parity


def _record_locations(path):
    locations = []
    with path.open("rb") as stream:
        size = path.stat().st_size
        while marker := stream.read(4):
            if len(marker) != 4:
                raise ValueError("Truncated CSF marker")
            length = struct.unpack("<I", marker)[0]
            if length > size - stream.tell() - 4:
                raise ValueError("CSF record extends beyond EOF")
            locations.append((stream.tell(), length))
            stream.seek(length, 1)
            if stream.read(4) != marker:
                raise ValueError("Mismatched CSF trailer")
    if len(locations) != 6:
        raise ValueError("Expected six native CSF records")
    return locations


def l2_determinant_transform(
    path: Path,
    bound_orbitals: tuple[int, ...],
    ncore: int,
    active_electrons: tuple[int, int],
    last_continuum_csf: int,
    *,
    max_terms: int = 2_000_000,
) -> csr_matrix:
    """Map trailing L2 CSFs into core-factored alpha/beta bitstring determinants.

    Bound orbitals are first within each native symmetry; frozen core orbitals
    are first in the concatenated bound inventory. The caller establishes the
    tail boundary/order from the pinned native contraction source and metadata.
    Packed records are memory mapped; no full continuum CSF expansion is built.
    """
    locations = _record_locations(path)
    with path.open("rb") as stream:
        stream.seek(locations[0][0])
        header = stream.read(locations[0][1])
    if len(header) != 304:
        raise ValueError("Unsupported native CSF header")
    h = struct.unpack("<qdddd9q6qdqq", header[120:])
    orbitals, spin_orbitals, dimension, electrons, terms, nsym, packed_size = (
        int(h[k]) for k in [5, 6, 7, 8, 9, 11, 13]
    )
    bound = sum(bound_orbitals)
    if (
        len(bound_orbitals) != nsym
        or any(type(n) is not int or n < 0 for n in bound_orbitals)
        or type(ncore) is not int
        or not 0 <= ncore < bound
        or len(active_electrons) != 2
        or any(type(n) is not int or not 0 <= n <= bound - ncore for n in active_electrons)
        or electrons != 2 * ncore + sum(active_electrons)
        or h[1:3] != (0.5 * (active_electrons[0] - active_electrons[1]),) * 2
        or spin_orbitals != 2 * orbitals
        or type(last_continuum_csf) is not int
        or not 0 <= last_continuum_csf < dimension
        or dimension - last_continuum_csf > 50_000
        or type(max_terms) is not int
        or not 1 <= max_terms <= 10_000_000
    ):
        raise ValueError("Unsupported bound inventory, spin, core or L2 tail")

    def mapped(record, dtype):
        offset, length = locations[record]
        if length % 8:
            raise ValueError("Non-eight-byte native payload")
        return np.memmap(path, mode="r", dtype=dtype, offset=offset, shape=(length // 8,))

    meta, offsets, packed, coefficients = (
        mapped(i, dtype) for i, dtype in [(2, "<i8"), (3, "<i8"), (4, "<i8"), (5, "<f8")]
    )
    if len(meta) < nsym + electrons + dimension:
        raise ValueError("Incomplete CSF metadata")
    inventory = meta[:nsym]
    reference = meta[nsym : nsym + electrons].tolist()
    counts = meta[nsym + electrons : nsym + electrons + dimension]
    if (
        sum(inventory) != orbitals
        or any(b > n for b, n in zip(bound_orbitals, inventory, strict=True))
        or len(offsets) != 2 * (dimension + 1)
        or len(packed) != packed_size
        or len(coefficients) != terms
        or len(set(reference)) != electrons
    ):
        raise ValueError("CSF record inventories disagree")
    co, po = np.split(offsets, 2)
    if (
        co[0] != 1
        or co[-1] != terms + 1
        or po[0] != 1
        or po[-1] != packed_size + 1
        or np.any(counts <= 0)
        or not np.array_equal(np.diff(co), counts)
        or np.any(np.diff(po) <= 0)
        or sum(counts[last_continuum_csf:]) > max_terms
    ):
        raise ValueError("Invalid native offsets or L2 term cap")
    orbital_map = {}
    native_offset = bound_offset = 0
    for total, kept in zip(inventory, bound_orbitals, strict=True):
        orbital_map.update({native_offset + i: bound_offset + i for i in range(kept)})
        native_offset += int(total)
        bound_offset += kept
    ncas = bound - ncore
    sizes = [comb(ncas, n) for n in active_electrons]
    if sizes[0] * sizes[1] > 1_000_000:
        raise ValueError("Active determinant dimension cap exceeded")
    addresses = [
        {
            mask: i
            for i, mask in enumerate(
                sorted(sum(1 << j for j in c) for c in combinations(range(ncas), n))
            )
        }
        for n in active_electrons
    ]
    core_mask = (1 << ncore) - 1
    rows, columns, values = [], [], []
    for column, csf in enumerate(range(last_continuum_csf, dimension)):
        cursor, stop = int(po[csf] - 1), int(po[csf + 1] - 1)
        seen = set()
        norm = 0.0
        for term in range(int(counts[csf])):
            degree = int(packed[cursor])
            end = cursor + 1 + 2 * degree
            if not 0 <= degree <= electrons or end > stop:
                raise ValueError("Packed L2 substitution bounds")
            holes = packed[cursor + 1 : cursor + degree + 1].tolist()
            particles = packed[cursor + degree + 1 : end].tolist()
            if len(set(holes)) != degree or any(hole not in reference for hole in holes):
                raise ValueError("Invalid L2 holes")
            occupied = reference.copy()
            for hole, particle in zip(holes, particles, strict=True):
                occupied[reference.index(hole)] = particle
            keys = []
            for address in occupied:
                spatial, spin = divmod(address - 1, 2)
                if spatial not in orbital_map:
                    raise ValueError("Continuum or invalid orbital in declared L2 tail")
                keys.append(spin * bound + orbital_map[spatial])
            sign = parity(keys)
            alpha = sum(1 << k for k in keys if k < bound)
            beta = sum(1 << (k - bound) for k in keys if k >= bound)
            if (
                alpha & core_mask != core_mask
                or beta & core_mask != core_mask
                or (alpha.bit_count(), beta.bit_count())
                != tuple(n + ncore for n in active_electrons)
            ):
                raise ValueError("L2 core occupation or spin projection differs")
            row = addresses[0][alpha >> ncore] * sizes[1] + addresses[1][beta >> ncore]
            if row in seen:
                raise ValueError("Duplicate determinant in L2 CSF")
            seen.add(row)
            value = float(coefficients[int(co[csf] - 1) + term])
            if not np.isfinite(value):
                raise ValueError("Nonfinite L2 coefficient")
            norm += value * value
            rows.append(row)
            columns.append(column)
            values.append(sign * value)
            cursor = end
        if cursor != stop or abs(norm - 1) > 1e-10:
            raise ValueError("Unconsumed or nonunit L2 CSF")
    return coo_matrix(
        (values, (rows, columns)), shape=(sizes[0] * sizes[1], dimension - last_continuum_csf)
    ).tocsr()
