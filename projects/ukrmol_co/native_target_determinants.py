"""Read-only CAS12 native-CSF expansion into frozen-core determinant amplitudes."""

from __future__ import annotations

import struct
from itertools import combinations
from math import comb

import numpy as np
from scipy.sparse import coo_matrix

from projects.ukrmol_co.native_csf_reconstruction import records as csf_records


def parity(values):
    """Parity of sorting unique nonnegative integer labels, using occupied bits."""
    seen = inversions = 0
    for value in values:
        value = int(value)
        if value < 0 or (seen >> value) & 1:
            raise ValueError("Invalid or repeated occupied spin orbital")
        inversions += (seen >> (value + 1)).bit_count()
        seen |= 1 << value
    return -1 if inversions % 2 else 1


def determinant_index(occupied, spin):
    """Native interleaved addresses -> alpha/beta active masks and sorting sign."""
    occupied = [int(a) for a in occupied]
    if len(occupied) != 14 or any(a < 1 or a > 28 for a in occupied):
        raise ValueError("Only14 electrons in the declared14 spatial orbitals")
    keys = [((a - 1) % 2) * 14 + (a - 1) // 2 for a in occupied]
    sign = parity(keys)
    bits = sum(1 << k for k in keys)
    alpha, beta = bits & ((1 << 14) - 1), bits >> 14
    if (
        alpha & 3 != 3
        or beta & 3 != 3
        or (alpha.bit_count(), beta.bit_count()) != (7 + spin, 7 - spin)
    ):
        raise ValueError("Core occupation or native spin projection differs")
    return alpha >> 2, beta >> 2, sign


def cidata(path):
    """Read all eight original8-root sets, retaining exact source phases."""
    sectors = {}
    with path.open("rb") as stream:

        def record():
            marker = stream.read(4)
            if len(marker) != 4:
                raise ValueError("Missing CIDATA record")
            size = struct.unpack("<I", marker)[0]
            data = stream.read(size)
            if len(data) != size or stream.read(4) != marker:
                raise ValueError("Truncated CIDATA record")
            return data

        for irrep in ["A1", "B1", "B2", "A2"]:
            for spin in [0, 1]:
                if record() != b"********                CIDATA  ":
                    raise ValueError("CIDATA marker")
                h = record()
                number, count = struct.unpack_from("<2q", h)
                centers, dimension, roots, symmetry = struct.unpack_from("<4q", h, 136)
                if not (
                    len(h) == 208
                    and number == len(sectors) + 1
                    and count == 11
                    and centers == 2
                    and roots == 8
                    and symmetry == ["A1", "B1", "B2", "A2"].index(irrep)
                    and struct.unpack_from("<2d", h, 168) == (float(spin), float(spin))
                    and struct.unpack_from("<q", h, 184)[0] == 14
                ):
                    raise ValueError("CIDATA layout")
                charges = [struct.unpack_from("<d", record(), 32)[0] for _ in range(2)]
                if charges != [6.0, 8.0]:
                    raise ValueError("CIDATA nuclei")
                raw = record()
                if len(raw) != 8 * (2 * dimension + 8):
                    raise ValueError("CIDATA mixed record")
                phases = np.frombuffer(raw, "<i8", count=dimension)
                energy = (
                    np.frombuffer(raw, "<f8", count=8, offset=8 * dimension)
                    + struct.unpack_from("<d", h, 192)[0]
                )
                vectors = []
                for root in range(1, 9):
                    raw = record()
                    if len(raw) != 8 * (dimension + 1) or struct.unpack_from("<q", raw)[0] != root:
                        raise ValueError("CIDATA root ordering")
                    vectors.append(np.frombuffer(raw, "<f8", offset=8))
                vectors = np.array(vectors)
                if (
                    not np.isfinite(vectors).all()
                    or not np.isfinite(energy).all()
                    or np.any(np.diff(energy) < 0)
                ):
                    raise ValueError("Invalid native vectors/energies")
                np.testing.assert_allclose(vectors @ vectors.T, np.eye(8), atol=1e-9, rtol=0)
                sectors[(spin, irrep)] = dict(phases=phases, energy=energy, vectors=vectors)
        if stream.read(1):
            raise ValueError("Trailing CIDATA set")
    return sectors


def expand(path, native, spin, irrep):
    """Sparse expansion; never normalize, solve, alter phases or write native files."""
    rows = csf_records(path)
    header = np.frombuffer(rows[0], "<i8", offset=120)
    dimension, electrons, terms, nsym, packed_size = map(int, header[[7, 8, 9, 11, 13]])
    if not (
        len(rows[0]) == 304
        and header[0] == ["A1", "B1", "B2", "A2"].index(irrep)
        and np.array_equal(np.frombuffer(rows[0], "<f8", count=2, offset=128), [spin, spin])
        and electrons == 14
        and header[5] == 14
        and header[6] == 28
        and header[10] == 1
        and nsym == 3
        and header[12] == 2
        and header[21] == dimension
        and header[22] == -1
    ):
        raise ValueError("Only the exact14-orbital CAS12 target-CSF layout is supported")
    phases = np.frombuffer(rows[1], "<i8")
    meta = np.frombuffer(rows[2], "<i8")
    if not np.array_equal(phases, native["phases"]) or not np.isin(phases, [-1, 1]).all():
        raise ValueError("CIDATA/source CSF phases differ")
    if len(meta) != dimension + 14 + 6 * nsym + 1 or not np.array_equal(meta[:3], [8, 3, 3]):
        raise ValueError("CAS12 native orbital inventory")
    reference = meta[3:17].tolist()
    counts = meta[17 : 17 + dimension]
    tail = meta[17 + dimension :]
    if not np.array_equal(tail, np.r_[0, [8, 3, 3], [8, 3, 3], [0, 0, 0], [8, 3, 3], [0, 0, 0]]):
        raise ValueError("CSF target orbital metadata")
    offsets = np.frombuffer(rows[3], "<i8")
    if len(offsets) != 2 * (dimension + 1):
        raise ValueError("CSF offsets")
    co, po = np.split(offsets, 2)
    packed = np.frombuffer(rows[4], "<i8")
    coeff = np.frombuffer(rows[5], "<f8")
    if not (
        len(packed) == packed_size
        and len(coeff) == terms
        and np.array_equal(co, np.r_[1, 1 + np.cumsum(counts)])
        and po[0] == 1
        and po[-1] == packed_size + 1
        and np.all(np.diff(po) > 0)
        and np.all(counts > 0)
        and np.isfinite(coeff).all()
    ):
        raise ValueError("CSF payload lengths")
    np.testing.assert_allclose(np.add.reduceat(coeff**2, co[:-1] - 1), 1, atol=1e-12, rtol=0)
    sizes = [comb(12, 5 + spin), comb(12, 5 - spin)]
    addresses = []
    for count in [5 + spin, 5 - spin]:
        masks = sorted(sum(1 << k for k in occupied) for occupied in combinations(range(12), count))
        addresses.append({mask: index for index, mask in enumerate(masks)})
    det_rows = np.empty(terms, np.int32)
    csf_columns = np.repeat(np.arange(dimension, dtype=np.int32), counts)
    values = np.empty(terms)
    slots = {a: i for i, a in enumerate(reference)}
    if len(slots) != 14 or reference != sorted(reference):
        raise ValueError("Invalid native reference determinant")
    for j, count in enumerate(counts):
        cursor, stop = int(po[j] - 1), int(po[j + 1] - 1)
        seen = set()
        for term in range(count):
            degree = int(packed[cursor])
            end = cursor + 1 + 2 * degree
            if degree < 0 or degree > 14 or end > stop:
                raise ValueError("Packed substitution bounds")
            holes = packed[cursor + 1 : cursor + 1 + degree]
            particles = packed[cursor + 1 + degree : end]
            if len(set(holes)) != degree or any(h not in slots for h in holes):
                raise ValueError("Invalid substitution holes")
            occupied = reference.copy()
            for hole, particle in zip(holes, particles, strict=True):
                occupied[slots[hole]] = int(particle)
            a, b, sign = determinant_index(occupied, spin)
            index = addresses[0][a] * sizes[1] + addresses[1][b]
            if index in seen:
                raise ValueError("Duplicate determinant within CSF")
            seen.add(index)
            offset = int(co[j] - 1 + term)
            if term == 0 and parity(occupied) * (1 if coeff[offset] > 0 else -1) != phases[j]:
                raise ValueError("First-term canonical phase differs")
            det_rows[offset] = index
            values[offset] = sign * coeff[offset]
            cursor = end
        if cursor != stop:
            raise ValueError("Unconsumed packed CSF data")
    transform = coo_matrix(
        (values, (det_rows, csf_columns)), shape=(sizes[0] * sizes[1], dimension)
    ).tocsr()
    result = (transform @ native["vectors"].T).T.reshape(8, *sizes)
    np.testing.assert_allclose(
        (transform.T @ result.reshape(8, -1).T).T, native["vectors"], atol=1e-10, rtol=0
    )
    np.testing.assert_allclose(
        result.reshape(8, -1) @ result.reshape(8, -1).T, np.eye(8), atol=1e-9, rtol=0
    )
    return result, dict(csfs=dimension, determinant_terms=terms, determinant_shape=sizes)
