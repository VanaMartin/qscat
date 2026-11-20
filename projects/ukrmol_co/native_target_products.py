"""Raw primary target roots and explicitly ordered native channel products.

No inner eigenvector recovery, guessed row layout or implicit Q state. Native
signed coupling and phase cancellation: docs/physics/co-native-continuum-product-map.md.
"""

from __future__ import annotations

import struct
from math import comb
from pathlib import Path

import numpy as np

from projects.ukrmol_co.core_sector_full_overlap import TargetProducts
from projects.ukrmol_co.core_sector_overlap import CoreSectorCI, _validate, spin_lower
from projects.ukrmol_co.native_l2_determinants import _record_locations, l2_determinant_transform


def read_primary_cidata(path: Path) -> tuple[dict, ...]:
    """Eight A1/B1/B2/A2 singlet/triplet sets, configuration-by-state raw C.

    Accept at most nine retained states and 50000 configurations per sector.
    No phase application, target selection, sorting, normalization or solve.
    """
    if not 0 < path.stat().st_size <= 64 * 1024**2:
        raise ValueError("Bounded primary CIDATA file required")
    sectors = []
    with path.open("rb") as stream:

        def record():
            marker = stream.read(4)
            if len(marker) != 4:
                raise ValueError("Missing CIDATA marker")
            size = struct.unpack("<I", marker)[0]
            if size > 4 * 1024**2:
                raise ValueError("CIDATA record cap exceeded")
            body = stream.read(size)
            if len(body) != size or stream.read(4) != marker:
                raise ValueError("Truncated CIDATA record")
            return body

        for irrep in range(4):
            for spin in range(2):
                if record() != b"********                CIDATA  ":
                    raise ValueError("CIDATA set marker")
                header = record()
                if len(header) != 208:
                    raise ValueError("CIDATA header length")
                number, count = struct.unpack_from("<2q", header)
                nuclei, dimension, roots, symmetry = struct.unpack_from("<4q", header, 136)
                if not (
                    number == len(sectors) + 1
                    and 2 <= nuclei <= 3
                    and count == nuclei + roots + 1
                    and 1 <= dimension <= 50_000
                    and 1 <= roots <= 9
                    and symmetry == irrep
                    and struct.unpack_from("<2d", header, 168) == (float(spin), float(spin))
                    and struct.unpack_from("<q", header, 184)[0] == 14
                ):
                    raise ValueError("Primary CIDATA sector inventory")
                centers = [record() for _ in range(nuclei)]
                if any(len(c) != 40 for c in centers):
                    raise ValueError("CIDATA center layout")
                charges = [struct.unpack_from("<d", c, 32)[0] for c in centers]
                if charges != [6.0, 8.0] + ([0.0] if nuclei == 3 else []):
                    raise ValueError("CIDATA nuclear inventory")
                mixed = record()
                if len(mixed) != 8 * (2 * dimension + roots):
                    raise ValueError("CIDATA mixed record length")
                phases = np.frombuffer(mixed, "<i8", count=dimension).copy()
                energies = np.frombuffer(mixed, "<f8", count=roots, offset=8 * dimension).copy()
                vectors = []
                for index in range(1, roots + 1):
                    body = record()
                    if (
                        len(body) != 8 * (dimension + 1)
                        or struct.unpack_from("<q", body)[0] != index
                    ):
                        raise ValueError("CIDATA state order or vector length")
                    vectors.append(np.frombuffer(body, "<f8", offset=8).copy())
                c = np.array(vectors).T
                if (
                    not np.isin(phases, [-1, 1]).all()
                    or not np.isfinite(c).all()
                    or not np.isfinite(energies).all()
                ):
                    raise ValueError("Finite raw target data and signed phases required")
                if np.max(abs(c.conj().T @ c - np.eye(roots))) > 1e-9:
                    raise ValueError("Raw target root Gram differs")
                sectors.append(
                    dict(spin=spin, irrep=irrep, phases=phases, coefficients=c, energies=energies)
                )
        if stream.read(1):
            raise ValueError("Trailing CIDATA set")
    return tuple(sectors)


def expand_primary_target(path: Path, coefficients, phases, spin: int) -> CoreSectorCI:
    """Contract all raw13orbital target CSFs with unchanged raw CIDATA columns.

    Reuse the sparse bound-only transformer with no discarded CSF prefix. Never
    apply scattering phase delta to this raw-target representation a second time.
    """
    c, p = np.asarray(coefficients), np.asarray(phases)
    if (
        type(spin) is not int
        or spin not in (0, 1)
        or c.ndim != 2
        or not 1 <= c.shape[0] <= 50_000
        or not 1 <= c.shape[1] <= 9
        or p.shape != (len(c),)
        or not np.isin(p, [-1, 1]).all()
        or not np.isfinite(c).all()
    ):
        raise ValueError("Raw configuration-by-state C, phases and singlet/triplet spin required")
    locations = _record_locations(path)
    with path.open("rb") as stream:
        stream.seek(locations[0][0])
        header = stream.read(locations[0][1])
        if len(header) != 304:
            raise ValueError("Native target CSF header length")
        h = struct.unpack("<qdddd9q6qdqq", header[120:])
        stream.seek(locations[1][0])
        raw_phases = stream.read(locations[1][1])
        stream.seek(locations[2][0])
        inventory = struct.unpack("<3q", stream.read(24))
    if not (
        h[5:9] == (13, 26, len(c), 14)
        and h[11:13] == (3, 2)
        and h[21:23] == (len(c), -1)
        and inventory == (7, 3, 3)
        and len(raw_phases) == 8 * len(c)
        and np.array_equal(np.frombuffer(raw_phases, "<i8"), p)
    ):
        raise ValueError("Raw target CSF inventory/phases differ from CIDATA")
    electrons = (5 + spin, 5 - spin)
    transform = l2_determinant_transform(path, (7, 3, 3), 2, electrons, 0, max_terms=2_000_000)
    roots = np.asarray(transform @ c).T
    if np.max(abs(transform.conj().T @ roots.T - c)) > 1e-9:
        raise ValueError("Target determinant/CSF roundtrip differs")
    if np.max(abs(roots.conj() @ roots.T - np.eye(c.shape[1]))) > 1e-9:
        raise ValueError("Expanded raw target Gram differs")
    return CoreSectorCI(
        2,
        11,
        (0, 1),
        (0, 1),
        electrons,
        roots.reshape(c.shape[1], comb(11, electrons[0]), comb(11, electrons[1])),
    )


def native_target_products(
    sectors: tuple[CoreSectorCI, ...], continuum_indices, amplitudes, ncontinuum: int
) -> tuple[TargetProducts, ...]:
    """Map declared native sector/state/continuum rows into coherent spin families.

    The caller supplies the source-resolved sector order and zero-based spatial
    continuum indices for each sector. Raw targets carry no loaded delta. Output
    factors are singletalpha+1, tripletalpha-1/sqrt3, tripletbeta+sqrt(2/3).
    This accepts a supplied channel vector, never manufactures the missing state.
    """
    if (
        not isinstance(sectors, tuple)
        or not 1 <= len(sectors) <= 8
        or len(continuum_indices) != len(sectors)
    ):
        raise ValueError("One to eight explicitly ordered target sectors required")
    if type(ncontinuum) is not int or not 1 <= ncontinuum <= 160:
        raise ValueError("Bounded common continuum inventory required")
    counts = []
    indices = []
    spins = []
    first = sectors[0]
    for state, selected in zip(sectors, continuum_indices, strict=True):
        c = _validate(state)
        if (
            (state.ncore, state.nactive) != (first.ncore, first.nactive)
            or state.core_alpha != tuple(range(state.ncore))
            or state.core_beta != state.core_alpha
        ):
            raise ValueError("Shared complete paired core and active partition required")
        projection = state.active_electrons[0] - state.active_electrons[1]
        if projection not in (0, 2) or sum(state.active_electrons) != sum(first.active_electrons):
            raise ValueError("SingletM0 or tripletM1 target sectors required")
        if (
            not isinstance(selected, tuple)
            or not selected
            or len(selected) > ncontinuum
            or len(set(selected)) != len(selected)
            or any(type(i) is not int or not 0 <= i < ncontinuum for i in selected)
        ):
            raise ValueError("Explicit unique continuum indices required")
        counts.append(len(c))
        indices.append(selected)
        spins.append(projection // 2)
    if sum(counts) > 64:
        raise ValueError("Combined target root cap exceeded")
    rows = np.asarray(amplitudes)
    if (
        rows.shape != (sum(n * len(i) for n, i in zip(counts, indices, strict=True)),)
        or not np.isfinite(rows).all()
    ):
        raise ValueError("Finite complete native sector/state/continuum row vector required")
    zero = []
    one = []
    alpha = []
    beta = []
    offset = 0
    for state, selected, spin, count in zip(sectors, indices, spins, counts, strict=True):
        size = count * len(selected)
        values = np.zeros((count, ncontinuum), dtype=rows.dtype)
        values[:, selected] = rows[offset : offset + size].reshape(count, len(selected))
        offset += size
        if spin:
            lowered = spin_lower(state)
            if lowered is None:
                raise ValueError("Triplet target has no spin-lowered component")
            zero.append(lowered.coefficients / np.sqrt(2))
            one.append(state.coefficients)
            alpha.append(-values / np.sqrt(3))
            beta.append(values * np.sqrt(2 / 3))
        else:
            zero.append(state.coefficients)
            alpha.append(values)
    core = tuple(range(first.ncore))
    total = sum(first.active_electrons)
    if total % 2:
        raise ValueError("Even target electron count required")
    m0 = CoreSectorCI(
        first.ncore, first.nactive, core, core, (total // 2,) * 2, np.concatenate(zero)
    )
    products = [TargetProducts(m0, 0, np.concatenate(alpha))]
    if one:
        m1 = CoreSectorCI(
            first.ncore,
            first.nactive,
            core,
            core,
            (total // 2 + 1, total // 2 - 1),
            np.concatenate(one),
        )
        products.append(TargetProducts(m1, 1, np.concatenate(beta)))
    return tuple(products)
