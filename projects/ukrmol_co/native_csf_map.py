"""Prove a signed diagonal correspondence between saved target-CSF address layouts.

This checks representations only; it never writes converted native target vectors.
See docs/physics/co-native-target-address-map.md for the phase-aware native route.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


@dataclass
class TargetCSFs:
    """The six-record electron-only C2v target format, with one-based offsets."""

    inventory: np.ndarray
    reference: np.ndarray
    counts: np.ndarray
    phases: np.ndarray
    coefficient_offsets: np.ndarray
    packed_offsets: np.ndarray
    packed: np.ndarray
    coefficients: np.ndarray
    symmetry: int
    spin: float


def sorting_phase(ordered):
    """Fermionic sign of sorting electron slots, without changing the input."""
    if len(set(ordered)) != len(ordered):
        raise ValueError("Repeated occupied spin orbital")
    inversions = sum(a > b for i, a in enumerate(ordered) for b in ordered[i + 1 :])
    return -1 if inversions % 2 else 1


def address_map(inventory, target=(7, 3, 3, 0)):
    """Map native spin addresses to common target addresses; zero marks continuum."""
    if len(inventory) > len(target) or any(n < 0 for n in inventory):
        raise ValueError("Invalid C2v orbital inventory")
    sizes = list(inventory) + [0] * (len(target) - len(inventory))
    if any(n < t for n, t in zip(sizes, target, strict=True)):
        raise ValueError("Target orbitals missing from native layout")
    mapping = np.zeros(2 * sum(sizes) + 1, dtype=np.int64)
    native = common = 0
    for n, t in zip(sizes, target, strict=True):
        mapping[2 * native + 1 : 2 * (native + t) + 1] = np.arange(
            2 * common + 1, 2 * (common + t) + 1
        )
        native += n
        common += t
    return mapping


def read_target_csfs(path: Path) -> TargetCSFs:
    """Read the pinned int64/float64 target format and require complete EOF."""
    with FortranFile(path, "r") as stream:
        records = [stream.read_record("u1").tobytes() for _ in range(6)]
        if stream._fp.tell() != path.stat().st_size:
            raise ValueError("Trailing target-CSF records")
    if len(records[0]) != 304 or any(len(r) % 8 for r in records):
        raise ValueError("Unknown target-CSF native word layout")
    header = np.frombuffer(records[0], "<i8", offset=120)
    spin = np.frombuffer(records[0], "<f8", count=2, offset=128)
    norb, nsrb, dimension, electrons, terms, idiag, nsym, symtyp, packed_size = header[5:14]
    if not (
        electrons == 14
        and idiag == 1
        and nsym in (3, 4)
        and symtyp == 2
        and nsrb == 2 * norb
        and header[21] == dimension
        and header[22] == -1
        and dimension > 0
        and terms > 0
        and packed_size > 0
        and spin[0] in (0, 1)
        and spin[0] == spin[1]
        and header[0] in range(4)
    ):
        raise ValueError("Not the supported C2v CO target-CSF header")
    phases = np.frombuffer(records[1], "<i8")
    metadata = np.frombuffer(records[2], "<i8")
    indices = np.frombuffer(records[3], "<i8")
    packed = np.frombuffer(records[4], "<i8")
    coefficients = np.frombuffer(records[5], "<f8")
    # wrnfto: nob, ndtrf, nodo, iposit, nob0, nobl(2*nsym), nob0l(2*nsym).
    if (len(phases), len(metadata), len(indices), len(packed), len(coefficients)) != (
        dimension,
        dimension + electrons + 6 * nsym + 1,
        2 * (dimension + 1),
        packed_size,
        terms,
    ):
        raise ValueError("Target-CSF payload dimensions disagree")
    inventory = metadata[:nsym]
    reference = metadata[nsym : nsym + electrons]
    counts = metadata[nsym + electrons : nsym + electrons + dimension]
    end = nsym + electrons + dimension
    tail = metadata[end:]
    target = np.array([7, 3, 3, 0][:nsym])
    if not (
        sum(inventory) == norb
        and tail[0] == 0
        and np.array_equal(tail[1 : 1 + nsym], target)
        and np.array_equal(tail[1 + nsym : 1 + 3 * nsym], np.r_[inventory, np.zeros(nsym)])
        and np.array_equal(tail[1 + 3 * nsym :], np.r_[target, np.zeros(nsym)])
        and np.all(np.isin(phases, [-1, 1]))
        and np.all(counts > 0)
        and np.all(np.isfinite(coefficients))
        and np.all(coefficients != 0)
    ):
        raise ValueError("Invalid target-CSF orbital/phase/count metadata")
    coefficient_offsets, packed_offsets = np.split(indices, 2)
    if not (
        np.array_equal(coefficient_offsets, np.r_[1, 1 + np.cumsum(counts)])
        and packed_offsets[0] == 1
        and packed_offsets[-1] == packed_size + 1
        and np.all(np.diff(packed_offsets) > 0)
    ):
        raise ValueError("Invalid target-CSF offsets")
    norms = np.add.reduceat(coefficients**2, coefficient_offsets[:-1] - 1)
    np.testing.assert_allclose(norms, 1, atol=1e-12, rtol=0)
    return TargetCSFs(
        inventory,
        reference,
        counts,
        phases,
        coefficient_offsets,
        packed_offsets,
        packed,
        coefficients,
        int(header[0]),
        float(spin[0]),
    )


def canonical_packed(csf: TargetCSFs):
    """Validate every substitution, relabel addresses, and reconstruct first-term phases."""
    mapping = address_map(csf.inventory).tolist()

    def map_addresses(values):
        if any(a <= 0 or a >= len(mapping) for a in values):
            raise ValueError("Invalid occupied spin-orbital address")
        mapped = [mapping[a] for a in values]
        if any(a == 0 for a in mapped):
            raise ValueError("Continuum occupied in a target CSF")
        return mapped

    reference = map_addresses(csf.reference)
    if np.any(np.diff(reference) <= 0):
        raise ValueError("Reference determinant is not strictly ordered")
    packed = csf.packed.tolist()
    mapped = packed.copy()
    reference_slots = {int(a): i for i, a in enumerate(reference)}
    for j, count in enumerate(csf.counts):
        cursor = int(csf.packed_offsets[j] - 1)
        stop = int(csf.packed_offsets[j + 1] - 1)
        for term in range(count):
            if cursor >= stop:
                raise ValueError("Missing packed determinant")
            degree = packed[cursor]
            end = cursor + 1 + 2 * degree
            if degree < 0 or degree > len(reference) or end > stop:
                raise ValueError("Invalid packed substitution degree")
            addresses = map_addresses(packed[cursor + 1 : end])
            holes, particles = addresses[:degree], addresses[degree:]
            if len(set(holes)) != degree or any(int(h) not in reference_slots for h in holes):
                raise ValueError("Substitution does not remove distinct reference slots")
            occupied = reference.copy()
            for h, p in zip(holes, particles, strict=True):
                occupied[reference_slots[int(h)]] = p
            if len(set(occupied)) != len(reference):
                raise ValueError("Substitution creates repeated occupied spin orbitals")
            if sum(a % 2 for a in occupied) != 7 + int(csf.spin):
                raise ValueError("Substitution changes target spin projection")
            mapped[cursor + 1 : end] = addresses
            if term == 0:
                coefficient = csf.coefficients[csf.coefficient_offsets[j] - 1]
                phase = sorting_phase(occupied) * (1 if coefficient > 0 else -1)
                if phase != csf.phases[j]:
                    raise ValueError("Canonical phase disagrees with first determinant")
            cursor = end
        if cursor != stop:
            raise ValueError("Trailing packed determinant data within CSF")
    return np.asarray(reference), np.asarray(mapped)


def signed_correspondence(source: TargetCSFs, destination: TargetCSFs):
    """Return Q's diagonal only after every reference, packed term and sign agrees."""
    if (source.symmetry, source.spin) != (destination.symmetry, destination.spin):
        raise ValueError("Different target-CSF sectors")
    for name in ("counts", "coefficient_offsets", "packed_offsets"):
        if not np.array_equal(getattr(source, name), getattr(destination, name)):
            raise ValueError(f"Different target-CSF {name}; no diagonal correspondence")
    for left, right in zip(canonical_packed(source), canonical_packed(destination), strict=True):
        if not np.array_equal(left, right):
            raise ValueError("Target determinant addresses differ after embedding")
    signs = source.phases * destination.phases
    expected = source.coefficients * np.repeat(signs, source.counts)
    if expected.tobytes() != destination.coefficients.tobytes():
        raise ValueError("Nonleading determinant coefficient violates the CSF phase rule")
    return signs
