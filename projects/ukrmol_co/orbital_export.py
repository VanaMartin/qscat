"""Read the bounded saved-basis export; stored overlaps are not a new cross metric."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from projects.ukrmol_co.saved_matrices import records


def read_orbital_export(path: Path) -> dict:
    """Validate complete little-endian double/int64 CGTO export framing and layout.

    Coefficients and shell normalization are returned unchanged. This reader
    performs no integral evaluation, orbital renormalization or qualification
    of a different geometry/domain.
    """
    if not 0 < path.stat().st_size <= 64 * 1024**2:
        raise ValueError("Saved orbital export file-size cap")
    raw = records(path)
    if len(raw) < 11 or raw[0] != b"CO_ORBITAL_EXPORT_V1".ljust(32):
        raise ValueError("Unsupported orbital export version")

    def array(record, dtype, shape):
        if len(record) != 8 * int(np.prod(shape)):
            raise ValueError("Orbital export record length differs from inventory")
        result = np.frombuffer(record, dtype).reshape(shape, order="F").copy()
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite orbital export payload")
        return result

    kind, nao, nmo, nuclei, symmetries, nshell = map(int, array(raw[1], "<i8", (6,)))
    if kind not in (0, 1) or not (
        0 < nao <= 10000
        and 0 < nmo <= 2000
        and 0 <= nuclei <= 64
        and 0 < symmetries <= 8
        and 0 < nshell <= 10000
    ):
        raise ValueError("Orbital export dimension/kind cap")
    if len(raw) != 11 + 5 * nshell:
        raise ValueError("Missing shell records or trailing orbital export data")
    radius = float(array(raw[2], "<f8", (1,))[0])
    counts = array(raw[3], "<i8", (symmetries,))
    mapping = array(raw[4], "<i8", (3, nmo))
    flags = array(raw[5], "<i8", (nmo,))
    coefficients = array(raw[6], "<f8", (nao, nmo))
    overlap = array(raw[7], "<f8", (nmo, nmo))
    ao_map = array(raw[8], "<i8", (2, nao))
    descriptors = array(raw[9], "<i8", (6, nshell))
    geometry = array(raw[10], "<f8", (4, nuclei))
    if np.any(counts < 0) or counts.sum() != nmo:
        raise ValueError("MO symmetry counts differ from the bounded inventory")
    wanted = [(g + 1, j + 1) for g, count in enumerate(counts) for j in range(count)]
    if (
        radius <= 0
        or np.any(counts < 0)
        or counts.sum() != nmo
        or list(zip(mapping[1], mapping[2], strict=True)) != wanted
        or set(mapping[0]) != set(range(1, nmo + 1))
        or not np.isin(flags, [0, 1]).all()
        or not np.all(descriptors[0] == 1)
        or set(descriptors[1]) != set(range(1, nshell + 1))
        or not np.isin(descriptors[2], [0, 1]).all()
    ):
        raise ValueError("MO/AO symmetry, continuum or shell inventory differs")
    shells = []
    for index in range(nshell):
        offset = 11 + 5 * index
        angular, primitives, continuum = map(int, array(raw[offset], "<i8", (3,)))
        if not 0 <= angular <= 20 or not 1 <= primitives <= 1000 or continuum not in (0, 1):
            raise ValueError("CGTO angular/primitive inventory cap")
        location = array(raw[offset + 1], "<f8", (4,))
        exponent, contraction, norm = [
            array(raw[offset + j], "<f8", (primitives,)) for j in (2, 3, 4)
        ]
        if np.any(exponent <= 0):
            raise ValueError("Nonpositive CGTO exponent")
        shells.append(
            dict(
                angular=angular,
                continuum=continuum,
                center=location[:3],
                normalization=location[3],
                exponents=exponent,
                contractions=contraction,
                primitive_norms=norm,
            )
        )
    covered = np.zeros(nao, int)
    for column in descriptors.T:
        _, shell_index, continuum, first, functions, angular = map(int, column)
        shell = shells[shell_index - 1]
        if not (
            1 <= first <= nao
            and functions == 2 * angular + 1
            and first + functions - 1 <= nao
            and angular == shell["angular"]
            and continuum == shell["continuum"]
        ):
            raise ValueError("CGTO descriptor and shell metadata differ")
        for ordinal in range(functions):
            covered[first - 1 + ordinal] += 1
            if tuple(ao_map[:, first - 1 + ordinal]) != (shell_index, ordinal + 1):
                raise ValueError("AO shell/function addressing differs")
    if not np.all(covered == 1):
        raise ValueError("AO shells do not cover every function exactly once")
    return dict(
        kind=kind,
        radius_bohr=radius,
        mo_counts=counts,
        mo_mapping=mapping,
        continuum_flags=flags,
        coefficients=coefficients,
        stored_overlap=overlap,
        ao_mapping=ao_map,
        shell_descriptors=descriptors,
        geometry=geometry,
        shells=shells,
    )
