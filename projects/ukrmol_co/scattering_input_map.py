"""Read-only audit of saved C2v CONGEN phases and lowest-five property IDs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import mmap
import re
import resource
import struct
import time
from itertools import pairwise
from pathlib import Path

SECTORS = tuple(
    f"{spin}.{irrep}" for irrep in ("A1", "B1", "B2", "A2") for spin in ("singlet", "triplet")
)
DIMS = (15548, 24366, 15196, 24540, 15196, 24540, 15044, 24564)
UNITS = (710, 730, 711, 731, 712, 732, 713, 733)
INVENTORY = (58, 37, 37, 22)
TARGET = (7, 3, 3, 0)
WITNESS = struct.Struct("<15hbb")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def record_spans(data):
    """Offsets exclude the four-byte little-endian sequential markers."""
    spans, cursor = [], 0
    while cursor < len(data):
        if cursor + 8 > len(data):
            raise ValueError("Truncated record marker")
        (size,) = struct.unpack_from("<I", data, cursor)
        stop = cursor + 4 + size
        if stop + 4 > len(data) or data[cursor : cursor + 4] != data[stop : stop + 4]:
            raise ValueError("Truncated or mismatched native record")
        spans.append((cursor + 4, size))
        cursor = stop + 4
    return spans


def integers(data, start, count):
    return struct.unpack_from(f"<{count}q", data, start)


def assignment(text, name):
    clean = re.sub(r"!.*", "", text)
    if len(re.findall(rf"\b{name}\s*=", clean, re.I)) != 1:
        raise ValueError(f"Missing or ambiguous native assignment: {name}")
    matches = re.findall(rf"\b{name}\s*=\s*([^/=&]+)", clean, re.I)
    if len(matches) != 1:
        raise ValueError(f"Missing or ambiguous native assignment: {name}")
    return [int(x) for x in re.findall(r"[-+]?\d+", matches[0])]


def continuum_phase(occupied, coefficient_sign, continuum):
    """Fermionic sorting followed by moving the marked continuum electron last."""
    if len(set(occupied)) != len(occupied) or coefficient_sign not in (-1, 1):
        raise ValueError("Invalid determinant or coefficient sign")
    if occupied.count(continuum) != 1:
        raise ValueError("Missing unique marked continuum electron")
    parity = sum(a > b for i, a in enumerate(occupied) for b in occupied[i + 1 :])
    parity += sum(a > continuum for a in occupied)
    return coefficient_sign * (-1 if parity % 2 else 1)


def orbital_irrep(address):
    if address < 1 or address > 2 * sum(INVENTORY):
        raise ValueError("Spin-orbital address outside native inventory")
    start = 0
    for irrep, (full, target) in enumerate(zip(INVENTORY, TARGET, strict=True)):
        if address <= 2 * (start + full):
            ordinal = (address - 2 * start + 1) // 2
            return irrep, ordinal, ordinal <= target
        start += full
    raise ValueError("Unreachable orbital address")


def audit_congen(path, total_irrep, target_phases, witness):
    """Reconstruct every dophz-defining term; do not diagonalize or evaluate CSFs."""
    with path.open("rb") as stream, mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
        spans = record_spans(data)
        if len(spans) != 6 or spans[0][1] != 304:
            raise ValueError("Unsupported electron-only CONGEN records")
        header = integers(data, spans[0][0] + 120, 23)
        ncsf, nelt, terms, nsym, packed_size = (header[i] for i in (7, 8, 9, 11, 13))
        if (
            header[0],
            header[5],
            header[6],
            nelt,
            header[10],
            nsym,
            header[12],
            header[21],
            header[22],
        ) != (total_irrep, 154, 308, 15, 1, 4, 2, sum(DIMS), 8):
            raise ValueError("Unexpected scattering-CSF header")
        if struct.unpack_from("<2d", data, spans[0][0] + 128) != (0.5, 0.5):
            raise ValueError("Not doublet highest-weight scattering CSFs")
        if [size for _, size in spans[1:]] != [
            (sum(DIMS) + 32) * 8,
            (ncsf + nelt + 6 * nsym + 1) * 8,
            2 * (ncsf + 1) * 8,
            packed_size * 8,
            terms * 8,
        ]:
            raise ValueError("CONGEN record dimensions disagree")
        phase_record = integers(data, spans[1][0], sum(DIMS) + 32)
        phases, tail = phase_record[: sum(DIMS)], phase_record[sum(DIMS) :]
        dims, prototypes, mcont, gucont = [tail[i : i + 8] for i in range(0, 32, 8)]
        want_mcont = tuple(total_irrep ^ (i // 2) for i in range(8))
        if dims != DIMS or prototypes != (2,) * 8 or mcont != want_mcont or gucont != (0,) * 8:
            raise ValueError("Changed target/prototype/continuum grouping")
        metadata = integers(data, spans[2][0], spans[2][1] // 8)
        if metadata[:4] != INVENTORY or metadata[4 + 15 + ncsf :] != (
            0,
            *TARGET,
            *INVENTORY,
            0,
            0,
            0,
            0,
            *TARGET,
            0,
            0,
            0,
            0,
        ):
            raise ValueError("Changed orbital or target inventory")
        reference = metadata[4:19]
        counts = metadata[19 : 19 + ncsf]
        indices = integers(data, spans[3][0], 2 * (ncsf + 1))
        coefficients, packed = indices[: ncsf + 1], indices[ncsf + 1 :]
        if (
            coefficients[0] != 1
            or coefficients[-1] != terms + 1
            or packed[0] != 1
            or packed[-1] != packed_size + 1
        ):
            raise ValueError("Invalid determinant offset endpoints")
        if any(
            c <= 0 or b - a != c for (a, b), c in zip(pairwise(coefficients), counts, strict=True)
        ) or any(b <= a for a, b in pairwise(packed)):
            raise ValueError("Invalid determinant counts or offset ordering")
        reports, phase_cursor, csf_cursor = [], 0, 0
        for sector, dimension, mult, cont_irrep, canonical in zip(
            SECTORS, DIMS, prototypes, mcont, target_phases, strict=True
        ):
            sector_phases = phases[phase_cursor : phase_cursor + dimension]
            for j, saved in enumerate(sector_phases):
                csf = csf_cursor + j * mult
                pos = spans[4][0] + 8 * (packed[csf] - 1)
                (degree,) = integers(data, pos, 1)
                if degree < 0 or degree > nelt or packed[csf] + 1 + 2 * degree > packed[csf + 1]:
                    raise ValueError("Invalid first packed determinant")
                holes = integers(data, pos + 8, degree)
                particles = integers(data, pos + 8 * (degree + 1), degree)
                if len(set(holes)) != degree or any(h not in reference for h in holes):
                    raise ValueError("Invalid first-term holes")
                occupied = list(reference)
                for h, p in zip(holes, particles, strict=True):
                    occupied[reference.index(h)] = p
                continuum = [a for a in occupied if not orbital_irrep(a)[2]]
                if len(continuum) != 1 or sum(a % 2 for a in occupied) != 8:
                    raise ValueError("Wrong continuum occupation or total spin projection")
                irrep, ordinal, _ = orbital_irrep(continuum[0])
                if irrep != cont_irrep or ordinal != TARGET[irrep] + 1:
                    raise ValueError("First prototype has unexpected continuum orbital")
                (coefficient,) = struct.unpack_from(
                    "<d", data, spans[5][0] + 8 * (coefficients[csf] - 1)
                )
                if not math.isfinite(coefficient) or coefficient == 0:
                    raise ValueError("Invalid first determinant coefficient")
                sign = 1 if coefficient > 0 else -1
                if saved != continuum_phase(occupied, sign, continuum[0]):
                    raise ValueError("Saved scattering phase fails independent reconstruction")
                witness.write(WITNESS.pack(*occupied, sign, saved))
            reports.append(
                dict(
                    sector=sector,
                    dimension=dimension,
                    continuum_irrep=cont_irrep,
                    prototype_continua=mult,
                    phase_sha256=hashlib.sha256(
                        struct.pack(f"<{dimension}q", *sector_phases)
                    ).hexdigest(),
                    different_from_full_target_phases=sum(
                        a != b for a, b in zip(sector_phases, canonical, strict=True)
                    ),
                )
            )
            phase_cursor += dimension
            csf_cursor += dimension * mult
        if csf_cursor >= ncsf:
            raise ValueError("Missing L2 CSFs")
        return dict(
            total_csfs=ncsf,
            phase_defining_terms=phase_cursor,
            l2_csfs=ncsf - csf_cursor,
            sector_phases=reports,
            native_record_bytes=[size for _, size in spans],
        )


def property_map(text, reference=None, expected_moments=None):
    """Property IDs are energy ordered, not sector ordinals; moments stay untouched."""
    states, moments, headers = [], {}, []
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "6":
            headers.append((int(fields[4]), int(fields[5])))
        elif fields[0] == "5":
            states.append(
                (
                    int(fields[1]),
                    fields[-2].lower() + "." + fields[-1],
                    float(fields[8].replace("D", "E")),
                )
            )
        elif fields[0] == "1":
            key = tuple(int(x) for x in fields[1:8])
            value = float(fields[8].replace("D", "E"))
            if len(key) != 7 or key in moments or not math.isfinite(value):
                raise ValueError("Invalid or duplicate moment")
            moments[key] = value
    if (
        headers != [(40, len(moments))]
        or [r[0] for r in states] != list(range(1, 41))
        or (expected_moments is not None and len(moments) != expected_moments)
    ):
        raise ValueError("Incomplete property inventory")
    energies = [r[2] for r in states]
    if any(not math.isfinite(e) for e in energies) or any(
        b - a < -1e-10 for a, b in pairwise(energies)
    ):
        raise ValueError("Invalid property energy ordering")
    mapping, errors, labels = [], [], {}
    if set(r[1] for r in states) != set(SECTORS):
        raise ValueError("Wrong property sectors")
    for sector in SECTORS:
        rows = [r for r in states if r[1] == sector]
        if len(rows) != 5:
            raise ValueError("Missing five-per-sector properties")
        for root, (state_id, _, energy) in enumerate(rows, 1):
            mapping.append(state_id)
            labels[state_id] = dict(sector=sector, root=root, energy_hartree=energy)
            if reference is not None:
                errors.append(abs(energy - reference[sector][root - 1]))
    if errors and max(errors) > 1e-7:
        raise ValueError("Property energies differ from original CIDATA")
    covered, ranks = set(), set()
    for left, li, right, ri, kind, rank, projection in moments:
        if (
            left not in labels
            or right not in labels
            or kind != 3
            or rank not in (1, 2)
            or abs(projection) > rank
        ):
            raise ValueError("Invalid moment state/operator")
        if (
            li != SECTORS.index(labels[left]["sector"]) // 2
            or ri != SECTORS.index(labels[right]["sector"]) // 2
        ):
            raise ValueError("Moment irrep does not match property state")
        covered.update((left, right))
        ranks.add(rank)
    if covered != set(range(1, 41)) or ranks != {1, 2}:
        raise ValueError("Incomplete moment coverage")
    return dict(
        idtarg=mapping,
        states_by_id=labels,
        moment_count=len(moments),
        maximum_cidata_energy_error_hartree=max(errors) if errors else None,
    )


def cidata_reference(path):
    data = path.read_bytes()
    spans = record_spans(data)
    cursor = 0
    reference = {}
    for sector, dimension in zip(SECTORS, DIMS, strict=True):
        marker, size = spans[cursor]
        if data[marker : marker + size] != b"********                CIDATA  ":
            raise ValueError("Wrong CIDATA marker")
        header, _ = spans[cursor + 1]
        centres, dim, roots = integers(data, header + 136, 3)
        if centres != 2 or dim != dimension or roots != (13 if sector.endswith("A1") else 8):
            raise ValueError("Wrong original CIDATA root inventory")
        (nuclear,) = struct.unpack_from("<d", data, header + 192)
        mixed, _ = spans[cursor + 2 + centres]
        relative = struct.unpack_from(f"<{roots}d", data, mixed + 8 * dimension)
        reference[sector] = [v + nuclear for v in relative[:5]]
        cursor += 3 + centres + roots
    if cursor != len(spans):
        raise ValueError("Trailing CIDATA sectors")
    return reference


def check_cards(congen, scatci, swinterf, total_irrep, old_map):
    expected_cont = [total_irrep ^ (i // 2) for i in range(8)]
    expanded = [INVENTORY[i] - TARGET[i] for i in expected_cont]
    checks = [
        (congen, "nob", list(INVENTORY)),
        (congen, "nob0", list(TARGET)),
        (congen, "iscat", [2]),
        (congen, "nelect", [15]),
        (scatci, "ntgsym", [8]),
        (scatci, "mcont", expected_cont),
        (scatci, "notgt", expanded),
        (scatci, "numtgt", [5] * 8),
        (scatci, "ntgtf", [s for s in range(1, 9) for _ in range(5)]),
        (scatci, "ntgts", list(range(1, 6)) * 8),
        (swinterf, "idtarg", old_map),
        (swinterf, "ntarg", [40]),
        (swinterf, "last_coeff_saved", [1410]),
    ]
    for text, name, want in checks:
        if assignment(text, name) != want:
            raise ValueError(f"Saved card mismatch: {name}")
    return dict(expanded_continua=expanded, channel_coefficients=sum(5 * n for n in expanded))


def execute(root, output, contract):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        saved_scattering_map_verdict=None,
        native_scattering_destination_reader_verdict=None,
        scattering_solves=0,
        qc_calls=0,
        target_diagonalizations=0,
        property_exports=0,
        physical_same_manifold_verdict=False,
    )

    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Input hash changed: {name}")

    try:
        verify()
        phases = []
        for unit in UNITS:
            data = (root / contract["full_target_csfs"] / f"fort.{unit}").read_bytes()
            spans = record_spans(data)
            phases.append(integers(data, spans[1][0], spans[1][1] // 8))
        old_map = property_map(
            (root / contract["old_properties"]).read_text(), expected_moments=679
        )
        new_map = property_map(
            (root / contract["new_properties"]).read_text(),
            cidata_reference(root / contract["cidata"]),
            expected_moments=676,
        )
        report.update(
            old_properties=old_map,
            new_properties=new_map,
            changed_idtarg_positions=sum(
                a != b for a, b in zip(old_map["idtarg"], new_map["idtarg"], strict=True)
            ),
        )
        report["symmetries"] = {}
        for irrep, index in (("B1", 1), ("B2", 2)):
            files = contract["scattering"][irrep]
            with (output / f"{irrep}-first-terms.bin").open("xb") as witness:
                row = audit_congen(root / files["csfs"], index, phases, witness)
            row["cards"] = check_cards(
                *[(root / files[key]).read_text() for key in ("congen", "scatci", "swinterf")],
                index,
                old_map["idtarg"],
            )
            report["symmetries"][irrep] = row
        verify()
        report["saved_scattering_map_verdict"] = True
        report["new_idtarg_card"] = (
            "  idtarg = " + ",".join(str(x) for x in new_map["idtarg"]) + ","
        )
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except BaseException as error:
        report["error"] = repr(error)
        (output / "failure.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - started
        report["peak_rss_platform_units"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        (output / "execution.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    execute(args.root, args.output, json.loads(args.contract.read_text()))
