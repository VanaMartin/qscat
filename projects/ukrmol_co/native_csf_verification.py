"""Check saved CSF correspondence and all CIDATA roots, never rewriting CIDATA."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

from projects.ukrmol_co.native_csf_map import read_target_csfs, signed_correspondence


def digest(data):
    return hashlib.sha256(data).hexdigest()


def verify_cidata(path, pairs):
    """Check integer phases, every root and an exact in-memory sign round trip."""
    sectors = []
    with FortranFile(path, "r") as stream:
        for irrep in ("A1", "B1", "B2", "A2"):
            for spin in ("singlet", "triplet"):
                label = f"{spin}.{irrep}"
                source, _, signs = pairs[label]
                if stream.read_record("u1").tobytes() != b"********                CIDATA  ":
                    raise ValueError("Wrong CIDATA block marker")
                header = stream.read_record("u1").tobytes()
                nset, nrec = np.frombuffer(header, "<i8", count=2)
                centers, dimension, roots, symmetry = np.frombuffer(
                    header, "<i8", count=4, offset=136
                )
                expected_roots = 13 if irrep == "A1" else 8
                if not (
                    len(header) == 208
                    and nset == len(sectors) + 1
                    and roots == expected_roots
                    and centers == 2
                    and nrec == centers + roots + 1
                    and dimension == len(signs)
                    and symmetry == source.symmetry
                    and header[16:136].decode().strip() == f"CO CI calculation - {spin} {irrep}"
                    and header[200:] == bytes(8)
                    and np.frombuffer(header, "<i8", count=1, offset=184)[0] == 14
                ):
                    raise ValueError("Unexpected CIDATA header/root inventory")
                np.testing.assert_allclose(
                    np.frombuffer(header, "<f8", count=2, offset=168),
                    [source.spin] * 2,
                    atol=0,
                    rtol=0,
                )
                core = np.frombuffer(header, "<f8", count=1, offset=192)[0]
                charges = []
                for _ in range(centers):
                    record = stream.read_record("u1").tobytes()
                    if (
                        len(record) != 40
                        or not np.isfinite(np.frombuffer(record, "<f8", offset=8)).all()
                    ):
                        raise ValueError("Invalid CIDATA physical centre")
                    charges.append(np.frombuffer(record, "<f8", count=1, offset=32)[0])
                np.testing.assert_allclose(charges, [6, 8], atol=0, rtol=0)
                record = stream.read_record("u1").tobytes()
                if len(record) != 8 * (2 * dimension + roots):
                    raise ValueError("Invalid mixed integer-phase/energy/diagonal record")
                phases = np.frombuffer(record, "<i8", count=dimension)
                if not np.array_equal(phases, source.phases):
                    raise ValueError("CIDATA phases differ from source CSF phases")
                energies = np.frombuffer(record, "<f8", count=roots, offset=8 * dimension)
                diagonals = np.frombuffer(record, "<f8", offset=8 * (dimension + roots))
                if not np.isfinite(np.r_[core, energies, diagonals]).all() or np.any(
                    np.diff(energies) < 0
                ):
                    raise ValueError("Invalid saved target energies/diagonals")
                vectors = []
                for root in range(roots):
                    record = stream.read_record("u1").tobytes()
                    if (
                        len(record) != 8 * (dimension + 1)
                        or np.frombuffer(record, "<i8", count=1)[0] != root + 1
                    ):
                        raise ValueError("Missing/reordered target root")
                    vectors.append(np.frombuffer(record, "<f8", offset=8))
                original = np.asarray(vectors)
                mapped = original * signs
                if (
                    not np.isfinite(original).all()
                    or (mapped * signs).tobytes() != original.tobytes()
                ):
                    raise ValueError("Target vector sign round trip failed")
                gram = original @ original.T
                mapped_gram = mapped @ mapped.T
                np.testing.assert_allclose(gram, np.eye(roots), atol=1e-9, rtol=0)
                np.testing.assert_allclose(mapped_gram, gram, atol=1e-14, rtol=0)
                sectors.append(
                    dict(
                        sector=label,
                        roots=int(roots),
                        dimension=int(dimension),
                        negative_csf_signs=int(np.count_nonzero(signs < 0)),
                        determinant_terms=len(source.coefficients),
                        sign_int64_sha256=digest(signs.astype("<i8").tobytes()),
                        original_vector_sha256=digest(original.tobytes()),
                        in_memory_mapped_vector_sha256=digest(mapped.tobytes()),
                        energy_bytes_sha256=digest(energies.tobytes()),
                        energies_hartree=(energies + core).tolist(),
                        vector_gram_error=float(np.max(abs(gram - np.eye(roots)))),
                        mapped_gram_change=float(np.max(abs(mapped_gram - gram))),
                        exact_coefficient_roundtrip=True,
                        source_centres=int(centers),
                        source_phase_match=True,
                    )
                )
        if stream._fp.tell() != path.stat().st_size:
            raise ValueError("Trailing CIDATA blocks")
    return sectors


def verify(root, contract):
    def check_inputs():
        for name, want in contract["input_sha256"].items():
            if digest((root / name).read_bytes()) != want:
                raise ValueError(f"Changed saved mapping input: {name}")
        for name, want in contract["native_source_sha256"].items():
            if digest((root / "ukrmol-in-3.3.0/source" / name).read_bytes()) != want:
                raise ValueError(f"Changed native source anchor: {name}")

    check_inputs()
    try:
        pairs = {}
        for label, unit in contract["sectors"].items():
            source = read_target_csfs(root / "target-only" / f"fort.{unit}")
            destination = read_target_csfs(root / "full-layout" / f"fort.{unit}")
            signs = signed_correspondence(source, destination)
            pairs[label] = (source, destination, signs)
        sectors = verify_cidata(root / "target-only/fort.26", pairs)
        return dict(
            signed_target_csf_correspondence_verdict=True,
            saved_cidata_phase_rule_verdict=True,
            native_consumer_execution_verdict=None,
            target_integral_equivalence_verdict=None,
            scattering_channel_map_verdict=None,
            retained_subspace_dominance_verdict=False,
            sectors=sectors,
            roots_verified=sum(s["roots"] for s in sectors),
            converted_native_files=0,
            native_calls=0,
            target_diagonalizations=0,
            property_exports=0,
            scattering_inner_solves=0,
            scope="Exact same-frame signed target-CSF correspondence and saved CIDATA phase rule. "
            "Retain original CIDATA unchanged; no live consumer, integral "
            "or scattering qualification.",
        )
    finally:
        check_inputs()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    report = {"signed_target_csf_correspondence_verdict": None}
    try:
        report.update(verify(args.root, json.loads(args.contract.read_text())))
    except BaseException as error:
        report["error"] = repr(error)
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - start
        report["maxrss_platform_units"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        with args.output.open("x") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")


if __name__ == "__main__":
    main()
