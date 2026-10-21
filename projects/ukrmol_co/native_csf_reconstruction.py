"""Separate scalar determinant-expansion reconstruction of saved native CSF pairs.

Uses struct record framing and canonical occupation keys, not the mapping reader
or its packed-array correspondence implementation. No native computation occurs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import time
from pathlib import Path


def records(path):
    """Read all sequential records without SciPy's Fortran reader."""
    result = []
    with path.open("rb") as stream:
        while marker := stream.read(4):
            if len(marker) != 4:
                raise ValueError("Truncated record marker")
            (length,) = struct.unpack("<I", marker)
            if length > path.stat().st_size - stream.tell() - 4:
                raise ValueError("Record extends beyond EOF")
            data = stream.read(length)
            if stream.read(4) != marker:
                raise ValueError("Mismatched record trailer")
            result.append(data)
    if len(result) != 6:
        raise ValueError("Wrong CSF record count")
    return result


def expansion(path):
    """Yield canonical determinant keys and coefficients for each physical CSF."""
    rows = records(path)
    h = struct.unpack("<qdddd9q6qdqq", rows[0][120:])
    dimension, electrons, terms, nsym, packed_size = h[7], h[8], h[9], h[11], h[13]

    def ints(data):
        return struct.unpack("<" + "q" * (len(data) // 8), data)

    meta = ints(rows[2])
    inventory = meta[:nsym]
    reference = meta[nsym : nsym + electrons]
    counts = meta[nsym + electrons : nsym + electrons + dimension]
    offsets = ints(rows[3])
    coefficient_offsets = offsets[: dimension + 1]
    packed_offsets = offsets[dimension + 1 :]
    substitutions = ints(rows[4])
    coefficients = struct.unpack("<" + "d" * terms, rows[5])
    if len(substitutions) != packed_size:
        raise ValueError("Packed length mismatch")

    # Decode native address directly into the ordered tuple (irrep, ordinal, spin).
    # This route does not construct or use the producer's target-address lookup.
    def orbital(address):
        spatial, spin = divmod(address - 1, 2)
        for irrep, count in enumerate(inventory):
            if spatial < count:
                if spatial >= (7, 3, 3, 0)[irrep]:
                    raise ValueError("Occupied continuum orbital")
                return irrep, spatial + 1, spin
            spatial -= count
        raise ValueError("Invalid spin orbital")

    for j, count in enumerate(counts):
        cursor = packed_offsets[j] - 1
        canonical = {}
        for term in range(count):
            degree = substitutions[cursor]
            holes = substitutions[cursor + 1 : cursor + degree + 1]
            particles = substitutions[cursor + degree + 1 : cursor + 2 * degree + 1]
            occupied = list(reference)
            for hole, particle in zip(holes, particles, strict=True):
                occupied[reference.index(hole)] = particle
            slots = [orbital(a) for a in occupied]
            inversions = sum(
                slots[a] > slots[b] for a in range(electrons) for b in range(a + 1, electrons)
            )
            phase = (-1) ** inversions
            key = tuple(sorted(slots))
            if len(set(key)) != electrons or key in canonical:
                raise ValueError("Repeated electron or duplicate canonical determinant")
            canonical[key] = phase * coefficients[coefficient_offsets[j] - 1 + term]
            cursor += 2 * degree + 1
        if cursor != packed_offsets[j + 1] - 1:
            raise ValueError("Unconsumed packed CSF terms")
        yield canonical


def reconstruct(root, contract, original):
    expected = {s["sector"]: s for s in original["sectors"]}
    outcomes = []
    for label, unit in contract["sectors"].items():
        count = negative = terms = 0
        signs = []
        streams = [
            expansion(root / directory / f"fort.{unit}")
            for directory in ("target-only", "full-layout")
        ]
        for source, destination in zip(*streams, strict=True):
            if source.keys() != destination.keys():
                raise ValueError("Different canonical determinant support")
            first = next(iter(source))
            sign = 1 if source[first] * destination[first] > 0 else -1
            for key, value in source.items():
                if struct.pack("<d", sign * value) != struct.pack("<d", destination[key]):
                    raise ValueError(
                        "Canonical determinant expansion violates signed correspondence"
                    )
            if abs(sum(v * v for v in source.values()) - 1) > 1e-12:
                raise ValueError("Nonunit canonical CSF norm")
            count += 1
            terms += len(source)
            negative += sign < 0
            signs.append(sign)
        sign_hash = hashlib.sha256(struct.pack("<" + "q" * len(signs), *signs)).hexdigest()
        want = expected[label]
        if (count, terms, negative, sign_hash) != (
            want["dimension"],
            want["determinant_terms"],
            want["negative_csf_signs"],
            want["sign_int64_sha256"],
        ):
            raise ValueError("Independent sign map differs from producer")
        outcomes.append(
            dict(
                sector=label,
                csfs=count,
                canonical_terms=terms,
                negative_signs=negative,
                sign_int64_sha256=sign_hash,
            )
        )
    return dict(
        canonical_expansion_reconstruction_verdict=True,
        sectors=outcomes,
        native_calls=0,
        converted_native_files=0,
        independence="Separate struct framing, scalar native-address to (irrep,ordinal,spin) "
        "decoding, full canonical determinant keys/parities and coefficient dictionaries. "
        "No producer reader or phase/sign map used to derive the answer.",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    result = {"canonical_expansion_reconstruction_verdict": None}
    contract = json.loads(args.contract.read_text())

    def verify_inputs():
        for name, want in contract["input_sha256"].items():
            if hashlib.sha256((args.root / name).read_bytes()).hexdigest() != want:
                raise ValueError(f"Changed saved input: {name}")

    try:
        verify_inputs()
        result.update(reconstruct(args.root, contract, json.loads(args.original.read_text())))
        verify_inputs()
    except BaseException as error:
        result["error"] = repr(error)
        raise
    finally:
        result["wall_seconds"] = time.monotonic() - start
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")


if __name__ == "__main__":
    main()
