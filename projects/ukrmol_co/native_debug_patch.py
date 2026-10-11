"""Add flushed pair/rank markers to a fresh checksum-pinned UKRmol build tree."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

CONTRACTED = "source/mpi-ci-diag/Modules/Hamiltonian/Contracted_Hamiltonian_module.f90"
ORIGINAL_SHA256 = "5dc113786e7a1e3cebe2dd6cf0cb453aa6f56fd83287b0f78b43b6df1186be34"


def instrument(text: str) -> str:
    """Bracket each class5/6/7 target pair; retain every original computational line."""
    for kind in (5, 6, 7):
        first = f"call this % evaluate_class_{kind}"
        lines = text.splitlines(keepends=True)
        indices = [i for i, line in enumerate(lines) if first in line]
        if len(indices) != 1:
            raise ValueError(f"No unique class{kind} call")
        start = indices[0]
        end = start
        while lines[end].rstrip().endswith("&"):
            end += 1

        def marker(status, kind=kind):
            return (
                "                        write(stdout, '(a,i0,a,i0,a,i0,a,i0,a,i0)') &\n"
                f"                            'CO_DIAG {status} class={kind} rank=', "
                "grid % grank, &\n"
                "                            ' pair=', i1, ',', i2, &\n"
                "                            ' targets=', "
                "this % options % num_target_state_sym(i1), &\n"
                "                            ',', this % options % num_target_state_sym(i2)\n"
                "                        flush(stdout)\n"
            )

        lines.insert(end + 1, marker("END"))
        lines.insert(start, marker("BEGIN"))
        text = "".join(lines)
    return text


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    path = args.source / CONTRACTED
    original = path.read_bytes()
    if hashlib.sha256(original).hexdigest() != ORIGINAL_SHA256:
        raise ValueError("Unrecognized or already patched native source")
    patched = instrument(original.decode()).encode()
    path.write_bytes(patched)
    with args.output.open("x") as stream:
        stream.write(
            json.dumps(
                {
                    "path": CONTRACTED,
                    "original_sha256": ORIGINAL_SHA256,
                    "instrumented_sha256": hashlib.sha256(patched).hexdigest(),
                    "scope": "Six flushed class5/6/7 pair/rank markers; no algorithm changes",
                },
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
