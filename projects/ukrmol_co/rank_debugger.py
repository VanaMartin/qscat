"""Launch one native MPI rank under GDB without treating debugger exits as native exits."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("program", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.program:
        parser.error("Native command required")
    rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", "0"))
    directory = Path(os.environ["CO_DEBUG_DIRECTORY"])
    command = [
        "gdb",
        "--batch",
        "--return-child-result",
        "--quiet",
        "-x",
        str(Path(__file__).with_name("rank_gdb.py")),
        "--args",
        *args.program,
    ]
    with (directory / f"rank-{rank}.gdb.log").open("x") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=False)
    path = directory / f"rank-{rank}.json"
    record = (
        json.loads(path.read_text())
        if path.exists()
        else {"rank": rank, "missing_inferior_event": True}
    )
    record.update(debugger_original_exit=result.returncode, command=command)
    path.write_text(json.dumps(record, indent=2) + "\n")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
