"""Replay native RESON fits from saved K-matrices without repeating scattering.

Vary the background order and automatic detection threshold to expose fit
sensitivity. Every replay preserves its input, output and fitted parameters.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

from projects.ukrmol_co.analyze import native_resonances


def main() -> None:
    """Fit the B1 eigenphase sum with a specified native background/detection setup."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--nback", type=int, choices=range(1, 5), default=2)
    parser.add_argument("--dphz", type=float, choices=[0.7, 1.0, 1.3], default=1.0)
    args = parser.parse_args()
    workdir = args.workdir.resolve()
    config = json.loads((workdir / "config.json").read_text())
    if not (workdir / "result.json").exists():
        parser.error("A completed, analyzed calculation is required")
    directory = workdir / "refits" / f"nback{args.nback}-dphz{args.dphz}"
    directory.mkdir(parents=True, exist_ok=False)
    data = workdir / "output/CO/collected_scattering_data"
    for unit, source in (
        (10, data / "channels/channels.geom1.doublet.B1"),
        (21, data / "rmat_amplitudes/ramps.geom1.doublet.B1"),
        (921, data / "K-matrices/K-matrix.geom1.doublet.B1"),
    ):
        if not source.exists():
            raise FileNotFoundError(source)
        (directory / f"fort.{unit}").symlink_to(source)
    template = workdir / "output/CO/geom1/inputs/scattering.reson.doublet.B1.inp"
    source = template.read_text().replace(
        "&res ", f"&res\n  nback = {args.nback},\n  dphz = {args.dphz},\n "
    )
    (directory / "reson.inp").write_text(source)
    precision = config["precision"]
    env = os.environ | {
        "LD_LIBRARY_PATH": f"/opt/ukrmolp/lib.{precision}:" + os.environ.get("LD_LIBRARY_PATH", "")
    }
    with (directory / "reson.inp").open() as inp, (directory / "reson.out").open("w") as out:
        status = subprocess.run(
            [f"/opt/ukrmolp/bin.{precision}/reson"],
            cwd=directory,
            env=env,
            stdin=inp,
            stdout=out,
            stderr=subprocess.STDOUT,
        ).returncode
    result = {
        "nback": args.nback,
        "dphz": args.dphz,
        "exit_code": status,
        "native_resonances": native_resonances(directory / "reson.out") if status == 0 else [],
    }
    (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(status)


if __name__ == "__main__":
    main()
