"""Package completed CO calibration evidence as a deterministic, indexed archive.

Use a collector snapshot to select completed batches and their runs. Large
integrals and K-matrices stay on the calculation host; orbital checkpoints,
inputs, logs, raw phases, cross sections and source snapshots travel with the
archive. Extract outside a checkout to avoid duplicate Python test packages.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

SUFFIXES = {
    ".json",
    ".tsv",
    ".out",
    ".err",
    ".log",
    ".inp",
    ".molden",
    ".py",
    ".pm",
    ".pl",
    ".csv",
    ".chk",
    ".md",
    ".f90",
}
NAMES = {"target.energies", "log_file.0", "UKRmol-in-COPYING"}


def archive(
    root: Path,
    snapshot: Path,
    attachments: list[Path],
    output: Path,
    *,
    diagnostics: list[str] | None = None,
) -> dict:
    """Select a completed snapshot and write its bytes with normalized tar metadata."""
    report = json.loads(snapshot.read_text())
    paths = {}
    for kind, names in (
        ("batches", [batch["name"] for batch in report["batches"]]),
        ("runs", list(report["runs"])),
        ("diagnostics", diagnostics or []),
    ):
        for name in names:
            if Path(name).name != name or name in (".", ".."):
                raise ValueError(f"Invalid {kind} directory name: {name}")
            directory = root / kind / name
            if not directory.is_dir():
                if kind == "runs" and report["runs"][name]["status"] == "setup_failed":
                    continue
                raise FileNotFoundError(directory)
            if kind == "diagnostics" and not (directory / "execution.json").is_file():
                raise ValueError(f"Diagnostic has no completed execution record: {name}")
            for path in directory.rglob("*"):
                if (
                    path.is_file()
                    and not path.is_symlink()
                    and (
                        path.suffix in SUFFIXES
                        or path.name in NAMES
                        or path.name.startswith(("eigenph.", "xsec."))
                    )
                ):
                    paths[path.relative_to(root).as_posix()] = path
    for path in [snapshot, *attachments]:
        if path.name in paths or path.name == "file-index.json":
            raise ValueError(f"Duplicate archive entry: {path.name}")
        paths[path.name] = path
    files = {
        name: {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
        for name, path in sorted(paths.items())
    }
    index = {
        "format": "ukrmol-calibration-evidence-v1",
        "snapshot": snapshot.name,
        "files": files,
    }
    index_bytes = (json.dumps(index, indent=2) + "\n").encode()
    # Exclusive creation protects a previous publication from accidental replacement.
    with (
        output.open("xb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
    ):
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as bundle:
            for name in sorted([*paths, "file-index.json"]):
                data = index_bytes if name == "file-index.json" else paths[name].read_bytes()
                if name != "file-index.json" and (
                    len(data) != files[name]["bytes"]
                    or hashlib.sha256(data).hexdigest() != files[name]["sha256"]
                ):
                    raise ValueError(f"Evidence changed while packaging: {name}")
                info = tarfile.TarInfo(f"state-averaged-evidence/{name}")
                info.size = len(data)
                info.mode = 0o644
                info.uid = info.gid = info.mtime = 0
                bundle.addfile(info, io.BytesIO(data))
    return {
        "archive": str(output),
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "bytes": output.stat().st_size,
        "indexed_files": len(files),
        "completed_attempts": len(report["runs"]),
        "completed_batches": len(report["batches"]),
    }


def main() -> None:
    """Package a copied evidence root, its aggregate and explicitly named attachments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--attachments", type=Path, nargs="*", default=[])
    parser.add_argument("--diagnostics", nargs="*", default=[], help="Completed diagnostic names")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = archive(
        args.root, args.snapshot, args.attachments, args.output, diagnostics=args.diagnostics
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
