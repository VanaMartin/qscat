"""Bounded hourly immutable status captures for an existing native owner."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import subprocess
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--owner", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument(
        "--models", nargs="*", default=[], help="Nested model directories to observe"
    )
    parser.add_argument("--interval", type=float, default=3600)
    parser.add_argument("--duration", type=float, default=259200)
    args = parser.parse_args()
    if not (0 < args.interval <= args.duration <= 259200):
        parser.error("Positive bounded interval and duration required")
    args.output.mkdir()
    started = time.time()
    record = {
        "supervisor_pid": os.getpid(),
        "started_unix": started,
        "captures": [],
        "scope": "Observe only; no follow-on calculations",
    }

    def save():
        temporary = args.output / "execution.json.tmp"
        temporary.write_text(json.dumps(record, indent=2) + "\n")
        temporary.replace(args.output / "execution.json")

    save()
    try:
        while time.time() - started < args.duration:
            now = datetime.datetime.now(datetime.UTC)
            owner = json.loads((args.owner / "execution.json").read_text())
            capture = {"observed_utc": now.isoformat(), "owner": owner}
            for label, cmd in {
                "docker_ps": ["docker", "ps", "--format", "{{.Names}} {{.Status}}"],
                "docker_stats": [
                    "docker",
                    "stats",
                    "--no-stream",
                    "--format",
                    "{{.Name}} CPU={{.CPUPerc}} MEM={{.MemUsage}}",
                ],
                "memory": ["free", "-b"],
                "disk": ["df", "-B1", str(args.owner)],
            }.items():
                result = subprocess.run(cmd, text=True, capture_output=True, timeout=30)
                capture[label] = {
                    "exit": result.returncode,
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                }
            run = args.owner / "native" / args.run
            for sector in ("B1", "B2"):
                path = run / sector / "scatci.out"
                if path.exists():
                    with path.open("rb") as stream:
                        stream.seek(max(0, path.stat().st_size - 6000))
                        capture[sector + "_scatci_tail"] = stream.read().decode(errors="replace")
            for name in (
                "pilot.json",
                "resources.json",
                "independent-comparison.json",
                "independent-root-decision.json",
            ):
                path = run / name
                if path.exists():
                    capture[name] = json.loads(path.read_text())
            for model in args.models:
                directory = run / model
                for sector in ("B1", "B2"):
                    path = directory / sector / "scatci.out"
                    if path.exists():
                        with path.open("rb") as stream:
                            stream.seek(max(0, path.stat().st_size - 6000))
                            capture[f"{model}_{sector}_scatci_tail"] = stream.read().decode(
                                errors="replace"
                            )
                for name in ("pilot.json", "resources.json"):
                    path = directory / name
                    if path.exists():
                        capture[f"{model}_{name}"] = json.loads(path.read_text())
            filename = now.strftime("%Y%m%dT%H%M%S.json")
            with (args.output / filename).open("x") as stream:
                stream.write(json.dumps(capture, indent=2) + "\n")
            record["captures"].append(filename)
            record["heartbeat_unix"] = time.time()
            save()
            if "finished_unix" in owner:
                record["finish_reason"] = "Finite owner ended"
                break
            time.sleep(min(args.interval, max(0, args.duration - (time.time() - started))))
        else:
            record["finish_reason"] = "Bounded observation deadline reached"
        record["original_exit"] = 0
    except BaseException as error:
        record.update(original_exit=1, error=repr(error))
        raise
    finally:
        record["finished_unix"] = time.time()
        save()


if __name__ == "__main__":
    main()
