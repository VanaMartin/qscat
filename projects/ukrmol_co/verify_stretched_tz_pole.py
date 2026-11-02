"""Strict own-TZ provenance adapter for the shared saved-only pole witness reader."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from projects.ukrmol_co.verify_stretched_dz_pole import verify as verify_saved

ALLOCATION = dict(
    source_commit="69609bf1dc137608dc0d4960cc09fdb0b4c52f6a",
    source_archive_sha256="3d51fff8c3da7ee337f69c90a8977e9fd6d101c47f70331e264e1e547dfc79db",
    execution_path="projects/ukrmol_co/stretched-tz-pole-execution-20261102.json",
    source_count=52,
    runtime_controls=31,
    case="tz",
)


def verify(archive, collection):
    return verify_saved(archive, collection, allocation=ALLOCATION)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("collection", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            verify(args.archive, json.loads(args.collection.read_text())), indent=2, allow_nan=False
        )
    )
