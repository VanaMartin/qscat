"""`config.resolved.yaml` must load back.

Every run writes the fully default-filled config beside its results, and the
artifact-store design leans on that file: expensive outputs live in object
storage precisely because the inputs needed to regenerate them stay in git. A
resolved config that cannot be handed back to `load_config` makes "re-run this
result" a manual reconstruction rather than a command, and nothing else in the
repository notices -- the file is written, looks right, and is never read.

The round trip is `load_config` -> `dataclasses.asdict` -> `yaml.safe_dump` ->
`load_config`, which is exactly what `write_artifacts` emits and what a person
re-running a published sweep would feed back in.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest
import yaml
from qscat_run.config import ExperimentConfig, load_config

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def _resolved_text(cfg: ExperimentConfig) -> str:
    """Byte-for-byte what `write_artifacts` writes as `config.resolved.yaml`."""
    return yaml.safe_dump(dataclasses.asdict(cfg), sort_keys=False)


def _comparable(cfg: ExperimentConfig) -> dict:
    d = dataclasses.asdict(cfg)
    # `config_dir` records where the file was read from, so it legitimately
    # differs between the original and its round-tripped copy.
    d.pop("config_dir", None)
    return d


@pytest.mark.parametrize("example", sorted(EXAMPLES.rglob("*.yaml")), ids=lambda p: p.stem)
def test_every_examples_resolved_config_loads_back_identically(
    example: Path, tmp_path: Path
) -> None:
    cfg = load_config(example)
    again = tmp_path / "config.resolved.yaml"
    again.write_text(_resolved_text(cfg))

    reloaded = load_config(again)

    assert _comparable(reloaded) == _comparable(cfg)


def test_a_resolved_config_survives_two_round_trips() -> None:
    """The second pass is not redundant: the first turns absent keys into
    explicit nulls, so a loader that tolerates `absent` but not `null` passes
    one round trip and fails the next."""
    import tempfile

    src = EXAMPLES / "n2-ve.yaml"
    cfg = load_config(src)
    with tempfile.TemporaryDirectory() as d:
        first = Path(d) / "first.yaml"
        first.write_text(_resolved_text(cfg))
        once = load_config(first)
        second = Path(d) / "second.yaml"
        second.write_text(_resolved_text(once))
        twice = load_config(second)
    assert _comparable(twice) == _comparable(once) == _comparable(cfg)
