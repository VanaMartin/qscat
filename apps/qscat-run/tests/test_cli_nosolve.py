"""`click.testing.CliRunner` coverage of the no-solve commands:
`list`, `init`, `validate`. `run` is a stub for now (Task 2/3)."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
from qscat_run.cli import main


def test_list_shows_molecules() -> None:
    r = CliRunner().invoke(main, ["list"])
    assert r.exit_code == 0
    assert "F2" in r.output and "H2P" in r.output


def test_init_scaffolds_valid_config(tmp_path: Path) -> None:
    out = tmp_path / "f2.yaml"
    r = CliRunner().invoke(main, ["init", "F2", "--observables", "ve,da", "-o", str(out)])
    assert r.exit_code == 0, r.output
    assert out.exists()
    r2 = CliRunner().invoke(main, ["validate", str(out)])
    assert r2.exit_code == 0, r2.output


def test_init_h2p_scaffolds_valid_config(tmp_path: Path) -> None:
    out = tmp_path / "h2p.yaml"
    r = CliRunner().invoke(
        main, ["init", "H2P", "--observables", "dr", "--methods", "ti", "-o", str(out)]
    )
    assert r.exit_code == 0, r.output
    r2 = CliRunner().invoke(main, ["validate", str(out)])
    assert r2.exit_code == 0, r2.output


def test_validate_rejects_bad(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("molecule: XX\nmethods: [ti]\nobservables: []\noutput_dir: o\n")
    assert CliRunner().invoke(main, ["validate", str(bad)]).exit_code != 0


def test_run_dry_run_prints_plan_and_solves_nothing(tmp_path: Path) -> None:
    out = tmp_path / "f2.yaml"
    CliRunner().invoke(main, ["init", "F2", "--observables", "ve", "-o", str(out)])
    output_dir = tmp_path / "out"
    r = CliRunner().invoke(main, ["run", str(out), "--output", str(output_dir), "--dry-run"])
    assert r.exit_code == 0, r.output
    assert "molecule: F2" in r.output
    assert "grid[ti]" in r.output
    # --dry-run never solves or writes any artifact.
    assert not output_dir.exists()
