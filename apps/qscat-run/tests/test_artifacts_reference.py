"""Regression coverage for the `reference:` artifact-writing wiring.

`test_reference.py` covers the parsing/loading layer (`load_config` /
`validate_config` / `load_reference`) but never drives `write_artifacts`
itself, so a defect in the resolve-and-merge glue -- e.g. reference values
leaking into `cross_section.csv`, or a relative path resolving against the
wrong base directory -- would pass every test committed with the feature.
This file closes that gap: it drives `write_artifacts` directly with a
minimal, hand-built `ExperimentResult`/`ExperimentConfig` (no real solve).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure
from qscat_run.artifacts import _disambiguated_labels, _write_cross_section_png, write_artifacts
from qscat_run.config import ArtifactSpec, ExperimentConfig, Observable, ReferenceSpec
from qscat_run.runner import ExperimentResult


def _fake_result(cfg: ExperimentConfig) -> ExperimentResult:
    return ExperimentResult(
        energies=np.array([0.03, 0.05]),
        cross_sections={"ti:ve:v0->0": np.array([1.0, 2.0])},
        wavefunctions=[],
        resolved_cfg=cfg,
    )


def _cfg(
    tmp_path: Path, *, ref_path: str, channels: tuple[int, ...] | None, label: str | None
) -> ExperimentConfig:
    return ExperimentConfig(
        molecule="N2",
        methods=("ti",),
        observables=(Observable(kind="ve", channels=1),),
        output_dir="unused",
        artifacts=ArtifactSpec(cross_section=True),
        reference=(ReferenceSpec(path=ref_path, format="houfek", label=label, channels=channels),),
        config_dir=str(tmp_path),
    )


def test_write_artifacts_keeps_the_references_own_energy_axis(tmp_path: Path) -> None:
    # The reference's OWN energies (0.1/0.2/0.3) are deliberately disjoint
    # from the run's (0.03/0.05) so any accidental merge/interpolation is
    # unmistakable in the assertions below.
    (tmp_path / "ref.dat").write_text("0.1 10.0\n0.2 20.0\n0.3 30.0\n")
    cfg = _cfg(tmp_path, ref_path="ref.dat", channels=(0,), label=None)
    result = _fake_result(cfg)
    out_dir = tmp_path / "out"

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    assert (out_dir / "reference.csv").exists()
    assert (out_dir / "reference.npz").exists()

    with (out_dir / "reference.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert [r["series"] for r in rows] == ["ref:ve:ch0"] * 3
    assert [float(r["energy"]) for r in rows] == [0.1, 0.2, 0.3]
    assert [float(r["sigma"]) for r in rows] == [10.0, 20.0, 30.0]

    npz = np.load(out_dir / "reference.npz")
    np.testing.assert_allclose(npz["ref:ve:ch0:energy"], [0.1, 0.2, 0.3])
    np.testing.assert_allclose(npz["ref:ve:ch0:sigma"], [10.0, 20.0, 30.0])


def test_the_run_sweep_has_no_reference_arrays(tmp_path: Path) -> None:
    """The own-axis rule, pinned: a reference keeps its own energy axis and
    must never be interleaved into the run's own sweep, whose rows are the
    run's energies. (This guarded `cross_section.csv` until sweeps moved to
    the machine tier and stopped being written as text at all.)"""
    (tmp_path / "ref.dat").write_text("0.1 10.0\n0.2 20.0\n")
    cfg = _cfg(tmp_path, ref_path="ref.dat", channels=(0,), label=None)
    result = _fake_result(cfg)
    out_dir = tmp_path / "out"

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    keys = list(np.load(out_dir / "cross_section.npz").files)
    assert keys == ["energy", "ti:ve:v0->0"]
    assert not any(k.startswith("ref:") for k in keys)


def test_relative_reference_path_resolves_against_config_dir_not_cwd(tmp_path: Path) -> None:
    """`ref_base` in `write_artifacts` must come from `cfg.config_dir`, not
    the process's current working directory -- this test's cwd (the repo
    root, per pytest's default) has no `sub/ref.dat`, so a wrong base would
    raise `FileNotFoundError` rather than silently succeed."""
    (tmp_path / "sub").mkdir()
    # `np.loadtxt` needs >=2 rows to infer a 2-D table -- a single-row file
    # collapses to 1-D and would raise for an unrelated reason.
    (tmp_path / "sub" / "ref.dat").write_text("0.5 7.0\n0.6 8.0\n")
    cfg = _cfg(tmp_path, ref_path="sub/ref.dat", channels=(0,), label=None)
    result = _fake_result(cfg)
    out_dir = tmp_path / "out"

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")  # must not raise

    npz = np.load(out_dir / "reference.npz")
    np.testing.assert_allclose(npz["ref:ve:ch0:energy"], [0.5, 0.6])


def test_no_reference_block_writes_no_reference_files(tmp_path: Path) -> None:
    cfg = ExperimentConfig(
        molecule="N2",
        methods=("ti",),
        observables=(Observable(kind="ve", channels=1),),
        output_dir="unused",
        artifacts=ArtifactSpec(cross_section=True),
        config_dir=str(tmp_path),
    )
    result = _fake_result(cfg)
    out_dir = tmp_path / "out"

    write_artifacts(result, cfg, out_dir, timestamp="2026-01-01T00:00:00")

    assert not (out_dir / "reference.csv").exists()
    assert not (out_dir / "reference.npz").exists()


def test_write_artifacts_wires_the_reference_label_into_the_png_legend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end through `write_artifacts` itself (not just the
    `_disambiguated_labels` helper in isolation): its own loop over
    `cfg.reference` must build and pass the label map through to the PNG
    writer."""
    closed: list[Figure] = []
    real_close = plt.close
    monkeypatch.setattr(plt, "close", lambda fig: closed.append(fig))

    (tmp_path / "ref.dat").write_text("0.1 10.0 1.0\n0.2 20.0 2.0\n")
    cfg = _cfg(tmp_path, ref_path="ref.dat", channels=(0, 1), label="Houfek (2006)")
    result = _fake_result(cfg)

    write_artifacts(result, cfg, tmp_path / "out", timestamp="2026-01-01T00:00:00")

    assert len(closed) == 1
    _, legend_labels = closed[0].axes[0].get_legend_handles_labels()
    assert "Houfek (2006) (ch0)" in legend_labels
    assert "Houfek (2006) (ch1)" in legend_labels
    for fig in closed:
        real_close(fig)


# --- legend labels -----------------------------------------------------------


def test_disambiguated_labels_none_falls_back_to_the_series_key() -> None:
    assert _disambiguated_labels(None, ["ref:ve:ch0", "ref:ve:ch1"]) == {}


def test_disambiguated_labels_single_channel_uses_label_verbatim() -> None:
    assert _disambiguated_labels("Houfek (2006)", ["ref:ve:ch0"]) == {"ref:ve:ch0": "Houfek (2006)"}


def test_disambiguated_labels_multi_channel_appends_the_channel_suffix() -> None:
    labels = _disambiguated_labels("Houfek (2006)", ["ref:ve:ch0", "ref:ve:ch1", "ref:ve:ch2"])
    assert labels == {
        "ref:ve:ch0": "Houfek (2006) (ch0)",
        "ref:ve:ch1": "Houfek (2006) (ch1)",
        "ref:ve:ch2": "Houfek (2006) (ch2)",
    }
    # unambiguous: no two channels share a label
    assert len(set(labels.values())) == 3


def test_cross_section_png_legend_uses_the_disambiguated_reference_labels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Drives `_write_cross_section_png` directly (bypassing `write_artifacts`
    so the figure can be inspected before it's closed) and reads the actual
    legend text back off the axes -- `label: "Houfek (2006)"` must show up as
    a real, disambiguated legend entry, not just sit unused in the config."""
    closed: list[Figure] = []
    real_close = plt.close
    monkeypatch.setattr(plt, "close", lambda fig: closed.append(fig))

    series = {"ti:ve:v0->0": np.array([1.0, 2.0])}
    reference = {
        "ref:ve:ch0": (np.array([0.1, 0.2]), np.array([10.0, 20.0])),
        "ref:ve:ch1": (np.array([0.1, 0.2]), np.array([1.0, 2.0])),
    }
    labels = _disambiguated_labels("Houfek (2006)", list(reference))

    _write_cross_section_png(tmp_path / "cs.png", np.array([0.03, 0.05]), series, reference, labels)

    assert len(closed) == 1
    fig = closed[0]
    _, legend_labels = fig.axes[0].get_legend_handles_labels()
    assert "Houfek (2006) (ch0)" in legend_labels
    assert "Houfek (2006) (ch1)" in legend_labels
    assert "ti:ve:v0->0" in legend_labels
    # not the raw, undisambiguated series keys
    assert "ref:ve:ch0" not in legend_labels
    assert "ref:ve:ch1" not in legend_labels
    for f in closed:
        real_close(f)
