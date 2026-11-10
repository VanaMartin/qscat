"""Synthetic QC replay contracts; no PySCF calculations or campaign data."""

import hashlib
import json
import sys
from unittest.mock import patch

import pytest

from projects.ukrmol_co import run, target
from projects.ukrmol_co.reuse_qc import FILES, authenticated_qc, replay_qc


def bundle(tmp_path):
    source = tmp_path / "accepted"
    source.mkdir()
    config = dict(
        bond_length=2.2,
        basis="test-basis",
        model="CAS-A",
        orbitals="state-averaged",
        frozen_orbitals=2,
        active_orbitals=[3, 1, 1, 0],
        sa_singlet_roots=[1] * 4,
        sa_triplet_roots=[1] * 4,
    )
    states = [
        dict(
            spin=spin,
            irrep=irrep,
            root=1,
            energy_hartree=-10 + i * 0.1,
            spin_square=0 if spin == "singlet" else 2,
            weight=1 / 8,
        )
        for i, (spin, irrep) in enumerate(
            (s, ir) for s in ("singlet", "triplet") for ir in ("A1", "B1", "B2", "A2")
        )
    ]
    energies = [s["energy_hartree"] for s in states]
    diagnostics = dict(
        orbital_converged=True,
        ci_converged=[True] * 8,
        fresh_ci_converged=[True] * 8,
        lowest_root_max_energy_difference_hartree=0,
        state_energies_hartree=energies,
        fresh_state_energies_hartree=energies,
        fresh_spin_square=[s["spin_square"] for s in states],
    )
    cz = -2.2 * 15.9994 / (12.0110 + 15.9994)
    molden = b"[Molden Format]\n[Atoms] AU\n[GTO]\n[MO]\n"
    orbitals = [
        dict(irrep=ir, space="core" if i < 2 else "active", energy_hartree=-1)
        for i, ir in enumerate(["A1"] * 5 + ["B1", "B2"])
    ]
    report = dict(
        converged=True,
        optimization=diagnostics,
        method="equal-weight multi-spin/multi-irrep SA-CASSCF",
        mo_orthogonality_max_error=0,
        ground_energy_hartree=-10,
        states=states,
        atoms_bohr=[["C", 0, 0, cz], ["O", 0, 0, cz + 2.2]],
        orbitals=orbitals,
        molden_sha256=hashlib.sha256(molden).hexdigest(),
    )
    for name, data in [
        ("config.json", config),
        ("target.json", report),
        ("target-diagnostics.json", diagnostics),
        ("resources.json", dict(exit_code=0)),
    ]:
        (source / name).write_text(json.dumps(data))
    (source / "co.molden").write_bytes(molden)
    (source / "co.casscf.chk").write_bytes(b"\x89HDF\r\n\x1a\nsynthetic-not-a-wavefunction")
    manifest = tmp_path / "manifest.json"

    def freeze():
        spec = dict(
            source_directory="accepted",
            files={n: hashlib.sha256((source / n).read_bytes()).hexdigest() for n in FILES},
        )
        manifest.write_text(json.dumps(spec))
        return hashlib.sha256(manifest.read_bytes()).hexdigest()

    sha = freeze()
    destination = tmp_path / "native"
    destination.mkdir()
    config.update(
        workdir=str(destination),
        reuse_qc_manifest=str(manifest),
        reuse_qc_sha256=sha,
        target_singlet_roots=[1, 1, 1, 3],
        target_triplet_roots=[1, 1, 1, 3],
    )
    return source, manifest, config, freeze


def test_replay_preserves_original_bytes_and_allows_more_native_roots(tmp_path):
    source, _, config, _ = bundle(tmp_path)
    before = {name: (source / name).read_bytes() for name in FILES}
    molden = tmp_path / "geometry.molden"
    result = replay_qc(config, molden)
    assert result == json.loads(before["target.json"])
    assert molden.read_bytes() == before["co.molden"]
    for name in ("target.json", "target-diagnostics.json", "co.casscf.chk"):
        assert (tmp_path / "native" / name).read_bytes() == before[name]
    assert {name: (source / name).read_bytes() for name in FILES} == before
    assert json.loads((tmp_path / "native/qc-reuse.json").read_text())["new_qc_calls"] == 0
    with pytest.raises(ValueError, match="overwrite"):
        replay_qc(config, molden)


@pytest.mark.parametrize("name", FILES)
def test_any_changed_input_rejects_without_outputs(tmp_path, name):
    source, _, config, _ = bundle(tmp_path)
    with (source / name).open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="digest"):
        replay_qc(config, tmp_path / "new.molden")
    assert not list((tmp_path / "native").iterdir())


@pytest.mark.parametrize(
    "key,value",
    [
        ("bond_length", 2.3),
        ("basis", "wrong"),
        ("frozen_orbitals", 4),
        ("active_orbitals", [4, 1, 1, 0]),
        ("sa_singlet_roots", [2, 1, 1, 1]),
    ],
)
def test_protocol_change_rejects(tmp_path, key, value):
    _, manifest, config, _ = bundle(tmp_path)
    config[key] = value
    with pytest.raises(ValueError, match="protocol"):
        authenticated_qc(manifest, config["reuse_qc_sha256"], config)


@pytest.mark.parametrize(
    "defect",
    ["exit", "convergence", "fresh_ci", "spin", "geometry", "checkpoint", "molden", "weight"],
)
def test_rehashed_but_invalid_bundle_rejects(tmp_path, defect):
    source, _, config, freeze = bundle(tmp_path)
    report = json.loads((source / "target.json").read_text())
    if defect == "exit":
        (source / "resources.json").write_text('{"exit_code":1}')
    elif defect == "checkpoint":
        (source / "co.casscf.chk").write_bytes(b"not-HDF5")
    elif defect == "molden":
        (source / "co.molden").write_bytes(b"truncated")
    else:
        if defect == "convergence":
            report["converged"] = False
        if defect == "fresh_ci":
            report["optimization"]["fresh_ci_converged"][0] = False
            (source / "target-diagnostics.json").write_text(json.dumps(report["optimization"]))
        if defect == "spin":
            report["states"][0]["spin_square"] = float("nan")
        if defect == "geometry":
            report["atoms_bohr"][0][3] += 0.01
        if defect == "weight":
            report["states"][0]["weight"] = 0.5
        (source / "target.json").write_text(json.dumps(report))
    config["reuse_qc_sha256"] = freeze()
    with pytest.raises(ValueError):
        replay_qc(config, tmp_path / "new.molden")
    assert not list((tmp_path / "native").iterdir())


def test_target_entrypoint_never_calls_builder_when_replaying(tmp_path, monkeypatch):
    source, _, config, _ = bundle(tmp_path)
    config_path = tmp_path / "native/config.json"
    config_path.write_text(json.dumps(config))
    geom = tmp_path / "geom"
    geom.mkdir()
    monkeypatch.chdir(geom)
    monkeypatch.setattr(
        sys, "argv", ["target", "--config", str(config_path), "--output", str(geom / "target.out")]
    )
    with patch.object(target, "build_target", side_effect=AssertionError("QC rerun")) as build:
        target.main()
    build.assert_not_called()
    assert (geom / "co.molden").read_bytes() == (source / "co.molden").read_bytes()


def test_cli_rejects_unpinned_replay_before_scripts_or_workdir(tmp_path, monkeypatch):
    monkeypatch.setattr(
        sys, "argv", ["run", "--workdir", str(tmp_path / "new"), "--reuse-qc-manifest", "missing"]
    )
    with patch.object(run, "acquire_scripts") as acquire, pytest.raises(SystemExit):
        run.main()
    acquire.assert_not_called()
    assert not (tmp_path / "new").exists()
