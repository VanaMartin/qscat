"""The two-geometry owner must fail closed before any automatic continuation."""

import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from projects.ukrmol_co import neutral_core_owner as owner


def capsule(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    files = {}
    for name in owner.REQUIRED_SOURCE:
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        data = (
            Path(owner.__file__).read_bytes() if name == owner.OWNER_PATH else b"generated source\n"
        )
        path.write_bytes(data)
        files[name] = hashlib.sha256(data).hexdigest()
    manifest = {
        "git_sha": "a" * 40,
        "source_files": files,
        "image": owner.IMAGE,
        "received_unix": time.time() - 5,
        "deadline_unix": time.time() + 600,
    }
    path = tmp_path / "source-manifest.json"
    path.write_text(json.dumps(manifest))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_container_command_is_fixed_and_bounded(tmp_path):
    command = owner.container_command(tmp_path, "owned-r2450", 2.45, 123)
    for flag in [
        "--pull=never",
        "--network=none",
        "--memory=24g",
        "--memory-swap=24g",
        "--cpuset-cpus=8-11",
        "--cpus=4",
        "--read-only",
        "--workdir=/work/source",
    ]:
        assert flag in command
    assert owner.IMAGE in command
    assert "123s" in command
    assert command[-1] == "20000"
    assert "projects.ukrmol_co.neutral_core_diagnostic" in command
    assert "--execute-reviewed" in command
    assert not any(x in command for x in ["pull", "build", "pip"])


@pytest.mark.parametrize("radius", [2.55, True, float("nan")])
def test_owner_rejects_unadmitted_geometry(tmp_path, radius):
    with pytest.raises(ValueError):
        owner.container_command(tmp_path, "owned", radius, 120)


def test_changed_frozen_source_stops_before_host_or_container(tmp_path, monkeypatch):
    pin = capsule(tmp_path)
    (tmp_path / "source/projects/ukrmol_co/neutral_core_diagnostic.py").write_text("changed")
    monkeypatch.setattr(owner, "host_ready", lambda root: pytest.fail("host should not be used"))
    assert owner.own(tmp_path, pin) != 0
    result = json.loads((tmp_path / "execution.json").read_text())
    assert result["pairs"] == [] and result["status"] == "stopped"


@pytest.mark.parametrize("failure_index", [0, 1, None])
def test_sequence_stops_on_first_pair_failure(tmp_path, monkeypatch, failure_index):
    pin = capsule(tmp_path)
    monkeypatch.setattr(owner, "host_ready", lambda root: {"generated_host": True})
    calls = []

    def pair(root, radius, name, budget, deadline):
        i = len(calls)
        calls.append(radius)
        return {
            "radius_bohr": radius,
            "original_command_exit": 1 if i == failure_index else 0,
            "container_exit": 1 if i == failure_index else 0,
            "OOMKilled": False,
            "container_running": False,
            "diagnostic_complete": i != failure_index,
        }

    monkeypatch.setattr(owner, "run_pair", pair)
    monkeypatch.setattr(owner, "owned_running", lambda names: "")
    result = owner.own(tmp_path, pin)
    assert (result == 0) == (failure_index is None)
    assert calls == ([2.45] if failure_index == 0 else [2.45, 2.5])


def test_expired_allocation_never_launches(tmp_path, monkeypatch):
    pin = capsule(tmp_path)
    path = tmp_path / "source-manifest.json"
    value = json.loads(path.read_text())
    value["deadline_unix"] = time.time() - 1
    path.write_text(json.dumps(value))
    pin = hashlib.sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(owner, "host_ready", lambda root: pytest.fail("expired"))
    assert owner.own(tmp_path, pin) != 0


def test_manifest_digest_is_checked(tmp_path):
    capsule(tmp_path)
    assert owner.own(tmp_path, "0" * 64) != 0


def test_missing_source_is_a_retained_stop(tmp_path, monkeypatch):
    pin = capsule(tmp_path)
    (tmp_path / "source/projects/ukrmol_co/run.py").unlink()
    monkeypatch.setattr(owner, "host_ready", lambda root: pytest.fail("missing source"))
    assert owner.own(tmp_path, pin) != 0
    assert json.loads((tmp_path / "execution.json").read_text())["source_integrity_after"] is False


def test_deadline_kills_only_cid_created_in_owned_directory(tmp_path, monkeypatch):
    calls = []
    state = {
        "Image": owner.IMAGE,
        "State": {"ExitCode": 137, "OOMKilled": False, "Running": False},
        "HostConfig": {"Memory": 24 * 2**30, "MemorySwap": 24 * 2**30, "CpusetCpus": "8-11"},
    }
    monkeypatch.setattr(owner, "_output", lambda command: json.dumps([state]))
    monkeypatch.setattr(owner, "_allocated", lambda root: 0)

    def run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=1 if "inspect" in command else 0)

    class Process:
        returncode = None

        def __init__(self, command, **kwargs):
            Path(command[command.index("--cidfile") + 1]).write_text("c" * 64)

        def poll(self):
            return self.returncode

        def wait(self, **kwargs):
            self.returncode = 137

    monkeypatch.setattr(owner.subprocess, "run", run)
    monkeypatch.setattr(owner.subprocess, "Popen", Process)
    report = owner.run_pair(tmp_path, 2.45, "generated-owned", 120, time.time() - 1)
    assert report["guard_stop"] == "wall_deadline"
    assert report["original_command_exit"] == report["container_exit"] == 137
    assert not report["diagnostic_complete"]
    assert all(c == ["docker", "kill", "c" * 64] for c in calls if "kill" in c)
    assert any("kill" in c for c in calls)
