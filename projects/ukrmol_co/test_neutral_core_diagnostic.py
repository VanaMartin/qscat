"""Generated electronic objects test orchestration without importing PySCF."""

import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
import pytest

from projects.ukrmol_co import neutral_core_diagnostic as diagnostic
from projects.ukrmol_co.analyze import neutral_properties


@pytest.fixture
def backend(monkeypatch):
    calls = []
    molecule_requests = []
    defect = {}

    class Molecule:
        nelectron, spin, charge, cart = 14, 0, 0, False

        def nao_nr(self):
            return defect.get("nao", 218)

    class Reference:
        def __init__(self, mol):
            self.mol = mol
            self.converged = True
            self.e_tot = -10.0
            self.mo_coeff = np.eye(218)
            self.mo_occ = np.r_[np.full(7, 2.0), np.zeros(211)]
            self.mo_energy = np.arange(218, dtype=float)

        def kernel(self):
            calls.append("rhf")
            self.converged = not defect.get("rhf")

        def stability(self, **kwargs):
            calls.append("stability")
            return None, None, not defect.get("internal"), not defect.get("external")

        def get_ovlp(self):
            return np.eye(218)

    class Correlation:
        def __init__(self, hf, frozen):
            calls.append(("cc", frozen, id(hf)))
            self._scf = hf
            self.frozen = frozen
            self._nmo = self._nocc = None
            self.nocc = 7 - (frozen or 0)
            self.nmo = 218 - (frozen or 0)
            if defect.get("nocc"):
                self.nocc += 1
            if defect.get("override"):
                self._nmo = self.nmo
            self.t1 = np.zeros((self.nocc, 211))
            self.t2 = np.broadcast_to(0.0, (self.nocc, self.nocc, 211, 211))
            self.e_tot = -11.0 - (0.125 if frozen is None else 0)
            self.converged = True
            self.converged_lambda = True

        def get_frozen_mask(self):
            mask = np.ones(218, dtype=bool)
            mask[: self.frozen or 0] = False
            if defect.get("mask"):
                mask[-1] = False
            return mask

        def kernel(self):
            calls.append(("kernel", self.frozen))
            self.converged = not defect.get("ccsd")
            if defect.get("nan"):
                self.e_tot = np.nan
            if defect.get("mutation"):
                self._scf.mo_coeff[0, 0] = -1
            if defect.get("amplitude"):
                self.t1[0, 0] = np.nan

        def dump_chk(self):
            pass

        def ccsd_t(self):
            calls.append(("triples", self.frozen))
            return np.nan if defect.get("triples") else -0.01

        def solve_lambda(self):
            calls.append(("lambda", self.frozen))
            self.converged_lambda = not defect.get("lambda")
            if defect.get("ae_only") and self.frozen is None:
                self.converged_lambda = False

        def make_rdm1(self, **kwargs):
            electrons = 10 if defect.get("density") else 14
            return np.eye(218) * electrons / 218

    def molecule(**kwargs):
        molecule_requests.append(kwargs)
        return Molecule()

    value = SimpleNamespace(
        __version__="2.11.0",
        gto=SimpleNamespace(M=molecule),
        scf=SimpleNamespace(RHF=Reference),
        cc=SimpleNamespace(CCSD=Correlation),
        lib=SimpleNamespace(num_threads=lambda n: None),
        calls=calls,
        defect=defect,
        molecule_requests=molecule_requests,
    )
    monkeypatch.setattr(diagnostic, "_basis_identity", lambda b: diagnostic.BASIS_SHA256)
    return value


def config(tmp_path):
    return {"workdir": str(tmp_path), "bond_length": 2.45, "ranks": 4, "memory_mb": 20000}


def test_pair_reuses_one_reference_and_stays_ineligible(tmp_path, backend):
    result = diagnostic.calculate_pair(config(tmp_path), backend=backend)
    assert result["status"] == "numerically_complete_diagnostic"
    assert result["eligible_for_delivered_curve"] is False
    assert result["diagnostic_only"] is True
    rows = result["core_treatments"]
    assert [r["correlated_electrons"] for r in rows] == [10, 14]
    assert rows[0]["rhf_reference_sha256"] == rows[1]["rhf_reference_sha256"]
    assert backend.calls.count("rhf") == backend.calls.count("stability") == 1
    constructions = [c for c in backend.calls if isinstance(c, tuple) and c[0] == "cc"]
    assert [c[1] for c in constructions] == [2, None]
    assert constructions[0][2] == constructions[1][2]
    request = backend.molecule_requests[0]
    assert request["basis"] == diagnostic.BASIS
    assert request["unit"] == "Bohr" and request["symmetry"] == "C2v"
    assert request["spin"] == request["charge"] == 0 and request["cart"] is False
    assert request["max_memory"] == 20000
    np.testing.assert_allclose(
        request["atom"][1][3] - request["atom"][0][3], 2.45, atol=1e-14, rtol=0
    )
    assert not (tmp_path / "neutral.json").exists()
    assert diagnostic._pair_complete(
        tmp_path / "diagnostic-pair.json", diagnostic.source_identity()
    )
    with pytest.raises(FileNotFoundError):
        neutral_properties(tmp_path, {"basis": diagnostic.BASIS, "bond_length": 2.45})
    np.testing.assert_allclose(rows[0]["total_energy_hartree"], -11.01, atol=1e-13, rtol=0)


@pytest.mark.parametrize(
    "defect",
    [
        "rhf",
        "internal",
        "external",
        "ccsd",
        "nan",
        "triples",
        "lambda",
        "density",
        "mask",
        "mutation",
        "amplitude",
        "nocc",
        "override",
    ],
)
def test_first_failure_stops_without_second_core_treatment(tmp_path, backend, defect):
    backend.defect[defect] = True
    with pytest.raises(ValueError):
        diagnostic.calculate_pair(config(tmp_path), backend=backend)
    result = json.loads((tmp_path / "diagnostic-pair.json").read_text())
    assert result["status"] == "failed_diagnostic"
    assert result["eligible_for_delivered_curve"] is False
    assert not any(c[:2] == ("cc", None) for c in backend.calls if isinstance(c, tuple))
    if defect in {"rhf", "internal", "external"}:
        assert not any(isinstance(c, tuple) and c[0] == "cc" for c in backend.calls)


def test_failed_all_electron_retains_successful_frozen_core(tmp_path, backend):
    backend.defect["ae_only"] = True
    with pytest.raises(ValueError):
        diagnostic.calculate_pair(config(tmp_path), backend=backend)
    rows = json.loads((tmp_path / "diagnostic-pair.json").read_text())["core_treatments"]
    assert rows[0]["status"] == "numerically_complete_diagnostic"
    assert rows[1]["status"] != "numerically_complete_diagnostic"


@pytest.mark.parametrize("radius", [2.55, float("nan"), True, 0, "2.45"])
def test_config_rejects_other_geometries_before_electronic_calls(tmp_path, backend, radius):
    c = config(tmp_path) | {"bond_length": radius}
    with pytest.raises(ValueError):
        diagnostic.calculate_pair(c, backend=backend)
    assert backend.calls == []


def test_preflight_has_no_scf_or_cc_calls(tmp_path, backend):
    result = diagnostic.preflight(config(tmp_path), backend=backend)
    assert result["spherical_AOs"] == 218
    assert result["eligible_for_delivered_curve"] is False
    assert backend.calls == []


@pytest.mark.parametrize("defect", ["nao", "version"])
def test_preflight_rejects_runtime_mismatch(tmp_path, backend, defect):
    if defect == "nao":
        backend.defect["nao"] = 219
    else:
        backend.__version__ = "other"
    with pytest.raises(ValueError):
        diagnostic.calculate_pair(config(tmp_path), backend=backend)
    assert backend.calls == []


def test_source_mismatch_stops_before_backend(tmp_path, backend):
    c = config(tmp_path) | {"source_identity": {"wrong": "source"}}
    with pytest.raises(ValueError, match="source identity"):
        diagnostic.calculate_pair(c, backend=backend)
    assert backend.calls == []


def test_basis_digest_is_checked_before_molecule_construction(tmp_path, monkeypatch):
    (tmp_path / "__init__.py").touch()
    data = b"generated basis identity fixture\n"
    (tmp_path / "aug-cc-pwcvqz.dat").write_bytes(data)
    backend = SimpleNamespace(
        gto=SimpleNamespace(basis=SimpleNamespace(__file__=str(tmp_path / "__init__.py")))
    )
    monkeypatch.setattr(diagnostic, "BASIS_SHA256", hashlib.sha256(data).hexdigest())
    assert diagnostic._basis_identity(backend) == diagnostic.BASIS_SHA256
    (tmp_path / "aug-cc-pwcvqz.dat").write_bytes(data + b"changed")
    with pytest.raises(ValueError, match="basis digest"):
        diagnostic._basis_identity(backend)


def test_bounded_launcher_preserves_original_failure(tmp_path):
    work = tmp_path / "pair"
    code = diagnostic.launch_pair(
        work,
        config(work),
        command=[sys.executable, "-c", "print('inert child'); raise SystemExit(7)"],
    )
    assert code == 7
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["original_child_exit"] == 7
    assert report["stop_reason"] is None
    assert "inert child" in (work / "run.log").read_text()


def test_bounded_launcher_stops_only_its_sleeping_child(tmp_path):
    work = tmp_path / "pair"
    code = diagnostic.launch_pair(
        work,
        config(work),
        command=[sys.executable, "-c", "import time; time.sleep(20)"],
        wall_seconds=0.03,
        poll_seconds=0.01,
    )
    assert code == 124
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["stop_reason"] == "overall_wall_cap"
    assert report["original_child_exit"] != 0


@pytest.mark.parametrize("phase", ["reference", "fc", "ae", "invalid"])
def test_bounded_launcher_enforces_phase_identity_and_caps(tmp_path, phase):
    work = tmp_path / "pair"
    script = (
        "import json,time; from pathlib import Path; "
        "Path('active-stage.json').write_text(json.dumps("
        f"dict(phase={phase!r}, started_monotonic=time.monotonic()-1600))); "
        "time.sleep(20)"
    )
    code = diagnostic.launch_pair(
        work, config(work), command=[sys.executable, "-c", script], poll_seconds=0.01
    )
    assert code == 124
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["stop_reason"] == (
        "invalid_stage_metadata" if phase == "invalid" else f"{phase}_wall_cap"
    )


def test_zero_exit_without_scientific_output_is_not_completion(tmp_path):
    work = tmp_path / "pair"
    assert diagnostic.launch_pair(work, config(work), command=[sys.executable, "-c", "pass"]) == 2
    result = json.loads((work / "diagnostic-resources.json").read_text())
    assert result["original_child_exit"] == 0
    assert result["stop_reason"] == "missing_or_invalid_pair_output"


def test_launcher_does_not_overwrite_prior_work(tmp_path):
    marker = tmp_path / "user.txt"
    marker.write_text("retained")
    with pytest.raises(FileExistsError):
        diagnostic.launch_pair(tmp_path, config(tmp_path), command=[sys.executable, "-c", "pass"])
    assert marker.read_text() == "retained"


@pytest.mark.parametrize("resource", ["disk", "memory"])
def test_launcher_enforces_resource_caps(tmp_path, monkeypatch, resource):
    if resource == "disk":
        monkeypatch.setattr(diagnostic, "allocated_bytes", lambda p: 21 * 2**30)
    else:
        monkeypatch.setattr(diagnostic, "cgroup_bytes", lambda p: 25 * 2**30)
    work = tmp_path / "pair"
    code = diagnostic.launch_pair(
        work, config(work), command=[sys.executable, "-c", "import time; time.sleep(20)"]
    )
    assert code == 124
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["stop_reason"] == (
        "allocated_disk_cap" if resource == "disk" else "container_memory_cap"
    )


def test_parent_rejects_source_changed_during_child(tmp_path, monkeypatch):
    identities = iter([{"source": "before"}, {"source": "after"}])
    monkeypatch.setattr(diagnostic, "source_identity", lambda: next(identities))
    work = tmp_path / "pair"
    assert diagnostic.launch_pair(work, config(work), command=[sys.executable, "-c", "pass"]) == 2
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["original_child_exit"] == 0
    assert report["stop_reason"] == "source_changed"


@pytest.mark.parametrize(
    "field", ["density_electrons", "total_energy_hartree", "retained_virtual_orbitals"]
)
def test_pair_reader_rejects_corrupted_completion(tmp_path, backend, field):
    value = diagnostic.calculate_pair(config(tmp_path), backend=backend)
    value["core_treatments"][1][field] = 123
    (tmp_path / "diagnostic-pair.json").write_text(json.dumps(value))
    assert not diagnostic._pair_complete(
        tmp_path / "diagnostic-pair.json", diagnostic.source_identity()
    )


def test_launcher_accepts_complete_generated_diagnostic(tmp_path, backend):
    fixture = tmp_path / "generated"
    fixture.mkdir()
    value = diagnostic.calculate_pair(config(fixture), backend=backend)
    script = (
        f"from pathlib import Path; Path('diagnostic-pair.json').write_text({json.dumps(value)!r})"
    )
    work = tmp_path / "launched"
    assert diagnostic.launch_pair(work, config(work), command=[sys.executable, "-c", script]) == 0
    report = json.loads((work / "diagnostic-resources.json").read_text())
    assert report["original_child_exit"] == 0 and report["stop_reason"] is None
    assert report["eligible_for_delivered_curve"] is False


def test_cli_requires_explicit_execution_before_creating_directory(tmp_path):
    work = tmp_path / "not-admitted"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.neutral_core_diagnostic",
            "--bond-length",
            "2.45",
            "--workdir",
            str(work),
        ],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 2
    assert "Explicit reviewed execution" in result.stderr
    assert not work.exists()


def test_terminating_launcher_cleans_owned_child_group(tmp_path):
    work = tmp_path / "pair"
    child_script = (
        "import os,time; from pathlib import Path; "
        "Path('child.pid').write_text(str(os.getpid())); time.sleep(20)"
    )
    parent_script = (
        "from projects.ukrmol_co.neutral_core_diagnostic import launch_pair; "
        f"launch_pair({str(work)!r}, {config(work)!r}, "
        f"command=[{sys.executable!r}, '-c', {child_script!r}])"
    )
    process = subprocess.Popen([sys.executable, "-c", parent_script], start_new_session=True)
    child_pid = None
    try:
        deadline = time.monotonic() + 5
        while not (work / "child.pid").exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        assert (work / "child.pid").exists()
        child_pid = int((work / "child.pid").read_text())
        process.terminate()
        process.wait(timeout=10)
        report = json.loads((work / "diagnostic-resources.json").read_text())
        assert report["original_child_exit"] != 0
        assert report["stop_reason"] == "launcher_exception"
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
    finally:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait(timeout=5)
        if child_pid is not None:
            try:
                os.killpg(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
