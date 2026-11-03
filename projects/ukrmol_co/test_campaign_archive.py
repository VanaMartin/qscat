"""Keep historical evidence separate from current source and test collection."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "projects/ukrmol_co/qualification-evidence/campaign-history"


def test_history_pointer_authenticates_the_tracked_original_manifest():
    pointer = json.loads((ARCHIVE / "artifacts.json").read_text())
    raw = (ARCHIVE / "manifest.json").read_bytes()
    pin = pointer["artifacts"]["manifest.json"]
    assert len(raw) == pin["bytes"]
    assert hashlib.sha256(raw).hexdigest() == pin["sha256"]
    manifest = json.loads(raw)
    assert manifest["git_sha"] == pointer["git_sha"]
    assert manifest["members"] == 363
    assert pointer["artifacts"]["campaign-history.tar.gz"]["sha256"] == (
        "1c09411aebf941fe5ecb8688af49749dee1600fea85f652bb6f310ec7eb91bea"
    )


def test_pytest_does_not_collect_historical_source_capsules(tmp_path):
    archived = tmp_path / "state-averaged-evidence" / "owner-source"
    archived.mkdir(parents=True)
    (archived / "test_historical.py").write_text(
        "raise RuntimeError('historical source imported')\n"
    )
    (tmp_path / "test_current.py").write_text("def test_current():\n    assert True\n")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "--collect-only",
            "-q",
            "-c",
            str(ROOT / "pyproject.toml"),
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        env=os.environ | {"PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"},
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "test_current" in result.stdout
    assert "test_historical" not in result.stdout
