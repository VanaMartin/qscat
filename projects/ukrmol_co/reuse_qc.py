"""Authenticate and replay accepted SA orbitals without invoking a QC solver."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

FILES = (
    "config.json",
    "target.json",
    "target-diagnostics.json",
    "resources.json",
    "co.molden",
    "co.casscf.chk",
)
PROTOCOL = (
    "bond_length",
    "basis",
    "model",
    "orbitals",
    "frozen_orbitals",
    "active_orbitals",
    "sa_singlet_roots",
    "sa_triplet_roots",
)
IRREPS = ("A1", "B1", "B2", "A2")


def authenticated_qc(manifest: Path, expected_sha256: str, config: dict) -> dict[str, bytes]:
    """Reject changed bytes/protocol or unsuccessful original QC before any writes.

    The manifest has a source_directory and exactly six file SHA256s. Relative
    source directories resolve beside the manifest. Native retained root counts
    are deliberately not the orbital ensemble. This is provenance replay, not
    a new wavefunction calculation or native target qualification.
    """
    payload = manifest.read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise ValueError("QC reuse manifest digest differs")
    spec = json.loads(payload)
    if set(spec["files"]) != set(FILES):
        raise ValueError("QC reuse needs the complete six-file bundle")
    source = Path(spec["source_directory"])
    if not source.is_absolute():
        source = manifest.parent / source
    data = {name: (source / name).read_bytes() for name in FILES}
    for name, contents in data.items():
        if hashlib.sha256(contents).hexdigest() != spec["files"][name]:
            raise ValueError(f"QC source digest differs: {name}")
    original = json.loads(data["config.json"])
    for key in PROTOCOL:
        if original[key] != config[key]:
            raise ValueError(f"QC orbital protocol differs: {key}")
    if (
        config["model"] != "CAS-A"
        or config["orbitals"] != "state-averaged"
        or config["frozen_orbitals"] != 2
    ):
        raise ValueError("QC reuse requires the two-core SA CAS target")
    if json.loads(data["resources.json"])["exit_code"] != 0:
        raise ValueError("Original QC process failed")
    report = json.loads(data["target.json"])
    diagnostics = json.loads(data["target-diagnostics.json"])
    if report["converged"] is not True or report["optimization"] != diagnostics:
        raise ValueError("QC convergence report/diagnostics differ")
    if report["method"] != "equal-weight multi-spin/multi-irrep SA-CASSCF":
        raise ValueError("Unexpected orbital objective")
    counts = [original["sa_singlet_roots"], original["sa_triplet_roots"]]
    sectors = sum(n > 0 for group in counts for n in group)
    if diagnostics["orbital_converged"] is not True or any(
        diagnostics[key] != [True] * sectors for key in ("ci_converged", "fresh_ci_converged")
    ):
        raise ValueError("Original orbital/fresh-CI convergence failed")

    def within(value, budget):
        return math.isfinite(value) and abs(value) <= budget

    if not within(diagnostics["lowest_root_max_energy_difference_hartree"], 1e-7):
        raise ValueError("Original fresh lowest-root audit failed")
    if not within(report["mo_orthogonality_max_error"], 1e-9):
        raise ValueError("Original MO orthogonality failed")
    atoms = report["atoms_bohr"]
    radius = config["bond_length"]
    cz = -radius * 15.9994 / (12.0110 + 15.9994)
    if [atom[0] for atom in atoms] != ["C", "O"] or any(
        not within(value - expected, 1e-10)
        for atom, z in zip(atoms, [cz, cz + radius], strict=True)
        for value, expected in zip(atom[1:], [0, 0, z], strict=True)
    ):
        raise ValueError("QC geometry/origin differs")
    expected_labels = [
        (spin, irrep, root + 1)
        for spin, group in zip(("singlet", "triplet"), counts, strict=True)
        for irrep, count in zip(IRREPS, group, strict=True)
        for root in range(count)
    ]
    states = report["states"]
    if [(s["spin"], s["irrep"], s["root"]) for s in states] != expected_labels:
        raise ValueError("QC ensemble inventory differs")
    if not states or any(
        len(diagnostics[key]) != len(states)
        for key in (
            "fresh_spin_square",
            "state_energies_hartree",
            "fresh_state_energies_hartree",
        )
    ):
        raise ValueError("Incomplete QC state diagnostics")
    for index, state in enumerate(states):
        spin_squared = 0 if state["spin"] == "singlet" else 2
        if not all(
            (
                within(state["weight"] - 1 / len(states), 1e-14),
                within(state["spin_square"] - spin_squared, 1e-6),
                within(diagnostics["fresh_spin_square"][index] - spin_squared, 1e-6),
                within(
                    state["energy_hartree"] - diagnostics["state_energies_hartree"][index], 1e-12
                ),
                within(
                    state["energy_hartree"] - diagnostics["fresh_state_energies_hartree"][index],
                    1e-7,
                ),
            )
        ):
            raise ValueError("QC state weight/spin/energy audit failed")
    ground = states[0]["energy_hartree"]
    if states[0]["irrep"] != "A1" or not within(report["ground_energy_hartree"] - ground, 1e-12):
        raise ValueError("QC ground total differs")
    if any(s["energy_hartree"] < ground - 1e-8 for s in states):
        raise ValueError("QC ground state is not lowest")
    orbitals = report["orbitals"]
    nactive = sum(config["active_orbitals"])
    if [o["space"] for o in orbitals] != ["core"] * 2 + ["active"] * nactive + ["virtual"] * (
        len(orbitals) - 2 - nactive
    ):
        raise ValueError("QC orbital block order differs")
    if [o["irrep"] for o in orbitals[:2]] != ["A1", "A1"] or any(
        sum(o["irrep"] == irrep for o in orbitals[2 : 2 + nactive]) != count
        for irrep, count in zip(IRREPS, config["active_orbitals"], strict=True)
    ):
        raise ValueError("QC core/active symmetry differs")
    if report["molden_sha256"] != spec["files"]["co.molden"]:
        raise ValueError("QC report references different Molden bytes")
    if not data["co.casscf.chk"].startswith(b"\x89HDF\r\n\x1a\n"):
        raise ValueError("QC checkpoint is not HDF5")
    if not all(
        section in data["co.molden"]
        for section in (b"[Molden Format]", b"[Atoms]", b"[GTO]", b"[MO]")
    ):
        raise ValueError("QC Molden sections are missing")
    return data


def replay_qc(config: dict, molden_path: Path) -> dict:
    """Copy exact accepted exports into a fresh native run; never recompute."""
    manifest = Path(config["reuse_qc_manifest"])
    data = authenticated_qc(manifest, config["reuse_qc_sha256"], config)
    root = Path(config["workdir"])
    destinations = {
        root / name: data[name]
        for name in (
            "target.json",
            "target-diagnostics.json",
            "co.casscf.chk",
        )
    }
    destinations[molden_path] = data["co.molden"]
    receipt = root / "qc-reuse.json"
    if receipt.exists() or any(path.exists() for path in destinations):
        raise ValueError("QC replay refuses to overwrite existing exports")
    for path, contents in destinations.items():
        with path.open("xb") as stream:
            stream.write(contents)
    record = dict(
        manifest_sha256=config["reuse_qc_sha256"],
        new_qc_calls=0,
        original_diagnostics_reused=True,
        files={name: hashlib.sha256(b).hexdigest() for name, b in data.items()},
    )
    with receipt.open("x") as stream:
        stream.write(json.dumps(record, indent=2) + "\n")
    return json.loads(data["target.json"])
