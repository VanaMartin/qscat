"""Check completed CO pilot outputs and collect native resonance fits.

Checks cover complete energy grids, finite nonnegative cross sections, agreement
of the two Pi components, and Psi4/UKRmol target energies. Native RESON fits are
reported in Hartree; a missing fit is reported explicitly, never as zero width.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np


def ground_cross_sections(path: Path) -> np.ndarray:
    """Join ground-state column blocks and exclude any excited-initial-state blocks."""
    rows = []
    blocks = []
    initial = None
    for line in path.read_text().splitlines():
        header = re.search(r"FOR INITIAL STATE\s+(\d+)", line)
        if header:
            if rows:
                blocks.append(np.loadtxt(io.StringIO("\n".join(rows)), ndmin=2))
                rows = []
            initial = int(header[1])
        elif initial == 1 and line.strip() and not line.lstrip().startswith("#"):
            rows.append(line)
    if rows:
        blocks.append(np.loadtxt(io.StringIO("\n".join(rows)), ndmin=2))
    if not blocks:
        raise ValueError(f"No ground-initial-state cross-section block: {path}")
    for block in blocks[1:]:
        np.testing.assert_allclose(block[:, :2], blocks[0][:, :2], atol=5e-5, rtol=0)
    return np.hstack([blocks[0], *(block[:, 2:] for block in blocks[1:])])


def native_resonances(path: Path) -> list[dict[str, float]]:
    """Read every native fitted position/width pair, converting Rydberg to Hartree."""
    fits = []
    pattern = (
        r"Fitted resonance parameters\s+Positions\s*/\s*Ryd\s*(.*?)"
        r"\s*Widths\s*/\s*Ryd\s*(.*?)\s*Background"
    )
    for positions, widths in re.findall(pattern, path.read_text(), flags=re.I | re.S):
        energies = [float(value.replace("D", "E")) for value in positions.split()]
        gammas = [float(value.replace("D", "E")) for value in widths.split()]
        for energy, gamma in zip(energies, gammas, strict=True):
            if energy > 0 and gamma > 0 and np.isfinite(energy + gamma):
                fits.append({"energy_hartree": energy / 2, "width_hartree": gamma / 2})
    return fits


def target_diagonalizer_properties(workdir: Path, config: dict) -> dict:
    """Reject a silent dense fallback or incomplete selected-root sector log."""
    if config.get("target_diagonalizer", "auto") == "auto":
        return {}
    sectors = {}
    for spin in ("singlet", "triplet"):
        counts = config.get(f"target_{spin}_roots", [config["target_roots"]] * 4)
        for irrep, count in zip(("A1", "B1", "B2", "A2"), counts, strict=True):
            if not count:
                continue
            log = workdir / f"output/CO/geom1/outputs/target.scatci.{spin}.{irrep}.out"
            text = log.read_text()
            requested = re.findall(r"Requested # of eigenpairs\s+(\d+)", text)
            if config["target_diagonalizer"] == "davidson-serial":
                completed = re.findall(
                    r"Davidson diagonalisation completed:\s+(\d+) iterations\s+"
                    r"(\d+) matrix vector multiplies IERR\s*=\s*(\d+)",
                    text,
                )
                valid = (
                    "Diagonalization done with Davidson" in text
                    and re.search(r"Sequential diagonalizations:\s*T", text) is not None
                    and len(completed) == 1
                    and int(completed[0][2]) == 0
                )
                diagnostics = (
                    {
                        "iterations": int(completed[0][0]),
                        "matrix_vector_multiplies": int(completed[0][1]),
                    }
                    if valid
                    else {}
                )
            else:
                completed = re.findall(r"stopped at it=(\d+)", text)
                valid = (
                    "KRYLOVSCHUR used as Diagonalizer" in text
                    and "Optimized SLEPC Matrix Format chosen" in text
                    and len(completed) == 1
                    and "Not all requested eigenpairs have converged" not in text
                    and "EIGEN-ENERGIES" in text
                )
                diagnostics = {"iterations": int(completed[0])} if valid else {}
            if not valid or requested != [str(count)]:
                raise ValueError(f"Incomplete or unexpected target diagonalizer in {log}")
            sectors[f"{spin}.{irrep}"] = {"requested_roots": count, **diagnostics}
    return {"target_diagonalizer": config["target_diagonalizer"], "target_solver_sectors": sectors}


def target_properties(workdir: Path, config: dict) -> dict:
    """Compare the independent target solvers, including every averaged QC root."""
    target_file = workdir / "output/CO/target.energies"
    target_table = np.loadtxt(target_file, ndmin=2)[0, 1:]
    target = target_table[0]
    header = next(line for line in target_file.read_text().splitlines() if "R_bohr" in line)
    labels = header.split()[2:]
    energies = dict(zip(labels, map(float, target_table), strict=True))
    excitations = {label: energy - target for label, energy in energies.items()}
    result = {
        "bond_length_bohr": config["bond_length"],
        "model": config["model"],
        "basis": config.get("basis", "aug-cc-pVDZ"),
        "orbitals": config.get("orbitals", "HF"),
        "active_orbitals": config.get("active_orbitals", [3, 1, 1, 0]),
        "frozen_orbitals": config.get("frozen_orbitals", 2),
        "neutral_energy_hartree": float(target),
        "target_excitations_hartree": excitations,
    }
    result.update(target_diagonalizer_properties(workdir, config))
    if config.get("orbitals") == "state-averaged":
        report = json.loads((workdir / "target.json").read_text())
        if not report["converged"]:
            raise ValueError("Unconverged state-averaged orbitals")
        molden = workdir / "output/CO/geom1/co.molden"
        if hashlib.sha256(molden.read_bytes()).hexdigest() != report["molden_sha256"]:
            raise ValueError("Imported Molden orbitals differ from the QC export")
        differences = {}
        for state in report["states"]:
            label = f"{state['spin']}.{state['irrep']}.{state['root']}"
            if label not in energies:
                raise ValueError(f"UKRmol did not compute averaged QC state {label}")
            differences[label] = energies[label] - state["energy_hartree"]
        np.testing.assert_allclose(list(differences.values()), 0, atol=1e-7, rtol=0)
        for spin in ("singlet", "triplet"):
            counts = config.get(f"target_{spin}_roots", [config["target_roots"]] * 4)
            for root in range(1, counts[1] + 1):
                np.testing.assert_allclose(
                    energies[f"{spin}.B1.{root}"],
                    energies[f"{spin}.B2.{root}"],
                    atol=1e-7,
                    rtol=0,
                )
        reference_energy = report["ground_energy_hartree"]
        # DENPROP's full-precision L=1,M=0 record is electronic position minus
        # nuclear charge position: the negative of the physical dipole.
        denprop = (workdir / "output/CO/geom1/outputs/target.denprop.out").read_text()
        dipole_match = re.findall(
            r"^\s*1\s+1\s+0\s+1\s+0\s+3\s+1\s+0\s+(\S+)\s*$",
            denprop,
            flags=re.M,
        )
        if len(dipole_match) != 1:
            raise ValueError("No unique DENPROP ground-state dipole record")
        native_dipole = -float(dipole_match[0].replace("D", "E"))
        dipole_difference = native_dipole - report["ground_dipole_au"][2]
        np.testing.assert_allclose(dipole_difference, 0, atol=1e-5, rtol=0)
        result.update(
            target_backend="pyscf",
            target_reference_energy_hartree=reference_energy,
            ground_state_dipole_au=report["ground_dipole_au"],
            ukrmol_ground_state_dipole_z_au=native_dipole,
            target_dipole_difference_au=dipole_difference,
            target_state_energy_differences_hartree=differences,
            target_state_max_energy_difference_hartree=max(map(abs, differences.values())),
        )
    else:
        psi4_output = (workdir / "output/CO/geom1/outputs/target.psi4.out").read_text()
        match = re.search(r"@(?:DF-)?RHF Final Energy:\s*(\S+)", psi4_output)
        if match is None or "Energy and wave function converged." not in psi4_output:
            raise ValueError("Psi4 did not report a converged RHF target energy")
        psi4_energy = float(match[1])
        reference_energy = psi4_energy
        dipoles = re.findall(r"Dipole Z\s*:\s*\S+\s+\S+\s+(\S+)", psi4_output)
        casscf = re.search(r"(?:MCSCF|CASSCF) Final Energy:\s*(\S+)", psi4_output)
        if config["model"] in ("SE", "SEP"):
            np.testing.assert_allclose(target, psi4_energy, atol=1e-7, rtol=0)
        elif config.get("orbitals", "HF") == "natural":
            if casscf is None:
                raise ValueError("Psi4 did not report a final CASSCF energy")
            reference_energy = float(casscf[1])
            np.testing.assert_allclose(target, reference_energy, atol=1e-7, rtol=0)
        elif not np.isfinite(target) or target >= psi4_energy:
            raise ValueError("CAS target energy must be finite and below the reference RHF energy")
        result.update(
            psi4_ground_state_dipole_au=float(dipoles[-1]) if dipoles else None,
            psi4_energy_hartree=reference_energy,
            psi4_rhf_energy_hartree=psi4_energy,
        )
    result["target_energy_difference_hartree"] = float(target - reference_energy)
    result["analyzer_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def qc_properties(workdir: Path, config: dict) -> dict:
    """Check a QC-only record without implying independent orbital-import agreement."""
    report = json.loads((workdir / "target.json").read_text())
    diagnostics = report["optimization"]
    if not (
        report["converged"]
        and diagnostics["orbital_converged"]
        and all(diagnostics["ci_converged"])
        and all(diagnostics["fresh_ci_converged"])
    ):
        raise ValueError("Unconverged QC-only target")
    expected = [
        (spin, irrep, root)
        for spin in ("singlet", "triplet")
        for irrep, count in zip(("A1", "B1", "B2", "A2"), config[f"sa_{spin}_roots"], strict=True)
        for root in range(1, count + 1)
    ]
    states = report["states"]
    if [(s["spin"], s["irrep"], s["root"]) for s in states] != expected:
        raise ValueError("QC-only roots differ from the requested orbital ensemble")
    energies = np.array([s["energy_hartree"] for s in states])
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(report["ground_dipole_au"])):
        raise ValueError("Nonfinite QC-only target properties")
    np.testing.assert_allclose([s["weight"] for s in states], 1 / len(states), atol=1e-14, rtol=0)
    np.testing.assert_allclose(
        energies.mean(), report["ensemble_energy_hartree"], atol=1e-8, rtol=0
    )
    np.testing.assert_allclose(
        diagnostics["fresh_state_energies_hartree"], energies, atol=1e-7, rtol=0
    )
    spins = [0 if spin == "singlet" else 2 for spin, _, _ in expected]
    np.testing.assert_allclose([s["spin_square"] for s in states], spins, atol=1e-6, rtol=0)
    np.testing.assert_allclose(diagnostics["fresh_spin_square"], spins, atol=1e-6, rtol=0)
    last = diagnostics["iterations"][-1]
    gradient = last.get("orbital_gradient_norm", last.get("orbital_ci_gradient_norm"))
    if (
        gradient is None
        or not np.isfinite(gradient)
        or gradient > config["target_gradient_tolerance"]
    ):
        raise ValueError("QC-only gradient exceeds its requested tolerance")
    if (
        not np.isfinite(report["mo_orthogonality_max_error"])
        or report["mo_orthogonality_max_error"] > 1e-9
    ):
        raise ValueError("QC-only orbitals are not orthonormal")
    values = {f"{s['spin']}.{s['irrep']}.{s['root']}": s["energy_hartree"] for s in states}
    np.testing.assert_allclose(
        values["singlet.A1.1"], report["ground_energy_hartree"], atol=1e-10, rtol=0
    )
    for spin in ("singlet", "triplet"):
        for root in range(1, config[f"sa_{spin}_roots"][1] + 1):
            np.testing.assert_allclose(
                values[f"{spin}.B1.{root}"], values[f"{spin}.B2.{root}"], atol=1e-7, rtol=0
            )
    if hashlib.sha256((workdir / "co.molden").read_bytes()).hexdigest() != report["molden_sha256"]:
        raise ValueError("QC-only Molden export differs from the target record")
    return {
        "qc_only": True,
        "bond_length_bohr": config["bond_length"],
        "basis": config["basis"],
        "target_backend": "pyscf",
        "ensemble_energy_hartree": report["ensemble_energy_hartree"],
        "ground_energy_hartree": report["ground_energy_hartree"],
        "ground_state_dipole_au": report["ground_dipole_au"],
        "target_state_energies_hartree": values,
        "analyzer_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "validation_scope": "QC solver/self-consistency; independent UKRmol import not performed",
    }


def neutral_properties(workdir: Path, config: dict) -> dict:
    """Check the neutral pilot's applicability and energy/density bookkeeping."""
    report = json.loads((workdir / "neutral.json").read_text())
    if report["basis"] != config["basis"]:
        raise ValueError("Neutral basis differs from its requested input")
    np.testing.assert_allclose(
        report["bond_length_bohr"], config["bond_length"], atol=1e-12, rtol=0
    )
    for flag in (
        "rhf_converged",
        "rhf_internally_stable",
        "rhf_externally_stable",
        "ccsd_converged",
        "lambda_converged",
    ):
        if not report[flag]:
            raise ValueError(f"Neutral pilot failed {flag}")
    energies = [
        report[k]
        for k in (
            "rhf_energy_hartree",
            "ccsd_energy_hartree",
            "triples_correction_hartree",
            "neutral_energy_hartree",
        )
    ]
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(report["ccsd_dipole_au"])):
        raise ValueError("Nonfinite neutral energy/dipole")
    if report["frozen_spatial_orbitals"] != 2 or report["correlated_electrons"] != 10:
        raise ValueError("Neutral pilot correlation space differs from its contract")
    diagnostics = [
        report["t1_frobenius_over_sqrt_correlated_electrons"],
        report["t1_largest_singular_value"],
    ]
    if not np.all(np.isfinite(diagnostics)) or min(diagnostics) < 0:
        raise ValueError("Invalid neutral amplitude diagnostics")
    np.testing.assert_allclose(report["density_electrons"], 14, atol=1e-7, rtol=0)
    np.testing.assert_allclose(
        report["ccsd_energy_hartree"] + report["triples_correction_hartree"],
        report["neutral_energy_hartree"],
        atol=1e-10,
        rtol=0,
    )
    if config["reference_check"]:
        for key in ("rhf_energy_hartree", "neutral_energy_hartree"):
            np.testing.assert_allclose(
                report[key], report["independent_reference"][key], atol=1e-7, rtol=0
            )
    return report | {
        "neutral_only": True,
        "analyzer_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "validation_scope": (
            "Correlated neutral pilot; basis and single-reference convergence not established"
        ),
    }


def analyze(workdir: Path) -> dict:
    """Validate a successful pilot and write a machine-readable result summary."""
    config = json.loads((workdir / "config.json").read_text())
    resources = json.loads((workdir / "resources.json").read_text())
    if resources["exit_code"] != 0:
        raise ValueError(f"Calculation failed; inspect {workdir / 'run.log'}")
    if config.get("qc_only", False) or config.get("neutral_only", False):
        result = (
            neutral_properties(workdir, config)
            if config.get("neutral_only")
            else qc_properties(workdir, config)
        )
        result["resources"] = resources
        (workdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        return result
    result = target_properties(workdir, config)
    if config.get("target_only", False):
        result.update(resources=resources, validation_scope="target solver/import consistency")
        (workdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        return result
    model = workdir / "output/CO"
    geometry = model / "geom1"
    energy = config["energy_start_ev"] + config["energy_step_ev"] * np.arange(config["energies"])
    phases = {}
    cross_sections = {}
    resonances = {}
    dimensions = {}
    for symmetry in ("B1", "B2"):
        phase = np.loadtxt(geometry / f"eigenph.doublet.{symmetry}", ndmin=2)
        cross = ground_cross_sections(geometry / f"xsec.doublet.{symmetry}")
        if phase.shape[0] != len(energy) or cross.shape[0] != len(energy):
            raise ValueError(f"Incomplete energy grid for {symmetry}")
        if not (np.all(np.isfinite(phase)) and np.all(np.isfinite(cross))):
            raise ValueError(f"Nonfinite scattering outputs for {symmetry}")
        np.testing.assert_allclose(phase[:, 0], energy, atol=5e-4, rtol=0)
        np.testing.assert_allclose(cross[:, 1], energy, atol=5e-5, rtol=0)
        if np.any(cross[:, 2:] < -1e-10):
            raise ValueError(f"Negative cross section for {symmetry}")
        np.testing.assert_allclose(cross[:, 2], cross[:, 3:].sum(axis=1), atol=1e-6, rtol=1e-6)
        phases[symmetry] = phase[:, 1]
        cross_sections[symmetry] = cross[:, 2]
        resonances[symmetry] = native_resonances(
            geometry / f"outputs/scattering.reson.doublet.{symmetry}.out"
        )
        scatci = (geometry / f"outputs/scattering.scatci.doublet.{symmetry}.out").read_text()
        size = re.search(r"MOCSF\s*=\s*(\d+)\s+dimension final Hamiltonian", scatci)
        if size is None:
            raise ValueError(f"No final Hamiltonian dimension reported for {symmetry}")
        dimensions[symmetry] = int(size[1])
    difference = phases["B1"] - phases["B2"]
    phase_error = np.remainder(difference + np.pi / 2, np.pi) - np.pi / 2
    np.testing.assert_allclose(phase_error, 0, atol=2e-3, rtol=0)
    np.testing.assert_allclose(cross_sections["B1"], cross_sections["B2"], atol=1e-6, rtol=1e-5)
    result.update(
        {
            "energy_points_per_symmetry": len(energy),
            "hamiltonian_dimensions": dimensions,
            "energy_range_ev": [float(energy[0]), float(energy[-1])],
            "pi_phase_max_difference_rad": float(np.max(np.abs(phase_error))),
            "pi_cross_section_max_absolute_difference_bohr2": float(
                np.max(np.abs(cross_sections["B1"] - cross_sections["B2"]))
            ),
            "native_resonances": resonances,
            "resources": resources,
            "validation_scope": (
                "pipeline and symmetry consistency; basis convergence not established"
            ),
        }
    )
    (workdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    """Analyze an existing run directory without launching another calculation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.workdir.resolve()), indent=2))


if __name__ == "__main__":
    main()
