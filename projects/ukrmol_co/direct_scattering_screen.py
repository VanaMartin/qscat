"""Necessary diagonal-family bounds for saved physical-channel CO S matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.saved_matrices import read_interface, read_k, scattering_matrix, unit


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def diagonal_bounds(s, phase_budget=0.05, margin=1e-10):
    """A rigorous off-diagonal lower bound and one constructive unitary upper bound."""
    s = np.asarray(s, dtype=complex)
    if s.ndim != 2 or s.shape[0] != s.shape[1] or not len(s) or not np.isfinite(s).all():
        raise ValueError("A finite square S matrix is required")
    if not 0 < phase_budget < np.pi / 2 or not 0 < margin < 2 * np.sin(phase_budget):
        raise ValueError("Invalid prospective phase budget or numerical margin")
    defect = float(np.max(abs(s.conj().T @ s - np.eye(len(s)))))
    if defect > 1e-10:
        raise ValueError("Nonunitary saved/control S matrix")
    off_diagonal = s.copy()
    np.fill_diagonal(off_diagonal, 0)
    lower = float(np.max(np.linalg.norm(off_diagonal, axis=0)))
    diagonal = np.diag(s)
    candidate = np.ones(len(s), dtype=complex)
    nonzero = abs(diagonal) > 0
    candidate[nonzero] = diagonal[nonzero] / abs(diagonal[nonzero])
    upper = float(np.linalg.norm(s - np.diag(candidate), ord=2))
    threshold = float(2 * np.sin(phase_budget))
    decision = (
        "reject"
        if lower > threshold + margin
        else "admit"
        if upper <= threshold - margin
        else "inconclusive"
    )
    return dict(
        lower_chord_bound=lower,
        constructive_upper_chord_bound=upper,
        lower_relative_eigenphase_bound_rad=float(np.arcsin(min(1, lower / 2))),
        constructive_relative_eigenphase_error_rad=float(np.arcsin(min(1, upper / 2))),
        chord_budget=threshold,
        decision=decision,
        unitarity_defect=defect,
        candidate_diagonal_real=candidate.real.tolist(),
        candidate_diagonal_imag=candidate.imag.tolist(),
        common_potential_verdict=None,
        width_identification_verdict=None,
    )


def execute(root, output, contract):
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    report = dict(
        cases=[],
        bare_diagonal_full_s_family_verdict=None,
        common_potential_verdict=None,
        width_identification_verdict=None,
        full_sweep_released=False,
        native_calls=0,
        source_scope="Saved K/Cayley algebra, not independent outer propagation",
    )
    raw = root / contract["raw_directory"]

    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed prerequisite: {name}")
        manifest = json.loads((raw / "input-manifest.json").read_text())
        if len(manifest) != 123:
            raise ValueError("Changed original input inventory")
        for name, want in manifest.items():
            if digest(raw / name) != want:
                raise ValueError(f"Changed original matrix input: {name}")
        return manifest

    try:
        manifest = verify()
        report["original_inputs_verified"] = len(manifest)
        previous = json.loads((root / contract["qualified_report"]).read_text())
        if previous["saved_matrix_reproduction_verdict"] is not False:
            raise ValueError("Original aggregate failure must remain false")
        for anchor in contract["anchors"]:
            for sector in contract["sectors"]:
                name = f"{anchor}/nine/{sector}"
                prior = next(row for row in previous["cases"] if row["case"] == name)
                if prior["saved_matrix_reproduction_verdict"] is not True:
                    raise ValueError("Unqualified nine-A2 prerequisite")
                directory = raw / name
                required = [
                    "rsolve.inp",
                    "fort.10",
                    "fort.21",
                    unit(directory / "rsolve.inp", "lukmt"),
                ]
                if any(f"{name}/{file}" not in manifest for file in required):
                    raise ValueError("Case is not completely checksum-bound")
                energy, matrices, header = read_k(directory / required[-1])
                interface = read_interface(directory / "fort.10", directory / "fort.21")
                if header != prior["header"] or interface != prior["interface"]:
                    raise ValueError("Metadata changed from qualified physical-channel case")
                if matrices.shape != (99, 6, 6):
                    raise ValueError("Expected99 six-open-channel records")
                labels = interface["channel_labels"][:6]
                if [row[2] for row in labels] != [1, 2, 3, 3, 4, 4] or any(
                    row[1] != 1 for row in labels
                ):
                    raise ValueError("Unexpected open partial waves or target ordering")
                requested = energy / contract["native_rydberg_per_requested_ev"]
                np.testing.assert_allclose(requested, 0.1 + 0.05 * np.arange(99), atol=5e-5, rtol=0)
                chosen = np.flatnonzero(np.arange(99) % 5 != 4)
                held = np.flatnonzero(np.arange(99) % 5 == 4)
                points = []
                saved_s = []
                for index in chosen:
                    k = matrices[index]
                    s = scattering_matrix(k)
                    values, vectors = np.linalg.eigh(k)
                    alternate = (vectors * np.exp(2j * np.arctan(values))) @ vectors.T
                    np.testing.assert_allclose(s, alternate, atol=1e-10, rtol=0)
                    np.testing.assert_allclose(s, s.T, atol=1e-10, rtol=0)
                    points.append(
                        dict(
                            index=int(index),
                            requested_energy_ev=float(requested[index]),
                            **diagonal_bounds(s, contract["phase_budget_rad"], contract["margin"]),
                        )
                    )
                    saved_s.append(s)
                np.savez(
                    output / f"{anchor}-{sector}-development.npz",
                    indices=chosen,
                    energy_rydberg=energy[chosen],
                    s=np.stack(saved_s),
                )
                row = dict(
                    case=name,
                    physical_open_channel_labels=labels,
                    development_indices=chosen.tolist(),
                    reserved_indices=held.tolist(),
                    points=points,
                    rejected_points=sum(p["decision"] == "reject" for p in points),
                    inconclusive_points=sum(p["decision"] == "inconclusive" for p in points),
                )
                report["cases"].append(row)
        decisions = [point["decision"] for row in report["cases"] for point in row["points"]]
        report["development_points_analyzed"] = len(decisions)
        report["reserved_points_not_analyzed"] = 6 * 19
        report["bare_diagonal_full_s_family_verdict"] = (
            False if "reject" in decisions else None if "inconclusive" in decisions else True
        )
        report["scientific_exit"] = (
            1 if "reject" in decisions else 2 if "inconclusive" in decisions else 0
        )
        verify()
        report["input_integrity_after"] = True
        report["wall_seconds"] = time.monotonic() - start
        report["peak_process_rss_platform_units"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        return report["scientific_exit"]
    except BaseException as error:
        report["error"] = repr(error)
        try:
            verify()
            report["input_integrity_after"] = True
        except BaseException as integrity_error:
            report["input_integrity_after"] = False
            report["integrity_error"] = repr(integrity_error)
        report["wall_seconds"] = time.monotonic() - start
        (output / "failure.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(execute(args.root, args.output, json.loads(args.contract.read_text())))
