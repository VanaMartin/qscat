"""Replay saved physical ensemble overlaps without a CI or native calculation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.ensemble_frame_overlaps import diagnose_overlap, full_minor_overlaps


def replay(root, comparison, contract):
    """Use full occupied minors instead of the producer's Schur contraction."""
    for name, want in contract["input_sha256"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != want:
            raise ValueError(f"Changed replay prerequisite: {name}")
    report = json.loads((comparison / "result.json").read_text())
    assert report["all_eight_sectors_complete"] is True and len(report["sectors"]) == 8

    def load(name):
        return np.load(comparison / (name + ".npy"), allow_pickle=False)

    a, b = load("old.mo"), load("new.mo")
    for name, c in (("old", a), ("new", b)):
        own = c.conj().T @ load(name + ".ao") @ c
        np.testing.assert_allclose(own, np.eye(13), atol=1e-9, rtol=0)
        np.testing.assert_allclose(own, load(name + ".self"), atol=1e-12, rtol=0)
    cross, back = load("cross.ao"), load("reverse.ao")
    np.testing.assert_allclose(cross, back.conj().T, atol=1e-12, rtol=0)
    metric = a.conj().T @ cross @ b
    np.testing.assert_allclose(metric, load("cross.mo"), atol=1e-12, rtol=0)
    np.testing.assert_allclose(b.conj().T @ back @ a, load("reverse.mo"), atol=1e-12, rtol=0)
    rows = []
    assert {s["sector"] for s in report["sectors"]} == {
        f"{spin}.{irrep}" for spin in ("singlet", "triplet") for irrep in ("A1", "B1", "B2", "A2")
    }
    for sector in report["sectors"]:
        label = sector["sector"]
        left = np.load(root / contract["old_cache"] / (label + ".ci.npy"), allow_pickle=False)
        right = load(label + ".ci")
        assert left.shape[0] == 12
        assert right.shape[0] == (13 if label.endswith("A1") else 8)
        for ci in (left, right):
            flat = ci.reshape(len(ci), -1)
            np.testing.assert_allclose(flat.conj() @ flat.T, np.eye(len(ci)), atol=1e-9, rtol=0)
        electrons = (5, 5) if label.startswith("singlet") else (6, 4)
        rebuilt = full_minor_overlaps(left, right, metric, 2, electrons)
        saved = load(label + ".overlap")
        np.testing.assert_allclose(rebuilt, saved, atol=1e-10, rtol=0)
        np.testing.assert_allclose(rebuilt, load(label + ".full-minors"), atol=1e-10, rtol=0)
        np.testing.assert_allclose(rebuilt, load(label + ".reverse").conj().T, atol=1e-10, rtol=0)
        retained = rebuilt[:5, :5]
        # Eigenspectrum of A†A is a second route to the producer's squared SVD.
        weights = np.linalg.eigvalsh(retained.conj().T @ retained)
        np.testing.assert_allclose(
            weights[::-1],
            sector["comparison"]["retained_squared_principal_overlaps"],
            atol=1e-10,
            rtol=0,
        )
        verdict = bool(weights[0] >= 0.5)
        assert verdict is sector["comparison"]["retained_subspace_dominance_passed"]
        diagnostic = diagnose_overlap(rebuilt)
        flags = {}
        for direction in ("old_to_new_assignments", "new_to_old_assignments"):
            recorded = sector["comparison"][direction]
            for actual, expected in zip(diagnostic[direction], recorded, strict=True):
                for key in (
                    "assigned_opposite_root",
                    "weak_assignment_trigger",
                    "ambiguous_row_trigger",
                    "retained_to_extra_trigger",
                ):
                    assert actual[key] == expected[key]
            flags[direction] = diagnostic[direction]
        rows.append(
            dict(
                sector=label,
                minimum_retained_squared_principal_overlap=float(weights[0]),
                retained_subspace_dominance_passed=verdict,
                maximum_rebuilt_overlap_error=float(np.max(abs(rebuilt - saved))),
                **flags,
            )
        )
    verdict = all(row["retained_subspace_dominance_passed"] for row in rows)
    assert verdict is report["retained_subspace_dominance_verdict"]
    assert report["original_seed_preservation_passed"] is False
    return dict(
        saved_overlap_reconstruction_verified=True,
        retained_subspace_dominance_verdict=verdict,
        sectors=rows,
        ci_solves=0,
        native_engine_invocations=0,
        original_seed_preservation_passed=False,
        scattering_ensemble_verdict=None,
        scope=(
            "Fresh full-minor contraction from authenticated saved CI/physical metrics; "
            "eigvalsh instead of squared SVD. Same frozen full-minor oracle as the producer "
            "cross-check. CI residual actions and AO integral generation are not repeated."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.root, args.comparison, json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
