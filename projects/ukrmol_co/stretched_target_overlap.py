"""One bounded saved-native DZ/TZ target-overlap screen; no target eigensolver."""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps
from projects.ukrmol_co.ensemble_frame_overlaps import full_minor_overlaps
from projects.ukrmol_co.native_target_determinants import cidata, expand
from projects.ukrmol_co.stretched_basis_properties import digest, save

CONTRACT = "projects/ukrmol_co/stretched-target-overlap-contract-20261102.json"
UNITS = [710, 730, 711, 731, 712, 732, 713, 733]


def inputs(root, source, c):
    for name, pin in c["source_sha256"].items():
        path = source / name
        if not path.resolve().is_relative_to(source.resolve()) or digest(path) != pin:
            raise ValueError("Changed source prerequisite " + name)
    collection = root / c["feasibility_collection"]
    if (
        not collection.resolve().is_relative_to(root.resolve())
        or digest(collection) != c["feasibility_collection_sha256"]
    ):
        raise ValueError("Original public feasibility collection")
    original = json.loads(collection.read_text())["files"]
    models = json.loads((source / c["property_contract"]).read_text())["models"]
    pins = {c["feasibility_collection"]: c["feasibility_collection_sha256"]}
    for model in models.values():
        for name, pin in model["input_sha256"].items():
            if name == "target.json" or (
                name.startswith("output/CO/geom1/")
                and Path(name).name
                in ["co.casscf.chk", "co.molden", "fort.26", *[f"fort.{u}" for u in UNITS]]
            ):
                pins[model["target_directory"] + "/" + name] = pin
    for spin in [0, 2]:
        for irrep in ["A1", "B1", "B2", "A2"]:
            for filename in ["result.json", "ci.npy"]:
                name = c["tz_reference_directory"] + f"/spin{spin}-{irrep}-space240/" + filename
                pins[name] = original[name]["sha256"]
    for name, pin in pins.items():
        path = root / name
        if not path.resolve().is_relative_to(root.resolve()) or digest(path) != pin:
            raise ValueError("Changed saved input " + name)
    if len(pins) != 41 or not c["execution_authorized"]:
        raise ValueError("Exact authorized41-input inventory")
    return models, pins


def canonical_order(target):
    labels = [(r["irrep"], r["index_in_irrep"]) for r in target["orbitals"][:14]]
    wanted = [
        (irrep, i)
        for irrep, count in [("A1", 8), ("B1", 3), ("B2", 3)]
        for i in range(1, count + 1)
    ]
    if len(set(labels)) != 14 or set(labels) != set(wanted) or labels[:2] != [("A1", 1), ("A1", 2)]:
        raise ValueError("Actual14-column orbital inventory")
    return [labels.index(k) for k in wanted]


def metric(checkpoints, targets, output):
    from pyscf import gto, lib, symm

    molecules = []
    coefficients = []
    own_metrics = []
    for case, path, target in zip(["dz", "tz"], checkpoints, targets, strict=True):
        mol = lib.chkfile.load_mol(str(path))
        full = np.asarray(lib.chkfile.load(str(path), "mcscf/mo_coeff"))
        if (
            int(lib.chkfile.load(str(path), "mcscf/ncore")) != 2
            or int(lib.chkfile.load(str(path), "mcscf/ncas")) != 12
        ):
            raise ValueError("Checkpoint core/active inventory")
        order = canonical_order(target)
        c = full[:, order]
        ao = mol.intor_symmetric("int1e_ovlp")
        own = c.conj().T @ ao @ c
        np.testing.assert_allclose(own, np.eye(14), atol=1e-9, rtol=0)
        actual = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, c)
        if list(actual) != ["A1"] * 8 + ["B1"] * 3 + ["B2"] * 3:
            raise ValueError("Checkpoint physical orbital symmetries")
        np.savez(output / f"{case}-metric.npz", mo=c, ao=ao, own=own, original_order=order)
        molecules.append(mol)
        coefficients.append(c)
        own_metrics.append(own)
    a, b = molecules
    np.testing.assert_array_equal(a.atom_charges(), b.atom_charges())
    np.testing.assert_allclose(a.atom_coords(), b.atom_coords(), atol=1e-12, rtol=0)
    np.testing.assert_allclose(
        np.linalg.norm(a.atom_coords()[1] - a.atom_coords()[0]), 2.5, atol=1e-12, rtol=0
    )
    forward = gto.intor_cross("int1e_ovlp", a, b)
    reverse = gto.intor_cross("int1e_ovlp", b, a)
    np.testing.assert_allclose(forward, reverse.conj().T, atol=1e-12, rtol=0)
    cross = coefficients[0].conj().T @ forward @ coefficients[1]
    back = coefficients[1].conj().T @ reverse @ coefficients[0]
    np.testing.assert_allclose(cross, back.conj().T, atol=1e-12, rtol=0)
    if np.linalg.svd(cross, compute_uv=False)[0] > 1 + 1e-9:
        raise ValueError("Noncontractive physical MO overlap")
    np.savez(output / "cross-metric.npz", ao=forward, reverse_ao=reverse, mo=cross, reverse_mo=back)
    return cross, back, own_metrics


def diagnose(overlap):
    singular = np.linalg.svd(overlap, compute_uv=False)
    if overlap.shape != (8, 8) or not np.isfinite(overlap).all() or singular[0] > 1 + 1e-9:
        raise ValueError("Physical eight-root overlap")
    squared = abs(overlap) ** 2
    rows, columns = linear_sum_assignment(-squared)
    return dict(
        singular_values=singular.tolist(),
        squared_principal_overlaps=(singular**2).tolist(),
        minimum_squared_principal_overlap=float(singular[-1] ** 2),
        retained_eight_dominance_passed=bool(singular[-1] ** 2 >= 0.5),
        dz_capture_in_tz=squared.sum(axis=1).tolist(),
        tz_capture_in_dz=squared.sum(axis=0).tolist(),
        assignment=[
            dict(dz_root=int(i + 1), tz_root=int(j + 1), squared_overlap=float(squared[i, j]))
            for i, j in zip(rows, columns, strict=True)
        ],
    )


def run(root, source, output, c):
    start = time.monotonic()
    report = dict(
        scientific_exit=2,
        retained_eight_dominance_verdict=None,
        sectors=[],
        native_calls=0,
        new_ci_solves=0,
        new_hamiltonian_actions=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
    )

    def persist():
        report["wall_seconds"] = time.monotonic() - start
        save(output / "result.json", report)

    try:
        models, pins = inputs(root, source, c)
        save(output / "input-manifest.json", pins)
        targets = [
            json.loads((root / models[k]["target_directory"] / "target.json").read_text())
            for k in ["dz", "tz"]
        ]
        bases = [root / models[k]["target_directory"] / "output/CO/geom1" for k in ["dz", "tz"]]
        natives = [cidata(base / "fort.26") for base in bases]
        cross, back, own = metric([b / "co.casscf.chk" for b in bases], targets, output)
        for index, (irrep, spin) in enumerate(
            (i, s) for i in ["A1", "B1", "B2", "A2"] for s in [0, 1]
        ):
            if time.monotonic() - start > c["science_cap_seconds"]:
                raise TimeoutError("Science cap")
            label = f"{'singlet' if spin == 0 else 'triplet'}.{irrep}"
            vectors = []
            details = []
            for base, native in zip(bases, natives, strict=True):
                vector, detail = expand(
                    base / f"fort.{UNITS[index]}", native[(spin, irrep)], spin, irrep
                )
                vectors.append(vector)
                details.append(detail)
            electrons = (5 + spin, 5 - spin)
            reference = root / c["tz_reference_directory"] / f"spin{2 * spin}-{irrep}-space240"
            saved = json.loads((reference / "result.json").read_text())
            ci = np.load(reference / "ci.npy", allow_pickle=False)
            if not (
                saved["physical_gates_pass"]
                and saved["ci_converged"] == [True] * 8
                and saved["ci_orthogonality_error"] <= 1e-9
            ):
                raise ValueError("Saved TZ CI qualification")
            for state in saved["states"]:
                if (
                    state["spin_square_error"] > 1e-6
                    or max(
                        state[k]
                        for k in [
                            "physical_residual_hartree",
                            "penalized_residual_hartree",
                            "spin_penalty_action_hartree",
                        ]
                    )
                    > 1e-9
                ):
                    raise ValueError("Saved TZ action/spin gate")
            np.testing.assert_allclose(
                [s["energy_hartree"] for s in saved["states"]],
                natives[1][(spin, irrep)]["energy"],
                atol=1e-7,
                rtol=0,
            )
            order = canonical_order(targets[1])
            permutation = np.eye(14)[order]
            conversion = frozen_core_overlaps(vectors[1], ci, permutation, 2, electrons)
            if min(abs(np.diag(conversion)) ** 2) < 1 - 1e-8:
                raise ValueError("Native determinant conversion differs from saved TZ CI")
            left, right = vectors
            overlap = frozen_core_overlaps(left, right, cross, 2, electrons)
            explicit = full_minor_overlaps(left, right, cross, 2, electrons)
            reverse = frozen_core_overlaps(right, left, back, 2, electrons)
            np.testing.assert_allclose(overlap, explicit, atol=1e-10, rtol=0)
            np.testing.assert_allclose(overlap, reverse.conj().T, atol=1e-10, rtol=0)
            for v, m in zip(vectors, own, strict=True):
                np.testing.assert_allclose(
                    frozen_core_overlaps(v, v, m, 2, electrons), np.eye(8), atol=1e-9, rtol=0
                )
            np.savez(
                output / (label + ".npz"),
                dz=left,
                tz=right,
                overlap=overlap,
                full_minors=explicit,
                reverse=reverse,
                conversion_overlap=conversion,
            )
            row = dict(
                sector=label,
                expansions=details,
                conversion_min_squared_overlap=float(min(abs(np.diag(conversion)) ** 2)),
                full_minor_error=float(abs(overlap - explicit).max()),
                reciprocity_error=float(abs(overlap - reverse.conj().T).max()),
                **diagnose(overlap),
            )
            report["sectors"].append(row)
            persist()
        verdict = all(s["retained_eight_dominance_passed"] for s in report["sectors"])
        report.update(
            retained_eight_dominance_verdict=verdict,
            scientific_exit=0 if verdict else 1,
            cross_basis_target_overlap_numerical_verdict=True,
            global_resonance_identity_verdict=None,
            scope=(
                "Finite retained-eight physical target-subspace dominance,"
                "not an N+1 pole identity or production-error bound"
            ),
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            inputs(root, source, c)
            report["inputs_unchanged"] = True
        except Exception as error:
            report.update(integrity_error=repr(error), scientific_exit=2)
        persist()
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    a = parser.parse_args()
    c = json.loads(a.contract.read_text())
    if a.profiled:

        def timeout(*_):
            raise TimeoutError("Authorized target-overlap science cap")

        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(c["science_cap_seconds"])
        return run(a.root, a.source, a.output, c)
    from projects.ukrmol_co.run import run_profiled

    a.output.mkdir(exist_ok=False)
    (a.output / "scratch").mkdir()
    return run_profiled(
        a.output,
        os.environ | dict(TMPDIR=str(a.output / "scratch")),
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.stretched_target_overlap",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
