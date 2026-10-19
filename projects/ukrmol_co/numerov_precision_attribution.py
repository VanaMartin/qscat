"""Fixed-grid free-matrix discrete-oracle and representation attribution."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar, exact_free_k
from projects.ukrmol_co.numerov_interior_analytic import boundary


def precision_difference(a, b):
    with mp.workdps(80):
        return float(max(abs(x - y) for x, y in zip(a, b, strict=True)))


def error_components(k, discrete, continuous):
    """Signed matrix decomposition, with high precision retained until presentation."""
    with mp.workdps(80):
        measured = mp.matrix(np.asarray(k).tolist())
        rounding, truncation = measured - discrete, discrete - continuous
        return dict(
            roundoff_k=float(max(abs(x) for x in rounding)),
            truncation_k=float(max(abs(x) for x in truncation)),
            continuous_k=float(max(abs(x) for x in measured - continuous)),
            roundoff_matrix=np.array(rounding.tolist(), float).tolist(),
            truncation_matrix=np.array(truncation.tolist(), float).tolist(),
            discrete_k_decimal=[[mp.nstr(x, 80) for x in row] for row in discrete.tolist()],
            continuous_k_decimal=[[mp.nstr(x, 80) for x in row] for row in continuous.tolist()],
        )


def oracles(q, thresholds, energy, endpoint, step):
    discrete = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 80)
    high = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 110)
    continuous = exact_free_k(q, thresholds, energy, 18.0, 80)
    refined = exact_free_k(q, thresholds, energy, 18.0, 110)
    return (
        high,
        refined,
        dict(
            discrete_precision=precision_difference(discrete, high),
            continuous_precision=precision_difference(continuous, refined),
        ),
    )


def propagate(q, thresholds, energy, endpoint, step):
    zero = np.zeros_like(q)
    start = time.monotonic()
    k, diagnostics = increment_numerov_k(
        q, thresholds, np.zeros(len(q), int), zero, zero, energy, 18.0, endpoint, step, 8
    )
    return k, dict(
        diagnostics=diagnostics, invariants=algebra(k), wall_seconds=time.monotonic() - start
    )


def study(prior, contract, output):
    report = dict(
        contract=contract,
        diagnostic_verdict=None,
        checks={},
        baseline=[],
        reexecutions=[],
        variants=[],
        hypotheses={},
        original_analytic_verdict=False,
        original_analytic_exits=[1, 1],
        candidate_qualified_identities=14,
        unique_energy_coverage=132,
        uncomputed_identities=462,
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
        native_invocations=0,
        new_inner_solves=0,
        full_sweep_released=False,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
    )
    checks = report["checks"]

    def save():
        (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    start = time.monotonic()
    save()
    try:
        if hashlib.sha256(prior.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
            raise ValueError("Prior report hash differs")
        old = json.loads(prior.read_text())
        if (
            old["analytic_interior_verdict"] is not False
            or sum(old["checks"].values()) != 64
            or len(old["checks"]) != 68
        ):
            raise ValueError("Original rejection differs")
        archived = [x for x in old["rows"] if x["model"] == "compressed"]
        expected = {
            (coupled, h, r, p)
            for coupled in [False, True]
            for h, r, p in [
                (0.00625, 400.0, 8),
                (0.003125, 400.0, 8),
                (0.0015625, 400.0, 8),
                (0.0015625, 200.0, 8),
                (0.0015625, 800.0, 8),
                (0.0015625, 400.0, 6),
            ]
        }
        if (
            len(archived) != 12
            or {
                (x["coupled_closed"], x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
                for x in archived
            }
            != expected
        ):
            raise ValueError("Archived compressed configuration inventory differs")
        energy = float.fromhex(contract["energy_rydberg_hex"])
        for i, row in enumerate(archived):
            if (
                row["energy_rydberg_hex"] != contract["energy_rydberg_hex"]
                or row["target_tangent"] != 320.0
            ):
                raise ValueError("Archived model identity differs")
            q, thresholds = boundary(energy, 320.0, row["coupled_closed"])
            discrete, continuous, precision = oracles(
                q, thresholds, energy, row["endpoint_bohr"], row["step_bohr"]
            )
            checks[f"baseline_{i}_precision"] = all(v < 1e-20 for v in precision.values())
            if not row["coupled_closed"]:
                scalar = exact_discrete_scalar(
                    float(q[0, 0]), energy, 18.0, row["endpoint_bohr"], row["step_bohr"], 110
                )
                with mp.workdps(80):
                    scalar_error = float(abs(discrete[0, 0] - scalar))
                precision["scalar_consistency"] = scalar_error
                checks[f"baseline_{i}_scalar_consistency"] = scalar_error < 1e-20
            report["baseline"].append(
                dict(
                    coupled_closed=row["coupled_closed"],
                    step_bohr=row["step_bohr"],
                    endpoint_bohr=row["endpoint_bohr"],
                    tail_order=row["tail_order"],
                    q=q.tolist(),
                    thresholds=thresholds.tolist(),
                    k=row["k"],
                    precision=precision,
                    **error_components(row["k"], discrete, continuous),
                )
            )
        save()
        if not all(checks.values()):
            report["diagnostic_verdict"] = False
            return report

        # These two reexecutions establish identity with the original binary64 controls.
        for coupled, endpoint in [(False, 800.0), (True, 400.0)]:
            q, thresholds = boundary(energy, 320.0, coupled)
            k, details = propagate(q, thresholds, energy, endpoint, 0.0015625)
            oldrow = next(
                x
                for x in archived
                if x["coupled_closed"] == coupled
                and (x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
                == (0.0015625, endpoint, 8)
            )
            checks[f"reexecution_{coupled}_exact"] = np.array_equal(k, np.asarray(oldrow["k"]))
            checks[f"reexecution_{coupled}_invariants"] = all(
                v <= 1e-6 for v in details["invariants"].values()
            )
            report["reexecutions"].append(
                dict(coupled_closed=coupled, endpoint_bohr=endpoint, k=k.tolist(), **details)
            )
            save()
            if not all(checks.values()):
                report["diagnostic_verdict"] = False
                return report

        q, thresholds = boundary(energy, 320.0, True)
        kap = np.sqrt(thresholds[-1] - energy)
        reduced = q[:2, :2] - np.outer(q[:2, 2], q[2, :2]) * kap / (1 + q[2, 2] * kap)
        eigenvalues, frame = np.linalg.eigh(reduced)
        original = exact_free_k(q, thresholds, energy, 18.0, 110)
        original_float = np.array(original.tolist(), float)
        report["eigenvectors"] = frame.tolist()
        for name, model, transport in [
            ("reduced_rotated", reduced, np.eye(2)),
            ("reduced_diagonal", np.diag(eigenvalues), frame),
        ]:
            for i, (step, endpoint) in enumerate(contract["variant_coordinates"]):
                discrete, continuous, precision = oracles(
                    model, np.zeros(2), energy, endpoint, step
                )
                reference = transport @ np.array(continuous.tolist(), float) @ transport.T
                equivalence = differences(reference, original_float)
                checks[f"{name}_{i}_precision"] = all(v < 1e-20 for v in precision.values())
                checks[f"{name}_{i}_equivalence"] = all(v <= 1e-8 for v in equivalence.values())
                if not all(checks.values()):
                    report["diagnostic_verdict"] = False
                    return report
                k, details = propagate(model, np.zeros(2), energy, endpoint, step)
                errors = error_components(k, discrete, continuous)
                transported_drift = transport @ np.asarray(errors["roundoff_matrix"]) @ transport.T
                checks[f"{name}_{i}_invariants"] = all(
                    v <= 1e-6 for v in details["invariants"].values()
                )
                report["variants"].append(
                    dict(
                        representation=name,
                        step_bohr=step,
                        endpoint_bohr=endpoint,
                        q=model.tolist(),
                        thresholds=[0.0, 0.0],
                        k=k.tolist(),
                        precision=precision,
                        continuous_equivalence=equivalence,
                        transported_k=(transport @ k @ transport.T).tolist(),
                        transported_roundoff_k=float(np.max(abs(transported_drift))),
                        **errors,
                        **details,
                    )
                )
                save()
                if not all(checks.values()):
                    report["diagnostic_verdict"] = False
                    return report
        points = [
            x
            for x in report["baseline"]
            if x["coupled_closed"]
            and x["step_bohr"] == 0.0015625
            and x["endpoint_bohr"] in [200.0, 400.0]
            and x["tail_order"] == 8
        ]
        ratios = [x["roundoff_k"] / x["truncation_k"] for x in points]
        effects = []
        for endpoint in [200.0, 400.0, 800.0]:
            pair = {
                x["representation"]: x["transported_roundoff_k"]
                for x in report["variants"]
                if x["step_bohr"] == 0.0015625 and x["endpoint_bohr"] == endpoint
            }
            effects.append(dict(endpoint_bohr=endpoint, **pair))
        report["hypotheses"] = dict(
            coupled_rejected_points_roundoff_dominates_5fold=all(v >= 5 for v in ratios),
            roundoff_to_truncation_ratios=ratios,
            all_finest_diagonal_roundoff_below_tenth_rotated_and_1em7=all(
                x["reduced_diagonal"] <= min(0.1 * x["reduced_rotated"], 1e-7) for x in effects
            ),
            finest_representation_roundoff=effects,
        )
        report["diagnostic_verdict"] = all(checks.values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
        ZeroDivisionError,
    ) as error:
        report.update(diagnostic_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        report["wall_seconds"] = time.monotonic() - start
        save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["prior", "contract", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    result = study(args.prior, json.loads(args.contract.read_text()), args.output)
    if result["diagnostic_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
