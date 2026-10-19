"""Fixed-grid arithmetic control with actual archived binary64 boundary inputs."""

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
from projects.ukrmol_co.numerov_increment_study import exact_free_k
from projects.ukrmol_co.numerov_typed_free import small_solve, typed_free_k


def decimal_matrix(k):
    return [[np.format_float_scientific(x, precision=39, unique=False) for x in row] for row in k]


def oracle_pair(q, thresholds, energy, step, endpoint):
    discrete = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 80)
    refined = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 110)
    continuous = exact_free_k(q, thresholds, energy, 18.0, 80)
    continuous_high = exact_free_k(q, thresholds, energy, 18.0, 110)
    with mp.workdps(80):
        precision = max(float(abs(a - b)) for a, b in zip(discrete, refined, strict=True))
        continuum_precision = max(
            float(abs(a - b)) for a, b in zip(continuous, continuous_high, strict=True)
        )
    return (
        refined,
        continuous_high,
        dict(discrete_precision=precision, continuous_precision=continuum_precision),
    )


def precision_capability():
    bits = {
        name: np.finfo(dtype).nmant + 1
        for name, dtype in [("binary64", np.float64), ("longdouble", np.longdouble)]
    }
    if bits["binary64"] != 53 or bits["longdouble"] < 64:
        raise ValueError(f"Unsupported arithmetic precision: {bits}")
    tiny = np.longdouble(2) ** -60
    a = np.array([[1, tiny], [tiny, 1]], dtype=np.longdouble)
    answer = small_solve(a, np.ones((2, 1), dtype=np.longdouble))
    if (
        answer.dtype != a.dtype
        or answer[0, 0] == np.longdouble(1)
        or abs(answer[0, 0] - (1 - tiny)) >= np.longdouble(2) ** -63
    ):
        raise ValueError("Small solver did not retain extended input precision")
    return bits


def study(prior, contract, output):
    report = dict(
        contract=contract,
        experiment_verdict=None,
        control_checks={},
        qualification_checks={},
        precision_bits={},
        oracle_controls=[],
        original_reexecutions=[],
        rows=[],
        hypotheses={},
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
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
    controls, checks = report["control_checks"], report["qualification_checks"]
    start = time.monotonic()

    def save():
        (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    save()
    try:
        report["precision_bits"] = precision_capability()
        controls["extended_precision_retained"] = True
        if hashlib.sha256(prior.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
            raise ValueError("Prior precision-attribution report hash differs")
        old = json.loads(prior.read_text())
        if old["diagnostic_verdict"] is not True or old["original_analytic_verdict"] is not False:
            raise ValueError("Prior diagnostic or original analytic verdict differs")
        energy = float.fromhex(contract["energy_rydberg_hex"])
        models = []
        for i, coordinate in enumerate(contract["coordinates"]):
            matches = [
                x
                for x in old["baseline"]
                if all(
                    x[k] == coordinate[k] for k in ["coupled_closed", "step_bohr", "endpoint_bohr"]
                )
                and x["tail_order"] == 8
            ]
            if len(matches) != 1:
                raise ValueError("Missing/duplicate archived coordinate")
            row = matches[0]
            q = np.asarray(row["q"])
            thresholds = np.asarray(row["thresholds"])
            discrete, continuous, precision = oracle_pair(
                q, thresholds, energy, coordinate["step_bohr"], coordinate["endpoint_bohr"]
            )
            controls[f"oracle_{i}_precision"] = all(v < 1e-20 for v in precision.values())
            report["oracle_controls"].append(
                dict(
                    coordinate=coordinate,
                    q=q.tolist(),
                    thresholds=thresholds.tolist(),
                    precision=precision,
                )
            )
            models.append((q, thresholds, discrete, continuous))
        save()
        if not all(controls.values()):
            report["experiment_verdict"] = False
            return report
        for coupled, endpoint in [(False, 800.0), (True, 400.0)]:
            source = next(
                x
                for x in old["baseline"]
                if x["coupled_closed"] == coupled
                and (x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
                == (0.0015625, endpoint, 8)
            )
            q = np.asarray(source["q"])
            thresholds = np.asarray(source["thresholds"])
            zero = np.zeros_like(q)
            k, diagnostics = increment_numerov_k(
                q,
                thresholds,
                np.zeros(len(q), int),
                zero,
                zero,
                energy,
                18.0,
                endpoint,
                0.0015625,
                8,
            )
            invariants = algebra(k)
            controls[f"original_{coupled}_exact"] = np.array_equal(k, np.asarray(source["k"]))
            controls[f"original_{coupled}_invariants"] = all(v <= 1e-6 for v in invariants.values())
            report["original_reexecutions"].append(
                dict(
                    coupled_closed=coupled,
                    endpoint_bohr=endpoint,
                    step_bohr=0.0015625,
                    k=k.tolist(),
                    diagnostics=diagnostics,
                    invariants=invariants,
                )
            )
            save()
            if not all(controls.values()):
                report["experiment_verdict"] = False
                return report
        for i, coordinate in enumerate(contract["coordinates"]):
            q, thresholds, discrete, continuous = models[i]
            for backend, dtype in [("binary64", np.float64), ("longdouble", np.longdouble)]:
                begin = time.monotonic()
                k, diagnostics = typed_free_k(
                    q,
                    thresholds,
                    energy,
                    18.0,
                    coordinate["endpoint_bohr"],
                    coordinate["step_bohr"],
                    dtype,
                )
                decimals = decimal_matrix(k)
                with mp.workdps(80):
                    measured = mp.matrix(decimals)
                    roundoff = float(
                        max(abs(a - b) for a, b in zip(measured, discrete, strict=True))
                    )
                    continuous_k = float(
                        max(abs(a - b) for a, b in zip(measured, continuous, strict=True))
                    )
                    discrete_decimal = [[mp.nstr(x, 80) for x in row] for row in discrete.tolist()]
                    continuous_decimal = [
                        [mp.nstr(x, 80) for x in row] for row in continuous.tolist()
                    ]
                k_float = np.asarray(k, float)
                reference = np.array(continuous.tolist(), float)
                errors = differences(k_float, reference)
                invariants = algebra(k_float)
                controls[f"{i}_{backend}_dtype"] = (
                    diagnostics["significand_bits"] == report["precision_bits"][backend]
                )
                controls[f"{i}_{backend}_invariants"] = all(v <= 1e-6 for v in invariants.values())
                if backend == "longdouble":
                    checks[f"{i}_extended_roundoff"] = roundoff <= 1e-8
                    if coordinate["step_bohr"] == 0.0015625:
                        checks[f"{i}_finest_continuous"] = continuous_k <= 1e-6 and all(
                            v <= 1e-6 for key, v in errors.items() if key != "k"
                        )
                report["rows"].append(
                    dict(
                        coordinate_index=i,
                        coordinate=coordinate,
                        backend=backend,
                        k_decimal=decimals,
                        k=k_float.tolist(),
                        roundoff_k=roundoff,
                        continuous_k=continuous_k,
                        continuous_errors=errors,
                        invariants=invariants,
                        discrete_k_decimal=discrete_decimal,
                        continuous_k_decimal=continuous_decimal,
                        diagnostics=diagnostics,
                        wall_seconds=time.monotonic() - begin,
                    )
                )
                save()
                if not all(controls.values()):
                    report["experiment_verdict"] = False
                    return report
        ratios = []
        for i in range(len(contract["coordinates"])):
            pair = {
                r["backend"]: r["roundoff_k"] for r in report["rows"] if r["coordinate_index"] == i
            }
            ratios.append(pair["binary64"] / pair["longdouble"])
        report["hypotheses"] = dict(
            all_six_roundoff_reduction_at_least_100fold=all(x >= 100 for x in ratios),
            roundoff_reduction_ratios=ratios,
        )
        report["experiment_verdict"] = all(controls.values()) and all(checks.values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
        ZeroDivisionError,
    ) as error:
        report.update(experiment_verdict=None, failure=f"{type(error).__name__}: {error}")
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
    if result["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
