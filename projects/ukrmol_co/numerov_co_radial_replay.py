"""Independent saved radial decisions and explicitly bounded full-row reexecution."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_full_reconstruction import digest, load_model, save, verify_inputs
from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.numerov_native_residual import native_residual
from projects.ukrmol_co.numerov_native_residual_control import canonical_hash
from projects.ukrmol_co.saved_matrices import read_k


def indices(endpoint, h):
    n = round((endpoint - 18) / h)
    result = {0, 2 * n - 2, 2 * n - 1}
    for radius in range(25, 301, 25):
        if radius in [25, 50, 100, 150, 200, 300] and radius < endpoint:
            i = round((radius - 18) / h)
            result.update([2 * i - 1, 2 * i])
    return sorted(result)


def condition_difference(left, right):
    if len(left["samples"]) != len(right["samples"]):
        raise ValueError("Condition sample coverage differs")
    values = []
    for a, b in zip(left["samples"], right["samples"], strict=True):
        x, y = a["reciprocal_condition_1"], b["reciprocal_condition_1"]
        if not np.isfinite([x, y]).all() or min(x, y) <= 0:
            raise ValueError("Invalid sampled reciprocal condition observation")
        values.append(abs(x - y) / max(x, y))
    return max(values)


def scientific_record(row):
    """Exclude measured timing and storage-byte identities, retaining scientific data."""
    result = copy.deepcopy(row)
    for key in ["kernel_wall_seconds", "role_solve_seconds", "unaudited_recurring_bins"]:
        result.pop(key)
    for sample in result["samples"]:
        sample.pop("reciprocal_condition_1")
    for entry in result["captures"] + result["oracles"]:
        entry.pop("sha256")
    return result


def captured_reference(path, sample, stored):
    """Full native/NumPy differential residual plus separately summed mp entries."""
    if digest(path) != stored["sha256"]:
        raise ValueError("Captured file digest differs")
    with np.load(path) as data:
        a, b, x = [data[key] for key in ["a", "b", "x"]]
    width = 6 if sample["role"] == "matching" else 298
    if a.shape != (298, 298) or b.shape != x.shape or b.shape != (298, width):
        raise ValueError("Captured shapes differ")
    for key, value in zip(["a", "b", "x"], [a, b, x], strict=True):
        if value.dtype != np.dtype(np.longdouble) or not np.isfinite(value).all():
            raise ValueError("Captured format or finiteness differs")
        if canonical_hash(value) != sample[key + "_sha256"]:
            raise ValueError("Capture does not match the actual solve")
    residual = b - a @ x
    if not np.array_equal(native_residual(a, b, x), residual):
        raise ValueError("Full captured residual differs from NumPy")
    positions = [[0, 0], [149, width // 2], [297, width - 1]]
    if stored["positions"] != positions:
        raise ValueError("Oracle entry selection differs")

    def exact(value):
        numerator, denominator = value.as_integer_ratio()
        return mp.mpf(numerator) / denominator

    references = []
    for digits in [80, 110]:
        with mp.workdps(digits):
            values = [
                exact(b[i, j]) - mp.fsum(exact(a[i, k]) * exact(x[k, j]) for k in range(298))
                for i, j in positions
            ]
            if [mp.nstr(v, digits) for v in values] != stored["references"][str(digits)]:
                raise ValueError("Independent residual references differ")
            references.append(values)
    with mp.workdps(110):
        norm = lambda value: max(mp.fsum(abs(exact(v)) for v in row) for row in value)  # noqa: E731
        denominator = norm(a) * norm(x) + norm(b)
        stability = max(abs(v - w) for v, w in zip(*references, strict=True)) / denominator
        error = (
            max(
                abs(exact(residual[i, j]) - value)
                for (i, j), value in zip(positions, references[1], strict=True)
            )
            / denominator
        )
        checks = dict(
            precision=stability <= mp.mpf("1e-60"), residual=error <= 32 * mp.mpf(2) ** -63
        )
    if checks != stored["checks"]:
        raise ValueError("Captured reference decisions differ")
    if (
        float(stability) != stored["reference_scaled_stability"]
        or float(error) != stored["residual_scaled_error"]
    ):
        raise ValueError("Captured reference metrics differ")
    return {key: bool(v) for key, v in checks.items()}


def saved_row(row, budget, directory):
    n = round((row["endpoint_bohr"] - 18) / row["step_bohr"])
    samples, stats = row["samples"], row["diagnostics"]
    expected = indices(row["endpoint_bohr"], row["step_bohr"])
    observed = [s["index"] for s in samples]
    if observed != expected:
        raise ValueError("Missing, duplicate or reordered radial sample")
    for sample in samples:
        index = sample["index"]
        role = (
            "initialization"
            if index == 0
            else "matching"
            if index == 2 * n - 1
            else "potential"
            if index % 2
            else "ratio"
        )
        radius = (
            row["endpoint_bohr"]
            if role == "matching"
            else 18 + ((index + 1) // 2) * row["step_bohr"]
        )
        if sample["role"] != role or sample["radius_bohr"] != radius:
            raise ValueError("Radial sample role/location differs")
    matrix = np.asarray(row["k"])
    if matrix.shape != (6, 6) or not np.isfinite(matrix).all():
        raise ValueError("Invalid unsymmetrized saved K")
    s = scattering(matrix)
    decisions = dict(
        solve_count=row["solve_calls"] == stats["large_system_solves"] == 2 * n,
        radial_intervals=stats["radial_intervals"] == n,
        real_and_complex_precision=stats["significand_bits"]
        == stats["complex_significand_bits"]
        == 64,
        backward_error=0 <= stats["maximum_solve_backward_error"] <= 32 * 2.0**-63,
        corrections=0 <= stats["maximum_refinement_corrections"] <= 3,
        sampled_indices=observed == expected,
        condition_observations=all(
            np.isfinite(v["reciprocal_condition_1"]) and v["reciprocal_condition_1"] > 0
            for v in samples
        ),
        unitarity=np.max(abs(s.conj().T @ s - np.eye(6))) <= budget,
        reciprocity=np.max(abs(s - s.T)) <= budget,
    )
    if sum(row["role_solve_seconds"].values()) > row["kernel_wall_seconds"]:
        raise ValueError("Disjoint solve times exceed complete propagation time")
    if sum(v["solves"] for v in row["unaudited_recurring_bins"].values()) != 2 * n - len(samples):
        raise ValueError("Unaudited recurring timing coverage differs")
    if row["backend"] == "native":
        a = row["residual_audit"]
        decisions.update(
            residuals_exact=a["evaluations"] == a["exact_residuals"],
            gate_decisions=a["gate_disagreements"] == 0,
            sampled_solves_audited=a["accepted_solves"] == len(samples),
            independent_backward_error=0 <= a["maximum_accepted_numpy_error"] <= 32 * 2.0**-63,
            precision=a["precision_preserved"],
            capture_count=len(row["captures"]) == 3,
            oracle_count=len(row["oracles"]) == 3,
        )
        if [v["index"] for v in row["captures"]] != [0, 2 * n - 2, 2 * n - 1]:
            raise ValueError("Captured solve coverage differs")
        for i, oracle in enumerate(row["oracles"]):
            capture = next(v for v in row["captures"] if v["file"] == oracle["file"])
            sample = next(v for v in samples if v["index"] == capture["index"])
            if capture["sha256"] != oracle["sha256"]:
                raise ValueError("Oracle/capture digest differs")
            decisions.update(
                {
                    f"oracle_{i}_{k}": v
                    for k, v in captured_reference(
                        directory / oracle["file"], sample, oracle
                    ).items()
                }
            )
    elif row["residual_audit"] is not None or row["captures"] or row["oracles"]:
        raise ValueError("Unexpected reference instrumentation")
    decisions = {key: bool(value) for key, value in decisions.items()}
    if decisions != row["checks"] or all(decisions.values()) is not row["experiment_verdict"]:
        raise ValueError("Independent row decisions differ")
    return decisions


def saved_comparisons(rows, native, old, budget):
    by_id = {row["id"]: row for row in rows}
    metrics, decisions = {}, {}

    def pair(label, left, right):
        metrics[label] = compare(np.asarray(left), np.asarray(right))
        decisions.update(
            {label + "_" + key: value <= budget for key, value in metrics[label].items()}
        )

    for name in ["native-200-fine", "reference-200-fine", "native-400-fine"]:
        row = by_id[name]
        previous = next(
            v
            for v in old["rows"]
            if (v["step_bohr"], v["endpoint_bohr"], v["tail_order"])
            == (0.00625, row["endpoint_bohr"], 8)
        )
        pair(name + "_native64", row["k"], native)
        pair(name + "_increment", row["k"], previous["k"])
    coarse, medium, fine = [by_id["native-200-" + key]["k"] for key in ["coarse", "medium", "fine"]]
    metrics["coarse_step"] = compare(np.asarray(coarse), np.asarray(medium))
    pair("fine_step", medium, fine)
    decisions["s_step_trend"] = metrics["fine_step"]["s"] <= max(
        metrics["coarse_step"]["s"] / 8, 1e-8
    )
    left, right = by_id["native-200-fine"], by_id["reference-200-fine"]
    pair("backend", left["k"], right["k"])
    decisions["extended_k_exact"] = left["extended_k_sha256"] == right["extended_k_sha256"]
    decisions["diagnostics_exact"] = left["diagnostics"] == right["diagnostics"]

    def clean(row):
        return [
            {k: v for k, v in s.items() if k != "reciprocal_condition_1"} for s in row["samples"]
        ]

    decisions["sampled_operands_exact"] = clean(left) == clean(right)
    metrics["auxiliary_rcond_relative_drift"] = condition_difference(left, right)
    decisions["auxiliary_rcond_rounding"] = metrics["auxiliary_rcond_relative_drift"] <= 32 * 2**-52
    pair("endpoint", fine, by_id["native-400-fine"]["k"])
    return {key: bool(v) for key, v in decisions.items()}, metrics


def audit(root, replay_contract):
    import qscat_extended_kernels as native_kernel

    if (
        tuple(native_kernel.capability()) != (16, 64, 16384, True)
        or np.finfo(np.longdouble).nmant != 63
    ):
        raise ValueError("Unqualified scalar capability")
    cpath = root / "source/projects/ukrmol_co/numerov-co-radial-pilot-contract.json"
    if digest(cpath) != replay_contract["producer_contract_sha256"]:
        raise ValueError("Producer contract digest differs")
    c = json.loads(cpath.read_text())
    manifest = root / "collection-sha256.json"
    if digest(manifest) != replay_contract["collection_manifest_sha256"]:
        raise ValueError("Collection manifest digest differs")
    for name, want in json.loads(manifest.read_text())["files"].items():
        if digest(root / "collected" / name) != want:
            raise ValueError("Collected original differs: " + name)
    if digest(root / "collected/result/report.json") != replay_contract["producer_report_sha256"]:
        raise ValueError("Producer report digest differs")
    report = json.loads((root / "collected/result/report.json").read_text())
    if report["contract"] != c:
        raise ValueError("Saved contract differs from frozen input")
    raw_contract = json.loads(
        (root / "source/projects/ukrmol_co/numerov-corner-contract.json").read_text()
    )
    if (
        digest(root / "source/projects/ukrmol_co/numerov-corner-contract.json")
        != c["raw_contract_sha256"]
    ):
        raise ValueError("Raw contract digest differs")
    pins = json.loads((root / "collected/source-manifest.json").read_text())
    for name, want in pins.items():
        if digest(root / "source" / name) != want:
            raise ValueError("Frozen producer source differs: " + name)
    for name, want in c["prior_sha256"].items():
        if digest(root / "collected/prior" / name) != want:
            raise ValueError("Prior evidence digest differs")
    for name, want in c["prerequisite_sha256"].items():
        if digest(root / "source/projects/ukrmol_co" / name) != want:
            raise ValueError("Prerequisite digest differs")
    verify_inputs(root / "raw", raw_contract)
    model = load_model(root / "raw" / c["case"])
    energies, native, header = read_k(root / "collected/prior/native64.bin")
    if (
        float(energies[0]).hex() != c["energy_rydberg_hex"]
        or float(model["energies"][c["energy_index"]]).hex() != c["energy_rydberg_hex"]
    ):
        raise ValueError("Exact archived/native energy differs")
    if native.shape != (1, 6, 6) or any(
        header[k] != v for k, v in model["header"].items() if k != "actual_energy_records"
    ):
        raise ValueError("Reference header/shape differs")
    old = json.loads((root / "collected/prior/corner.json").read_text())
    rows = report["rows"]
    expected = {r["id"]: r for r in c["configurations"]}
    if (
        len(rows) != 5
        or len({r["id"] for r in rows}) != 5
        or set(expected) != {r["id"] for r in rows}
    ):
        raise ValueError("Configuration coverage differs")
    if rows[-1]["id"] != "native-400-fine":
        raise ValueError("Endpoint stage was not last")
    row_decisions = []
    for row in rows:
        if (
            any(row[k] != v for k, v in expected[row["id"]].items())
            or row["case"] != c["case"]
            or row["energy_rydberg_hex"] != c["energy_rydberg_hex"]
        ):
            raise ValueError("Saved row coordinate differs")
        directory = root / "collected/result" / row["id"]
        if json.loads((directory / "report.json").read_text()) != row:
            raise ValueError("Row/aggregate disk records differ")
        row_decisions.append(saved_row(row, c["absolute_budget"], directory))
    if not all(value for decision in row_decisions for value in decision.values()):
        raise ValueError("Saved row qualification did not pass")
    decisions, metrics = saved_comparisons(rows, native[0], old, c["absolute_budget"])
    decisions.update(
        configuration_coverage=True,
        solve_count=sum(r["solve_calls"] for r in rows) == c["mixed_solves"],
    )
    if (
        decisions != report["checks"]
        or report["experiment_verdict"] is not True
        or not all(decisions.values())
    ):
        raise ValueError("Aggregate decisions differ or fail")
    execution = json.loads((root / "collected/runtime-execution.json").read_text())
    owner = json.loads((root / "collected/original-owner-exit.json").read_text())
    if [
        execution["original_command_exit"],
        execution["original_container_exit"],
        owner["exit"],
    ] != [0, 0, 0]:
        raise ValueError("Original exits differ")
    count = len(decisions) + sum(map(len, row_decisions))
    return dict(
        saved_decisions_reconstructed=True,
        experiment_verdict=True,
        rows=5,
        checks_passed=count,
        checks_total=count,
        full_captured_operand_triples=12,
        paired_oracle_entries=36,
        independently_recomputed_full_residuals=12,
        independent_metrics=metrics,
        original_runtime_exits=[0, 0],
        original_owner_exit=0,
        numerical_replay_passed=None,
    )


def reexecute(root, output, replay_contract):
    from projects.ukrmol_co.numerov_co_radial_pilot import worker

    report = json.loads((root / "collected/result/report.json").read_text())
    config = next(
        c for c in report["contract"]["configurations"] if c["id"] == replay_contract["row_id"]
    )
    output.mkdir()
    fresh = worker((root / "raw", output, report["contract"], config))
    old = next(row for row in report["rows"] if row["id"] == config["id"])
    save(output / "fresh-row.json", fresh)
    if fresh.get("experiment_verdict") is not True:
        raise ValueError("Fresh full-row propagation did not pass")
    if scientific_record(old) != scientific_record(fresh):
        raise ValueError("Fresh full-row scientific record differs")
    relative = condition_difference(old, fresh)
    if relative > replay_contract["rcond_relative_budget"]:
        raise ValueError("Fresh auxiliary condition observation exceeds rounding budget")
    for capture in old["captures"]:
        with (
            np.load(root / "collected/result" / old["id"] / capture["file"]) as before,
            np.load(output / fresh["id"] / capture["file"]) as after,
        ):
            if any(not np.array_equal(before[key], after[key]) for key in ["a", "b", "x"]):
                raise ValueError("Fresh full captured operand entries differ")
    return dict(
        numerical_replay_passed=True,
        full_rows_reexecuted=1,
        fresh_mixed_solves=fresh["solve_calls"],
        full_scientific_record_exact=True,
        full_captured_operand_entries_exact=True,
        maximum_rcond_relative_drift=relative,
        fresh_kernel_wall_seconds=fresh["kernel_wall_seconds"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["root", "contract", "output"]:
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--replay-output", type=Path)
    a = parser.parse_args()
    contract = json.loads(a.contract.read_text())
    summary = audit(a.root, contract)
    save(a.output, summary)
    if a.replay_output is not None:
        summary.update(reexecute(a.root, a.replay_output, contract))
        save(a.output, summary)


if __name__ == "__main__":
    main()
