"""Finite primitive differential/oracle qualification on archived solve operands."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.numerov_native_residual import native_residual, numpy_residual


def exact(value):
    numerator, denominator = value.as_integer_ratio()
    return mp.mpf(int(numerator)) / int(denominator)


def reference(a, b, x, positions, digits):
    with mp.workdps(digits):
        return [
            exact(b[i, j]) - mp.fsum(exact(a[i, k]) * exact(x[k, j]) for k in range(a.shape[1]))
            for i, j in positions
        ]


def canonical_hash(value):
    """Hash scalar bytes excluding host padding; qualified little-endian x87 only."""
    raw = np.ascontiguousarray(value).view(np.uint8).reshape(-1, value.dtype.itemsize)
    width = 10 if value.dtype == np.dtype(np.longdouble) else value.dtype.itemsize
    return hashlib.sha256(raw[:, :width].tobytes()).hexdigest()


def small_cases():
    dtype = np.longdouble
    a = ((np.arange(77).reshape(7, 11) * 7 % 23) - 11).astype(dtype) / 8
    x = ((np.arange(55).reshape(11, 5) * 11 % 19) - 9).astype(dtype) / 16
    a += dtype(2) ** -60
    b = np.full((7, 5), dtype(3) / 7)
    return {
        "dense_contiguous": (a, b, x),
        "dense_fortran": tuple(np.asfortranarray(v) for v in (a, b, x)),
        "dense_transpose": tuple(v.T.copy().T for v in (a, b, x)),
        "dense_reversed": (a[:, ::-1], b[::-1, ::-1], x[::-1, ::-1]),
        "low_bit_signal": (
            np.ones((1, 1), dtype=dtype) + dtype(2) ** -60,
            np.ones((1, 1), dtype=dtype),
            np.ones((1, 1), dtype=dtype),
        ),
        "large_cancellation": (
            np.array([[dtype(2) ** 60, 1, -(dtype(2) ** 60)]], dtype=dtype),
            np.zeros((1, 1), dtype=dtype),
            np.ones((3, 1), dtype=dtype),
        ),
        "zero": (
            np.zeros((3, 5), dtype=dtype),
            np.zeros((3, 2), dtype=dtype),
            np.zeros((5, 2), dtype=dtype),
        ),
    }


def norm(value):
    return np.max(np.sum(abs(value), axis=1))


def qualify(name, a, b, x, positions, contract, native):
    before = [v.copy() for v in (a, b, x)]
    outputs = {"numpy": numpy_residual(a, b, x)}
    if native:
        outputs["native"] = native_residual(a, b, x)
    denominator = norm(a) * norm(x) + norm(b)
    low, high = [reference(a, b, x, positions, digits) for digits in contract["mp_digits"]]
    eps = np.finfo(a.dtype).eps
    with mp.workdps(contract["mp_digits"][-1]):
        scale = exact(denominator) if denominator else mp.mpf(1)
        pair_error = max(abs(v - w) for v, w in zip(low, high, strict=True)) / scale
        errors = {
            backend: float(
                max(
                    abs(exact(out[i, j]) - want)
                    for (i, j), want in zip(positions, high, strict=True)
                )
                / scale
            )
            for backend, out in outputs.items()
        }
        decimal_reference = [mp.nstr(v, 100) for v in high]
    checks = dict(
        inputs_unchanged=all(
            np.array_equal(v, old) for v, old in zip((a, b, x), before, strict=True)
        ),
        shape_dtype_owned=all(
            out.shape == b.shape
            and out.dtype == a.dtype
            and out.flags.owndata
            and np.isfinite(out).all()
            for out in outputs.values()
        ),
        independent_oracle_pair=pair_error <= contract["oracle_pair_scaled_error_maximum"],
        oracle_accuracy=all(
            error <= contract["scaled_error_eps_factor"] * eps for error in errors.values()
        ),
    )
    difference = (
        float(norm(outputs["native"] - outputs["numpy"]) / denominator)
        if native and denominator
        else 0.0
    )
    if native:
        checks["full_matrix_differential"] = difference <= contract["scaled_error_eps_factor"] * eps
    backwards = {
        backend: float(norm(value) / denominator) if denominator else float(norm(value))
        for backend, value in outputs.items()
    }
    if name.startswith("operands-"):
        checks["accepted_X_backward_gate"] = all(v <= 32 * eps for v in backwards.values())
    return dict(
        name=name,
        dtype=a.dtype.name,
        shapes=[list(v.shape) for v in (a, b, x)],
        reference_positions=positions,
        reference_decimal=decimal_reference,
        oracle_pair_scaled_error=float(pair_error),
        oracle_scaled_errors=errors,
        full_matrix_differential_scaled_error=difference,
        backward_errors=backwards,
        output_numeric_hashes={
            backend: canonical_hash(value) for backend, value in outputs.items()
        },
        checks={key: bool(value) for key, value in checks.items()},
    )


def benchmark(a, b, x, contract):
    functions = dict(numpy=numpy_residual, native=native_residual)
    for function in functions.values():
        for _ in range(contract["timing"]["warmups_per_backend_per_capture"]):
            function(a, b, x)
    rows = []
    for repetition, order in enumerate(contract["timing"]["order"]):
        for backend in order:
            begin = time.perf_counter()
            result = functions[backend](a, b, x)
            elapsed = time.perf_counter() - begin
            rows.append(
                dict(
                    repetition=repetition,
                    backend=backend,
                    wall_seconds=elapsed,
                    output_numeric_hash=canonical_hash(result),
                )
            )
    means = {
        backend: float(np.mean([row["wall_seconds"] for row in rows if row["backend"] == backend]))
        for backend in functions
    }
    speedup = means["numpy"] / means["native"]
    return dict(
        rows=rows,
        mean_seconds=means,
        speedup=speedup,
        performance_hypothesis_passed=speedup
        >= contract["timing"]["minimum_speedup_each_extended_capture"],
    )


def study(root, contract, output):
    import qscat_extended_kernels as kernel

    report = dict(contract=contract, experiment_verdict=None, rows=[], timings=[], checks={})
    begin = time.monotonic()
    try:
        capability = list(kernel.capability())
        report.update(
            native_capability=capability,
            numpy_significand_bits=np.finfo(np.longdouble).nmant + 1,
            extension_sha256=hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest(),
        )
        if capability != contract["native_capability"] or report["numpy_significand_bits"] != 64:
            raise ValueError("Unsupported extended format")
        profile = json.loads((root / "report.json").read_text())
        captures = [item for row in profile["rows"] for item in row["captures"]]
        if len(captures) != contract["captured_triples"]:
            raise ValueError("Capture inventory differs")
        for item in captures:
            path = root / item["file"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
                raise ValueError("Capture bytes differ")
            with np.load(path) as operands:
                a, b, x = [operands[key] for key in ("a", "b", "x")]
            positions = [
                (i, j)
                for i in contract["capture_reference_rows"]
                for j in contract["capture_reference_columns"]
            ]
            native = a.dtype == np.dtype(np.longdouble)
            row = qualify(path.stem, a, b, x, positions, contract, native)
            row["operand_sha256"] = item["sha256"]
            report["rows"].append(row)
            if native:
                report["timings"].append(dict(name=path.stem, **benchmark(a, b, x, contract)))
        cases = small_cases()
        if list(cases) != contract["small_cases"]:
            raise ValueError("Small case inventory differs")
        for name, (a, b, x) in cases.items():
            positions = [(i, j) for i in range(b.shape[0]) for j in range(b.shape[1])]
            row = qualify(name, a, b, x, positions, contract, True)
            if name == "low_bit_signal":
                row["checks"]["signal_not_downcast"] = native_residual(a, b, x)[0, 0] == -(
                    np.longdouble(2) ** -60
                )
                row["checks"]["signal_not_downcast"] = bool(row["checks"]["signal_not_downcast"])
            report["rows"].append(row)
        for row in report["rows"]:
            report["checks"].update(
                {row["name"] + "_" + key: value for key, value in row["checks"].items()}
            )
        report["experiment_verdict"] = all(report["checks"].values())
        report["all_timing_hypotheses_passed"] = all(
            row["performance_hypothesis_passed"] for row in report["timings"]
        )
    except (ValueError, OSError, RuntimeError) as error:
        report["failure"] = f"{type(error).__name__}: {error}"
    finally:
        report["wall_seconds"] = time.monotonic() - begin
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inputs", "contract", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    report = study(args.inputs, json.loads(args.contract.read_text()), args.output)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
