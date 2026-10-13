"""Check an audited engine against the preserved same-root scientific baseline."""

from __future__ import annotations

import argparse
import ctypes
import json
import re
import subprocess
from pathlib import Path

from projects.ukrmol_co.sparse_comparison import compare_runs
from projects.ukrmol_co.sparse_comparison_owner import sha


def label_engine_reproduction(result: dict) -> dict:
    """Same-root agreement qualifies reproduction, not omitted-spectrum convergence."""
    result["matched_engine_baseline_verdict"] = result["root_pair_stability_verdict"]
    result["root_pair_stability_verdict"] = None
    for row in result["sectors"].values():
        row["engine_reproduction_verdict"] = row["root_pair_stability_verdict"]
        row["root_pair_stability_verdict"] = None
    result["scope"] = (
        "Matched-root engine reproduction only; no root-count, omitted-spectrum, "
        "electronic-model or pole qualification"
    )
    return result


def audit_baseline(reference: Path, candidate: Path, contract: dict) -> dict:
    """Require byte-identical native inputs before independently reading both runs."""
    if contract["baseline_roots"] != 2048 or contract["candidate_roots"] != 2048:
        raise ValueError("Matched engine control requires 2048 roots on both sides")
    inputs = audit_native_inputs(reference, candidate, contract, same_roots=True)
    result = label_engine_reproduction(compare_runs(reference, candidate, contract))
    result["matched_native_input_sha256"] = inputs
    return result


def audit_native_inputs(
    reference: Path, candidate: Path, contract: dict, *, same_roots: bool
) -> dict:
    """Hash actual Hamiltonian inputs; only nstat may change for a root comparison."""
    inputs = {}
    for sector in contract["sectors"]:
        left, right = reference / sector, candidate / sector
        decks = [(p / "scatci.inp").read_text() for p in (left, right)]
        if same_roots and (left / "scatci.inp").read_bytes() != (right / "scatci.inp").read_bytes():
            raise ValueError("Matched native deck changed")
        if not same_roots:
            for text, roots in zip(
                decks, (contract["baseline_roots"], contract["candidate_roots"]), strict=True
            ):
                counts = re.findall(r"\bnstat\s*=\s*(\d+)", text)
                if counts != [str(roots)]:
                    raise ValueError("Native requested root count changed")
            if len({re.sub(r"\bnstat\s*=\s*\d+", "nstat = ROOTS", t) for t in decks}) != 1:
                raise ValueError("Native Hamiltonian input changed beyond root count")
        unit = int(re.search(r"megul\s*=\s*(\d+)", (left / "scatci.inp").read_text())[1])
        inputs[sector] = {}
        names = ["fort.16", "fort.24", "fort.26", f"fort.{unit}"]
        if same_roots:
            names.append("scatci.inp")
        for name in names:
            digest = sha(left / name)
            if sha(right / name) != digest:
                raise ValueError(f"Changed matched {sector} input: {name}")
            inputs[sector][name] = digest
    return inputs


def audit_root_pair(reference: Path, candidate: Path, contract: dict) -> dict:
    """Read the declared 2048-to-4096 pair without granting an omission bound."""
    if (contract["baseline_roots"], contract["candidate_roots"]) != (2048, 4096):
        raise ValueError("Audited scientific control requires 2048-to-4096 roots")
    inputs = audit_native_inputs(reference, candidate, contract, same_roots=False)
    result = compare_runs(reference, candidate, contract)
    result["matched_native_input_sha256"] = inputs
    result["engine_correction_verdict"] = None
    return result


def runtime_gate(contract: dict) -> dict:
    """Verify numerical-library bytes and load the retained telemetry in this runtime."""
    hashes = {}
    for name, expected in contract["numerical_library_sha256"].items():
        hashes[name] = sha(Path(name))
        if hashes[name] != expected:
            raise ValueError(f"Changed numerical library: {name}")
    telemetry = Path("/work") / contract["telemetry_directory"] / "scattering_telemetry.so"
    if sha(telemetry) != contract["telemetry_sha256"]:
        raise ValueError("Changed retained telemetry")
    library = ctypes.CDLL(str(telemetry))
    _ = library.EPSSolve
    return {
        "runtime_gate_passed": True,
        "numerical_library_sha256": hashes,
        "telemetry_sha256": sha(telemetry),
        "telemetry_ldd": subprocess.check_output(["ldd", str(telemetry)], text=True),
        "scope": "Dependency byte and symbol-loading gate; native residual audit follows",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("runtime-gate", "compare", "root-compare"))
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--candidate", type=Path)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    if args.action == "runtime-gate":
        result = runtime_gate(contract)
    else:
        if args.reference is None or args.candidate is None:
            parser.error("compare requires --reference and --candidate")
        reader = audit_root_pair if args.action == "root-compare" else audit_baseline
        result = reader(args.reference, args.candidate, contract)
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.action == "compare" and result["matched_engine_baseline_verdict"] is not True:
        raise SystemExit(1)
    if args.action == "root-compare" and result["root_pair_stability_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
