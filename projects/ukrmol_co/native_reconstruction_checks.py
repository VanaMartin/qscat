"""Frozen real-axis reconstruction configurations and per-energy qualification gates."""

from __future__ import annotations

from itertools import pairwise

import numpy as np

from projects.ukrmol_co.native_exterior_attribution import algebra, differences


def configurations(contract, probe):
    h = contract["steps_bohr"]
    radius = contract["comparison_endpoint_bohr"]
    order = contract["tail_orders"][-1]
    if not probe:
        return [(step, radius, order) for step in h[-2:]]
    return list(
        dict.fromkeys(
            [(step, radius, order) for step in h]
            + [(h[-1], r, order) for r in contract["endpoints_bohr"]]
            + [(h[-1], radius, p) for p in contract["tail_orders"]]
        )
    )


def energy_checks(rows, reference, contract, probe):
    """Keep K, S, determinant phase and physical invariants as separate vetoes."""
    expected = configurations(contract, probe)
    keys = [(r["step_bohr"], r["endpoint_bohr"], r["tail_order"]) for r in rows]
    if keys != expected or not all(r["completed"] for r in rows):
        raise ValueError("Incomplete, duplicated or reordered independent configuration family")
    matrices = {key: np.asarray(row["k"]) for key, row in zip(keys, rows, strict=True)}
    if any(k.shape != reference.shape or not np.isfinite(k).all() for k in matrices.values()):
        raise ValueError("Invalid independent matrix")
    h = contract["steps_bohr"]
    radius = contract["comparison_endpoint_bohr"]
    order = contract["tail_orders"][-1]
    finest = matrices[h[-1], radius, order]
    metrics = {"native_comparison": differences(finest, reference)}
    families = {"step": [(step, radius, order) for step in h[-(3 if probe else 2) :]]}
    if probe:
        families.update(
            endpoint=[(h[-1], r, order) for r in contract["endpoints_bohr"]],
            order=[(h[-1], radius, p) for p in contract["tail_orders"]],
        )
    checks = {}
    budget = contract["matrix_and_phase_absolute_budget"]
    for family, sequence in families.items():
        values = [differences(matrices[a], matrices[b]) for a, b in pairwise(sequence)]
        metrics[family + "_changes"] = values
        for name, value in values[-1].items():
            checks[family + "_" + name] = bool(np.isfinite(value) and value <= budget)
    for name, value in metrics["native_comparison"].items():
        checks["native_" + name] = bool(np.isfinite(value) and value <= budget)
    if probe:
        coarse, fine = metrics["step_changes"]
        checks["spatial_trend"] = fine["s"] <= max(
            coarse["s"] / contract["spatial_improvement_factor"], contract["roundoff_floor"]
        )
    invariants = [algebra(k) for k in matrices.values()]
    metrics["invariant_maxima"] = {
        name: max(row[name] for row in invariants) for name in invariants[0]
    }
    for name, value in metrics["invariant_maxima"].items():
        checks[name] = bool(np.isfinite(value) and value <= contract["invariant_absolute_budget"])
    return checks, metrics


def complete_coverage(report, contract):
    """A matching count alone cannot substitute for unique case/energy identities."""
    expected = {
        (case, index)
        for case in contract["cases"]
        for index in range(contract["expected_energy_points"])
    }
    rows = report["probes"] + report["remainder"]
    identities = [(r["case"], r["energy_index"]) for r in rows]
    return (
        len(rows) == contract["full_energy_identities"]
        and len(set(identities)) == len(rows)
        and set(identities) == expected
        and all(r["verdict"] is True for r in rows)
        and sum(len(r["rows"]) for r in rows) == contract["independent_solves_if_complete"]
    )
