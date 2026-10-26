"""Narrow case/loader adapter and decisions for the saved compressed A2 contrast.

No CLI or automatic propagation: a future owner must persist and pin both real
development checks before the existing pole runner can use a five-model case.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co.compressed_five_view import SELECTED_PROPERTY_IDS, load_compressed_five
from projects.ukrmol_co.native_full_reconstruction import _model_from_interface, load_model
from projects.ukrmol_co.pole_reference import pole_difference
from projects.ukrmol_co.saved_matrices import read_interface


def load_pole_model(source, contract):
    """Dispatch the explicitly qualified finite model views at the compressed anchor."""
    source = Path(source)
    case = contract["case"]
    if "/".join(source.parts[-3:]) != case:
        raise ValueError("Case and source directory differ")
    if case in ["compressed/nine/B1", "compressed/nine/B2"]:
        return load_model(source)
    if (
        case == "compressed/cas12/B1"
        and contract.get("model_view") == "compressed-cas12-native40-v1"
    ):
        interface = read_interface(source / "fort.10", source / "fort.21")
        model = _model_from_interface(source, interface, (40, 250))
        model["interface"] = interface
        return model
    if (
        case in ["compressed/five/B1", "compressed/five/B2"]
        and contract.get("model_view") == "compressed-five-property-ids-v1"
        and tuple(contract.get("selected_property_ids", [])) == SELECTED_PROPERTY_IDS
    ):
        return load_compressed_five(source)
    raise ValueError("An explicit supported compressed model view is required")


def contrast_contract(proposal, base, sector, guard_path, guard_sha256, *, seed=None):
    """Translate the prospective allocation, not authorization to start a job."""
    view_gate = json.loads(
        Path(__file__).with_name("compressed-five-view-gate-contract-20261025.json").read_text()
    )
    proposal_bytes = Path(__file__).with_name(view_gate["source_proposal"]).read_bytes()
    if hashlib.sha256(proposal_bytes).hexdigest() != view_gate[
        "source_proposal_sha256"
    ] or proposal != json.loads(proposal_bytes):
        raise ValueError("The frozen contrast proposal changed")
    if (
        sector not in ["B1", "B2"]
        or proposal["anchor"]["name"] != "compressed"
        or tuple(proposal["five_selected_property_ids_in_written_row_order"])
        != SELECTED_PROPERTY_IDS
    ):
        raise ValueError("Only the declared compressed five-model contrast is supported")
    if sector == "B1":
        z = complex(*proposal["b1_seed_incident_hartree"])
    elif seed is None or not np.isfinite(seed) or complex(seed).imag >= 0:
        raise ValueError("B2 requires the newly qualified five/B1 seed")
    else:
        z = complex(seed)
    result = deepcopy(base)
    half = proposal["b1_box_halfwidth_hartree" if sector == "B1" else "b2_box_halfwidth_hartree"]
    result.update(
        case="compressed/five/" + sector,
        model_view="compressed-five-property-ids-v1",
        selected_property_ids=list(SELECTED_PROPERTY_IDS),
        real_guard_result=guard_path,
        real_guard_budgets={
            key: proposal[key]
            for key in [
                "real_s_absolute_budget",
                "real_invariant_absolute_budget",
                "real_phase_sum_budget_rad",
                "real_flux_absolute_budget",
            ]
        },
        search_box_incident_hartree=[z.real - half, z.real + half, z.imag - half, z.imag + half],
        seeds_incident_hartree=[
            [z.real, z.imag],
            [z.real + proposal["second_seed_delta_hartree"], z.imag],
        ],
        root_seed_refinement_delta_hartree=proposal["second_seed_delta_hartree"],
        settings=proposal["settings"],
        root_evaluation_caps=proposal["root_caps"][sector],
        claim_boundary=(
            "One finite compressed CAS11/five-A2 outgoing pole; no electronic accuracy or curve"
        ),
    )
    for key in [
        "reference_energy_rydberg",
        "same_root_position_full_width_budget_ev",
        "contour_radius_hartree",
        "contour_nodes",
        "projection_columns",
        "projection_seed",
        "derivative_steps_hartree",
        "derivative_residue_relative_budget",
    ]:
        result[key] = proposal[key]
    resources = proposal["future_resources"]
    result["resources"] = {
        key: resources[key]
        for key in [
            "cpuset_cpus",
            "threads_per_worker",
            "contour_workers",
            "memory_gib",
            "per_map_seconds",
            "rhs_cap",
        ]
    }
    result["resources"]["total_map_cap"] = (
        resources["b1_map_cap" if sector == "B1" else "b2_map_cap"] - 2
    )
    result["input_sha256"].update(view_gate["input_sha256"])
    result["input_sha256"][guard_path] = guard_sha256
    return result


def require_real_guard(root, contract):
    """No complex five-model use without two pinned, measured real-axis checks."""
    path = contract.get("real_guard_result")
    expected = contract["input_sha256"].get(path)
    if not path or not expected:
        raise ValueError("A pinned five-model real guard is required")
    payload = (Path(root) / path).read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected:
        raise ValueError("Five-model real guard digest differs")
    data = json.loads(payload)
    if (
        data.get("case") != contract["case"]
        or data.get("scientific_exit") != 0
        or [row.get("index") for row in data.get("real_rows", [])] != [60, 80]
    ):
        raise ValueError("Real guard does not cover the declared five-model development case")
    budgets = contract["real_guard_budgets"]
    for row in data["real_rows"]:
        bounds = [
            (row["s_absolute_error"], budgets["real_s_absolute_budget"]),
            (row["phase_sum_error_mod_pi_rad"], budgets["real_phase_sum_budget_rad"]),
        ]
        for name in ["flux_defect", "wronskian_residual", "unitarity_defect", "reciprocity_defect"]:
            bounds.append(
                (
                    row["invariants"][name],
                    budgets["real_flux_absolute_budget"]
                    if name == "flux_defect"
                    else budgets["real_invariant_absolute_budget"],
                )
            )
        if row.get("completed") is not True or any(
            not np.isfinite(value) or not 0 <= value <= limit for value, limit in bounds
        ):
            raise ValueError("Five-model real-axis prerequisite failed or is incomplete")


def pole_insensitivity(own_pole, reference_pole, own_screen, reference_screen, proposal):
    """Position/FULL-width bands after a necessary, phase-free local identity check."""
    if (
        not np.isfinite([own_pole, reference_pole]).all()
        or own_pole.imag >= 0
        or reference_pole.imag >= 0
    ):
        raise ValueError("Finite outgoing poles with positive full widths required")
    own, reference = own_screen["reference"], reference_screen["reference"]
    if own["canonical_labels"] != reference["canonical_labels"]:
        raise ValueError("Cross-model canonical channel labels differ")
    tv = float(np.sum(abs(np.asarray(own["profile"]) - reference["profile"])) / 2)
    difference = pole_difference(own_pole, reference_pole)
    policy = proposal["model_budgets"]
    tolerance = dict(
        position_ev=policy["position_ev"],
        full_width_ev=max(
            policy["full_width_floor_ev"],
            policy["full_width_relative_to_nine"] * hartree_to_ev(-2 * reference_pole.imag),
        ),
    )
    bands = dict(
        position_ev=policy["position_boundary_band_ev"],
        full_width_ev=policy["full_width_boundary_band_ev"],
    )
    decisions = {
        key: "pass"
        if value < tolerance[key] - bands[key]
        else "reject"
        if value > tolerance[key] + bands[key]
        else "inconclusive"
        for key, value in difference.items()
    }
    identity = (
        own_screen["disposition"] == reference_screen["disposition"] == "pass"
        and tv <= proposal["cross_model_profile_tv_identity_max"]
    )
    disposition = (
        (
            "reject"
            if "reject" in decisions.values()
            else "pass"
            if all(value == "pass" for value in decisions.values())
            else "inconclusive"
        )
        if identity
        else "inconclusive"
    )
    return dict(
        disposition=disposition,
        identity_screen_passed=identity,
        identity_profile_tv=tv,
        difference=difference,
        tolerances=tolerance,
        boundary_bands=bands,
        components=decisions,
    )
