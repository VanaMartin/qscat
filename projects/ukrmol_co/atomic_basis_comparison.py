"""Finite QZ/5Z endpoint diagnostics on fixed neutral references,not a De/CBS claim."""

from __future__ import annotations

import math

from qscat.units import HARTREE_TO_EV


def compare(qz, five_z, neutral, budgets):
    """Use qualified saved atoms and matching-basis molecular reference energies."""
    energies = {}
    for basis in ("aug-cc-pVQZ", "aug-cc-pV5Z"):
        candidates = [
            r
            for r in (run.get("result") or {} for run in neutral["runs"].values())
            if r.get("basis") == basis and r.get("bond_length_bohr") == 2.1323
        ]
        if len(candidates) != 1:
            raise ValueError("One unique fixed molecular reference required")
        r = candidates[0]
        if not (
            r["frozen_spatial_orbitals"] == 2
            and r["rhf_internally_stable"]
            and r["rhf_externally_stable"]
            and r["ccsd_converged"]
            and r["lambda_converged"]
            and math.isfinite(r["neutral_energy_hartree"])
        ):
            raise ValueError("Unqualified fixed molecular reference")
        energies[basis] = r["neutral_energy_hartree"]
    offsets = {}
    for label, result, basis in [("QZ", qz, "aug-cc-pVQZ"), ("5Z", five_z, "aug-cc-pV5Z")]:
        if not result["finite_method_endpoint_reconstruction_verdict"]:
            raise ValueError("Both finite atomic results must independently reconstruct")
        offsets[label] = {
            k: (v - energies[basis]) * HARTREE_TO_EV
            for k, v in result["fragment_sum_hartree"].items()
        }
    ea_changes = {
        k: five_z["electron_affinity_ev"][k] - qz["electron_affinity_ev"][k] for k in ("O", "C")
    }
    offset_changes = {k: offsets["5Z"][k] - offsets["QZ"][k] for k in offsets["QZ"]}
    if not all(math.isfinite(v) for v in [*ea_changes.values(), *offset_changes.values()]):
        raise ValueError("Finite endpoint differences required")
    ea_pass = max(map(abs, ea_changes.values())) <= budgets["ea_basis_change_ev"]
    offset_pass = (
        max(map(abs, offset_changes.values()))
        <= budgets["atomic_to_fixed_molecular_reference_offset_basis_change_ev"]
    )
    return dict(
        ea_basis_change_ev=ea_changes,
        atomic_to_fixed_reference_offset_ev=offsets,
        atomic_to_fixed_reference_offset_basis_change_ev=offset_changes,
        ea_basis_diagnostic_passed=ea_pass,
        reference_offset_basis_diagnostic_passed=offset_pass,
        finite_basis_diagnostic_passed=ea_pass and offset_pass,
        molecular_reference_energy_hartree=energies,
        molecular_reference_R_bohr=2.1323,
        reference_offset_is_dissociation_depth=False,
        complete_basis_convergence_qualified=False,
        atomic_absolute_accuracy_qualified=False,
        molecular_dissociation_qualified=False,
        no_automatic_larger_basis=True,
        scope="Two finite bases under the same atomic method;"
        "molecular offsets use saved restricted-neutral references,"
        "not a demonstrated dissociation path",
    )
