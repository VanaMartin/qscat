"""Explicit property-ID view of the saved compressed250-channel calculation."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_full_reconstruction import _model_from_interface
from projects.ukrmol_co.saved_matrices import read_interface

SELECTED_PROPERTY_IDS = (*range(1, 39), 40, 42)


def selected_ids(source):
    """Require the original pilot and actual OUTER_INTERFACE selection, not a guess."""
    source = Path(source)
    if tuple(source.parts[-3:-1]) != ("compressed", "five") or source.name not in ["B1", "B2"]:
        raise ValueError("Only the saved compressed/five B1/B2 scope is supported")
    pilot = json.loads((source.parent / "pilot.json").read_text())
    deck = source.parents[3] / "writer-evidence/compressed/five" / source.name / "scatci.inp"
    text = "\n".join(line.split("!", 1)[0] for line in deck.read_text().splitlines())
    blocks = re.findall(r"&outer_interface\b(.*?)/", text, re.I | re.S)
    if len(blocks) != 1:
        raise ValueError("One original OUTER_INTERFACE selection required")
    matches = re.findall(r"\bidtarg\s*=\s*([\d,\s]+)", blocks[0], re.I)
    counts = re.findall(r"\bntarg\s*=\s*(\d+)\b", blocks[0], re.I)
    declared = pilot.get("target_map", [])
    actual = [int(value) for value in re.findall(r"\d+", matches[0])] if len(matches) == 1 else []
    if (
        pilot.get("model") != "five"
        or counts != ["40"]
        or declared != actual
        or len(actual) != len(set(actual))
        or tuple(sorted(actual)) != SELECTED_PROPERTY_IDS
    ):
        raise ValueError("Pilot/deck selection differs from the declared compressed property IDs")
    return actual


def load_compressed_five(source):
    """Return the original five-model W/D/CF with a verified dense target lookup."""
    source = Path(source)
    selection = selected_ids(source)
    interface = read_interface(
        source / "fort.10", source / "fort.21", selected_property_ids=selection
    )
    reference = source.parents[1] / "nine" / source.name
    full = read_interface(reference / "fort.10", reference / "fort.21")
    if (full["target_components"], full["boundary_channels"], full["boundary_poles"]) != (
        48,
        298,
        2048,
    ):
        raise ValueError("The original48-target subset reference is missing")
    if (
        interface["target_labels"][0] != [1, 0, 1, 0]
        or interface["channel_labels"][:6] != full["channel_labels"][:6]
        or any(interface[key] != full[key] for key in ["symmetry_charge", "geometry_mass"])
    ):
        raise ValueError("Compressed five/nine ground/open-channel conventions differ")
    ids = np.asarray(interface["selected_property_ids"]) - 1
    expected_energies = np.asarray(full["target_energies_hartree"])[ids]
    error = float(np.max(abs(np.asarray(interface["target_energies_hartree"]) - expected_energies)))
    own_labels = np.asarray(interface["target_labels"])[:, 1:]
    full_labels = np.asarray(full["target_labels"])[ids, 1:]
    if error > 1e-10 or not np.array_equal(own_labels, full_labels):
        raise ValueError("Selected target-table subset differs from the original48-state reference")
    model = _model_from_interface(source, interface, (40, 250))
    model.update(
        interface=interface,
        view_evidence=dict(
            selection_matches_pilot_and_deck=True,
            target_subset_energy_error_hartree=error,
            channels_remapped=sum(
                row[1] != dense
                for row, dense in zip(
                    interface["channel_labels"], interface["channel_target_rows"], strict=True
                )
            ),
            original_native_bytes_edited=False,
            full_reference_used_for_target_metadata_only=True,
        ),
    )
    return model
