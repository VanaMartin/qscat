"""Read-only deterministic selection from all saved native32 matrices."""

from __future__ import annotations

import json

import numpy as np

from projects.ukrmol_co.native_exterior_attribution import differences
from projects.ukrmol_co.saved_matrices import read_k


def select(saved, history, qualified, cases):
    """Largest interior entrywise K per sector, plus the fixed index48 midpoint.

    Selection observes only saved native matrices, not candidate residuals.
    Previously qualified maxima are retained as reuse decisions, never replaced
    by a more favourable new point. Equal maxima choose the lowest index.
    """
    old = json.loads(history.read_text())
    visited = {
        (x["case"], x["energy_index"])
        for stage in ["probes", "remainder"]
        for x in old[stage]
        if x["rows"]
    }
    prior = json.loads(qualified.read_text())
    if prior["experiment_verdict"] is not True:
        raise ValueError("Previous finite candidate qualification did not pass")
    passed = {(x["case"], x["energy_index"]) for x in prior["cases"] if x["verdict"] is True}
    if len(visited) != 132 or len(passed) != 14:
        raise ValueError("Historical or candidate identity inventory differs")
    peaks, new, reused = [], [], []
    for case in cases:
        unit = "fort.921" if case.endswith("B1") else "fort.922"
        folder = saved / "result/native" / case
        energies, k, header = read_k(folder / "32" / unit)
        e16, k16, h16 = read_k(folder / "16" / unit)
        if k.shape != (99, 6, 6) or not np.array_equal(energies, e16) or header != h16:
            raise ValueError("Saved native reference inventory differs")
        magnitude = np.max(abs(k), axis=(1, 2))
        index = int(np.argmax(magnitude[1:98])) + 1
        eigenvalues = np.linalg.eigvalsh(k[index])
        amplitude = float(np.max(abs(np.eye(6) + k[index] @ k[index])))
        peaks.append(
            dict(
                case=case,
                energy_index=index,
                energy_rydberg_hex=float(energies[index]).hex(),
                requested_energy_ev=float(energies[index] / 0.0735),
                maximum_k_entry=float(magnitude[index]),
                eigenvalues=eigenvalues.tolist(),
                uniform_phase_k_amplification=amplitude,
                leading_dispersion_score=amplitude * float(energies[index]) ** 2.5,
                native_reference_change=differences(k16[index], k[index]),
                already_candidate_qualified=(case, index) in passed,
            )
        )
        for selected, role in [(index, "peak"), (48, "midpoint")]:
            item = dict(
                case=case,
                energy_index=selected,
                role=role,
                energy_rydberg_hex=float(energies[selected]).hex(),
                historically_computed=(case, selected) in visited,
            )
            if (case, selected) in passed:
                reused.append(item)
            elif not any((x["case"], x["energy_index"]) == (case, selected) for x in new):
                new.append(item)
    return dict(
        selection_rule=(
            "Largest max-absolute-entry native32 K on indices1..97, first tie; "
            "plus index48 in every case; retain already-qualified maxima as reuse."
        ),
        native_matrices_examined=594,
        historical_unique_identities=len(visited),
        previous_candidate_qualified_identities=len(passed),
        peaks=peaks,
        new_identities=new,
        reused_qualified_peaks=reused,
        new_unique_identities_if_all_computed=sum(not x["historically_computed"] for x in new),
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
    )
