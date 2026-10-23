"""Phase-free support diagnostics of finite-model open-channel pole residues."""

from __future__ import annotations

import numpy as np

from projects.ukrmol_co.pole_reference import pole_difference


def canonical_labels(labels):
    """Validate the six ground-target channels; discard ordinals and the m sign."""
    labels = np.asarray(labels)
    if labels.shape != (6, 4) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("Six native integer (channel,target,l,m) labels required")
    keys = [
        (int(target), int(angular), abs(int(projection)))
        for _, target, angular, projection in labels
    ]
    expected = {(1, 1, 1), (1, 2, 1), (1, 3, 1), (1, 3, 3), (1, 4, 1), (1, 4, 3)}
    if set(labels[:, 0]) != set(range(1, 7)) or set(keys) != expected:
        raise ValueError("Unique verified ground-target/l/m channel labels required")
    return keys


def _band(value, passed, rejected):
    if value <= passed:
        return "pass"
    if value >= rejected:
        return "reject"
    return "inconclusive"


def _combine(dispositions):
    if "reject" in dispositions:
        return "reject"
    return "pass" if all(value == "pass" for value in dispositions) else "inconclusive"


def character(residue, labels, full_width, policy):
    """Support of an S residue, not a probability, branching ratio or Lambda proof."""
    keys = canonical_labels(labels)
    residue = np.asarray(residue, complex)
    if residue.shape != (6, 6) or not np.isfinite(residue).all():
        raise ValueError("A finite6x6 open-channel S residue is required")
    if not np.isfinite(full_width) or full_width <= 0:
        raise ValueError("Finite positive full width in the residue's energy unit required")
    scale = float(np.max(abs(residue)))
    scaled = residue / scale if scale else residue
    squared = abs(scaled) ** 2
    total = float(squared.sum())
    norm = float(np.sqrt(total) * scale)
    if not np.isfinite(norm):
        raise ValueError("Residue norm overflows")
    singular = np.linalg.svd(scaled, compute_uv=False)
    rank_ratio = float(singular[1] / singular[0]) if scale else None
    transpose = float(np.linalg.norm(scaled.T - scaled) / np.sqrt(total)) if scale else None
    guards = dict(
        signal=norm / full_width > policy["signal_over_full_width_min"],
        rank=rank_ratio is not None and rank_ratio <= policy["rank_ratio_max"],
        transpose=transpose is not None and transpose <= policy["transpose_relative_max"],
    )
    complement = np.array([key[2] != 1 for key in keys])
    outgoing = float(squared[complement, :].sum() / total) if scale else None
    incoming = float(squared[:, complement].sum() / total) if scale else None
    leakage = max(incoming, outgoing) if scale else None
    profile = (squared.sum(axis=0) + squared.sum(axis=1)) / (2 * total) if scale else np.zeros(6)
    order = sorted(range(6), key=lambda i: keys[i])
    disposition = (
        _band(leakage, policy["leakage_pass"], policy["leakage_reject"])
        if all(guards.values())
        else "inconclusive"
    )
    return dict(
        disposition=disposition,
        guards=guards,
        residue_norm=norm,
        rank_ratio=rank_ratio,
        transpose_relative_error=transpose,
        incoming_leakage=incoming,
        outgoing_leakage=outgoing,
        leakage=leakage,
        canonical_labels=[list(keys[i]) for i in order],
        profile=[float(profile[i]) for i in order],
    )


def screen_witnesses(residues, labels, full_width, policy):
    """Require agreement of16/32-node and two derivative residue estimates."""
    if len(residues) != 4:
        raise ValueError("Exactly four frozen residue witnesses required")
    rows = [character(value, labels, full_width, policy) for value in residues]
    reference = rows[1]  # 32-node residue, never a selected best-scoring estimate.
    if not all(all(row["guards"].values()) for row in rows):
        return dict(disposition="inconclusive", reference=reference, witnesses=rows, stability=None)
    leakage_change = max(abs(row["leakage"] - reference["leakage"]) for row in rows)
    profile_change = max(
        float(np.sum(abs(np.array(row["profile"]) - reference["profile"])) / 2) for row in rows
    )
    norm_change = max(
        abs(row["residue_norm"] - reference["residue_norm"])
        / max(row["residue_norm"], reference["residue_norm"])
        for row in rows
    )
    checks = dict(
        dispositions=len({row["disposition"] for row in rows}) == 1,
        leakage=leakage_change <= policy["witness_leakage_change_max"],
        profile=profile_change <= policy["witness_profile_tv_max"],
        norm=norm_change <= policy["witness_norm_relative_max"],
    )
    return dict(
        disposition=reference["disposition"] if all(checks.values()) else "inconclusive",
        reference=reference,
        witnesses=rows,
        stability=dict(
            leakage_change=leakage_change,
            profile_tv_change=profile_change,
            norm_relative_change=norm_change,
            checks=checks,
        ),
    )


def compare_pair(left, right, left_pole, right_pole, policy, energy_budget_ev):
    """Necessary counterpart checks; phase-free profiles are not full residues/S."""
    a, b = left["reference"], right["reference"]
    if a["canonical_labels"] != b["canonical_labels"]:
        raise ValueError("Canonical counterpart labels differ")
    if not np.isfinite([left_pole, right_pole]).all() or energy_budget_ev <= 0:
        raise ValueError("Finite poles and a positive energy budget required")
    errors = pole_difference(left_pole, right_pole)
    tv = float(np.sum(abs(np.array(a["profile"]) - b["profile"])) / 2)
    largest = max(a["residue_norm"], b["residue_norm"])
    strength = abs(a["residue_norm"] - b["residue_norm"]) / largest if largest else None
    dispositions = dict(
        b1=left["disposition"],
        b2=right["disposition"],
        energy="pass" if max(errors.values()) <= energy_budget_ev else "reject",
        profile=_band(tv, policy["pair_profile_tv_pass"], policy["pair_profile_tv_reject"]),
        norm=(
            _band(strength, policy["pair_norm_relative_pass"], policy["pair_norm_relative_reject"])
            if strength is not None
            else "inconclusive"
        ),
    )
    # A deficient signal/rank cannot be classified by its noisy profile or norm.
    valid = all(
        all(row["guards"].values()) for group in [left, right] for row in group["witnesses"]
    )
    disposition = _combine(list(dispositions.values())) if valid else "inconclusive"
    return dict(
        disposition=disposition,
        components=dispositions,
        pole_difference=errors,
        profile_tv=tv,
        norm_relative_difference=strength,
    )
