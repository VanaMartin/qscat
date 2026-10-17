"""Render the dated CO qualification overview from its small evidence snapshot.

Run with ``python -m projects.ukrmol_co.progress_figure`` and matplotlib installed.
Literature comparisons and diagnostic gate ratios retain distinct meanings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from io import StringIO
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SNAPSHOT = Path(__file__).with_name("progress-snapshot-20261017T1340.json")
BLUE = "#2468ac"
GREEN = "#158477"
ORANGE = "#cf7b25"
RED = "#b64b50"
INK = "#233649"
GRAY = "#73808b"


def _panel(ax, letter, title, subtitle):
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=14, pad=38)
    ax.text(0, 1.035, subtitle, transform=ax.transAxes, fontsize=9.5, color=GRAY)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e3e8ed", linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)


def _note(ax, text):
    ax.text(0, -0.27, text, transform=ax.transAxes, fontsize=9.2, color=INK, va="top")


def _resonances(ax, data):
    _panel(
        ax,
        "A",
        "Where do our resonances sit?",
        "Equilibrium R = 2.1323 bohr; native fit candidates",
    )
    reference_labels = [
        ("Laporta 2012\nSEP", (1.20, 0.91)),
        ("Dora 2016\nTZ / CAS10 / 50", (1.70, 0.63)),
        ("Dora 2016\nDZ / CAS10 / 40", (1.90, 0.52)),
        ("Dora 2016\nDZ / CAS11 / 40", (2.30, 0.83)),
        ("Dora 2020\n6Z / 41", (1.49, 1.37)),
    ]
    for row, (label, position) in zip(data["literature"], reference_labels, strict=True):
        xy = row["position_ev"], row["full_width_ev"]
        ax.scatter(*xy, marker="*", s=150, c=INK, zorder=5)
        ax.annotate(
            label,
            xy,
            xytext=position,
            fontsize=8.5,
            color=INK,
            arrowprops={"arrowstyle": "-", "color": GRAY, "lw": 0.7},
        )
    our_labels = [
        ("Our SA-CAS10", (2.30, 1.02)),
        ("Our SA-CAS11", (2.35, 1.27)),
        ("Our SEP39", (1.20, 0.73)),
        ("Our SEP60", (1.16, 0.47)),
    ]
    for i, (row, (label, position)) in enumerate(
        zip(data["equilibrium_candidates"], our_labels, strict=True)
    ):
        color = BLUE if i < 2 else ORANGE
        xy = row["position_ev"], row["full_width_ev"]
        ax.scatter(*xy, marker="o" if i < 2 else "s", s=65, c=color, zorder=6)
        ax.annotate(
            label,
            xy,
            xytext=position,
            fontsize=9,
            fontweight="bold",
            color=color,
            arrowprops={"arrowstyle": "-", "color": color, "lw": 0.8},
        )
    first, second = data["equilibrium_candidates"][2:]
    ax.annotate(
        "",
        (second["position_ev"], second["full_width_ev"]),
        (first["position_ev"], first["full_width_ev"]),
        arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 1.6, "linestyle": "--"},
    )
    ax.set(
        xlim=(1.12, 2.82), ylim=(0.4, 1.53), xlabel="Position relative to neutral threshold / eV"
    )
    ax.set_ylabel("Full width / eV")
    _note(
        ax,
        "Stars: published model values, not experimental error bars.\n"
        "SEP39 → SEP60 drifts by 0.294 eV despite nearby literature.\n"
        "Next: balance basis, active space, ensemble and channels.",
    )


def _diagnostics(ax, data):
    _panel(
        ax,
        "B",
        "What still limits reliability?",
        "Each sensitivity divided by its own declared tolerance",
    )
    rows = data["diagnostic_gate_ratios"]
    values = np.asarray([row["value"] for row in rows])
    y = np.arange(len(rows))
    ax.axvspan(1e-7, 1, color=GREEN, alpha=0.05)
    ax.axvspan(1, 30, color=RED, alpha=0.05)
    ax.barh(
        y,
        values - 1e-7,
        left=1e-7,
        color=[GREEN if row.get("passed", i < 3) else RED for i, row in enumerate(rows)],
        height=0.56,
    )
    ax.axvline(1, color=INK, lw=1.1, linestyle="--")
    labels = data.get(
        "diagnostic_labels",
        [
            "Sparse vs dense\nphase",
            "Propagation\nphase",
            "Angular cutoff\nphase",
            "Active space\nphase",
            "Stretched\nfit width",
            "Equilibrium\nfit width",
            "Compressed\nfit width",
        ],
    )
    ax.set_yticks(y, labels, fontsize=9)
    for i, value in enumerate(values):
        label = f"{value:.2g}×"
        ax.text(value * 1.5, i, label, va="center", fontsize=9, color=INK)
    ax.set_xscale("log")
    ax.set_xlim(1e-7, 30)
    ax.set_xticks([1e-6, 1e-4, 1e-2, 1, 10], ["10⁻⁶", "10⁻⁴", "10⁻²", "1", "10"])
    ax.invert_yaxis()
    ax.set_xlabel("Sensitivity / tolerance   (log scale; ≤1 passes)")
    _note(
        ax,
        data.get(
            "diagnostic_note",
            "Numerical controls pass for the tested models.\n"
            "Active-space and all-background extraction tests reject.\n"
            "These are diagnostic ratios, not common uncertainty bars.",
        ),
    )


def _neutral(ax, data):
    from scipy.constants import c, physical_constants

    _panel(
        ax,
        "C",
        data.get("neutral_title", "A better neutral energy curve"),
        "Frozen-core CCSD(T): 26 geometries per basis",
    )
    for basis, color, marker, label in (
        ("aug-cc-pVQZ", ORANGE, "o", "aug-QZ"),
        ("aug-cc-pV5Z", BLUE, ".", "aug-5Z"),
    ):
        xy = np.asarray(data["neutral"][basis]["points_R_bohr_relative_E_ev"])
        ax.plot(xy[:, 0], xy[:, 1], marker=marker, ms=3.8, color=color, lw=1.5, label=label)
    ax.axvline(2.1323, color=GRAY, lw=0.8, linestyle=":")
    ax.legend(loc="upper center", frameon=False, ncol=2, fontsize=10)
    ax.set(xlim=(1.88, 2.52), ylim=(-0.10, 1.65), xlabel="Internuclear distance R / bohr")
    ax.set_ylabel("E(R) − E(2.1323), each basis's own zero / eV")
    checks = data["neutral_checks"]
    # Presentation boundary: one Debye is 1e-21/c coulomb metres.
    dipole = (
        abs(data["neutral"]["aug-cc-pV5Z"]["equilibrium_dipole_au"])
        * physical_constants["atomic unit of electric dipole mom."][0]
        * c
        * 1e21
    )
    ax.text(
        0.37,
        0.67,
        f"QZ → 5Z: ≤{checks['maximum_basis_relative_energy_difference_ev'] * 1e3:.2f} meV\n"
        f"Cubic held points: ≤{checks['maximum_cubic_held_error_ev'] * 1e3:.3f} meV\n"
        "Pilot budgets: 20 meV / 1 meV",
        transform=ax.transAxes,
        fontsize=9.5,
        color=GREEN,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#eef7f4", "edgecolor": "none"},
    )
    _note(
        ax,
        data.get(
            "neutral_note",
            f"CCSD density dipole: {dipole:.3f} D; "
            f"quoted experiment: {data['reference_dipole']['magnitude_debye']:.3f} D.\n"
            "Near-equilibrium basis / held-point checks pass.\n"
            "Next: correlation checks; unstable RHF at R = 3–4 bohr.",
        ),
    )


def _memory(ax, data):
    _panel(
        ax,
        "D",
        data.get("memory_title", "Larger models within local RAM"),
        "Measured kernel peaks; CAS12 dense is an allocation floor",
    )
    costs = data["sparse_costs"] + data.get("additional_sparse_costs", [])
    values = [row["kernel_peak_gib"] for row in costs]
    values += [
        data["dense_cas11"]["kernel_peak_gib"],
        data["memory_reference"]["cas12_dense_workspace_floor_gib"],
    ]
    count = len(costs)
    x = np.arange(count + 2)
    colors = [GREEN if row.get("qualified", True) else ORANGE for row in costs]
    bars = ax.bar(x, values, color=[*colors, BLUE, RED], width=data.get("memory_bar_width", 0.64))
    bars[-1].set_hatch("///")
    if "selected_cas11_roots" in data:
        selected = next(
            i for i, row in enumerate(costs) if row["roots"] == data["selected_cas11_roots"]
        )
        bars[selected].set_edgecolor(INK)
        bars[selected].set_linewidth(2)
        ax.text(selected, 5.2, "REUSE", ha="center", color="white", fontweight="bold", fontsize=9)
    ram = data["memory_reference"]["host_ram_gib"]
    ax.axhline(ram, color=RED, lw=1, linestyle="--")
    ax.text(0.015, ram, f"Sadaharu total RAM: {ram:.1f} GiB", fontsize=9, color=RED, va="bottom")
    ax.set_yscale("log")
    ax.set_ylim(3, 450)
    ax.set_yticks([5, 10, 25, 50, 100, 250], ["5", "10", "25", "50", "100", "250"])
    ax.set_xticks(
        x,
        [f"{row.get('model', 'CAS11')}\n{row['roots']}" for row in costs]
        + ["CAS11\ndense", "CAS12\ndense"],
        fontsize=data.get("memory_tick_fontsize", 10),
    )
    ax.set_ylabel("Peak / lower-bound memory / GiB (log scale)")
    for i, value in enumerate(values):
        ax.text(
            i, value * 1.14, f"{value:.1f}" if i < count + 1 else ">222", ha="center", fontsize=10
        )
    time_fontsize = data.get(
        "memory_time_fontsize", 8 if data.get("additional_sparse_costs") else 9
    )
    for i, row in enumerate(costs):
        ax.text(
            i, 3.4, f"{row['wall_hours']:.2f} h", ha="center", fontsize=time_fontsize, color="white"
        )
    ax.text(
        count,
        3.4,
        f"{data['dense_cas11']['wall_hours']:.2f} h*",
        ha="center",
        fontsize=time_fontsize,
        color="white",
    )
    _note(
        ax,
        data.get(
            "memory_note",
            f"2048–{data['sparse_costs'][-1]['roots']} roots reproduce dense scattering "
            "within gates.\n"
            "*Dense includes setup; sparse timings are replay costs.\n"
            + data.get(
                "memory_next_step",
                "Next: finish 16384 roots, then qualify a CAS12 sparse pilot.",
            ),
        ),
    )


def _tracking(ax, data):
    _panel(
        ax,
        "E",
        "Can we follow the same states?",
        "AO-following frame diagnostic; minimum over eight sectors",
    )
    pairs = data["overlaps"]["pair_order_R_bohr"]
    x = np.arange(len(pairs))
    for key, color, marker, label in (
        ("eight_root", ORANGE, "o", "8 roots: independently reconstructed"),
        (
            "twelve_root",
            BLUE,
            "s",
            data.get("twelve_root_label", "12 roots: producer checks pass"),
        ),
    ):
        values = data["overlaps"][key]["min_ao_following_singular_values"]
        ax.plot(x, values, marker=marker, ms=7, color=color, lw=1.6, label=label)
    ax.axhline(1, color=GRAY, linestyle="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(1e-13, 5)
    ax.set_xticks(x, ["1.900 → 2.1323", "2.1323 → 2.500", "1.900 → 2.500"], fontsize=9)
    ax.set_xlabel("Geometry pairs R / bohr   (lines guide the eye)")
    ax.set_ylabel("Minimum root-manifold singular value (log scale)")
    ax.legend(loc="center left", frameon=False, fontsize=9)
    _note(
        ax,
        "More roots recover missing overlap-manifold directions.\n"
        "Retained states still match outside the first five.\n"
        + data.get(
            "tracking_next_step", "Next: independent 12-root replay and finer geometry tracking."
        ),
    )


def _roadmap(ax, data):
    from matplotlib.patches import FancyBboxPatch

    _panel(
        ax,
        "F",
        data.get("roadmap_title", "Why the next experiments matter"),
        "Goal: a trustworthy neutral + resonance potential for dynamics",
    )
    ax.axis("off")
    entries = [
        (
            GREEN,
            "ESTABLISHED",
            "Reproducible engine, targets and sparse replay",
            "Independent roots / dipoles / phases; preserved failures",
        ),
        (
            RED if data.get("roadmap_import_status") == "STOPPED" else BLUE,
            data.get("roadmap_import_status", "ACTIVE"),
            data.get("roadmap_import_title", "Finish CAS11 → verify three CAS12 imports"),
            data.get("roadmap_import_detail", "16384 roots: B1 done, B2 running; imports queued"),
        ),
        (
            ORANGE,
            "NEXT",
            "Resolve state identity and extraction sensitivity",
            data.get(
                "roadmap_tracking_detail",
                "12-root replay, intermediate R, justified pole extraction",
            ),
        ),
        (
            ORANGE,
            "NEXT",
            "Balance the electronic model and neutral treatment",
            "Basis / active space / 40–50 ensemble / channels",
        ),
        (
            GRAY,
            "GATED",
            "Geometry sweep → fit → held-out scattering tests",
            "Then assess the 2-D potential's predictive accuracy",
        ),
    ]
    if "roadmap_entries" in data:
        colors = {"reuse": GREEN, "active": BLUE, "next": ORANGE, "gated": GRAY, "reject": RED}
        entries = [
            (colors[row["kind"]], row["status"], row["title"], row["detail"])
            for row in data["roadmap_entries"]
        ]
    for i, (color, status, title, detail) in enumerate(entries):
        y = 0.98 - i * 0.185
        ax.add_patch(
            FancyBboxPatch(
                (0, y - 0.155),
                1,
                0.15,
                boxstyle="round,pad=0.01,rounding_size=0.02",
                transform=ax.transAxes,
                facecolor=color,
                alpha=0.075,
                edgecolor="none",
            )
        )
        ax.text(
            0.025,
            y - 0.032,
            status,
            color=color,
            fontsize=8.5,
            fontweight="bold",
            transform=ax.transAxes,
        )
        ax.text(
            0.025,
            y - 0.078,
            title,
            color=INK,
            fontsize=10,
            fontweight="bold",
            transform=ax.transAxes,
        )
        ax.text(0.025, y - 0.124, detail, color=GRAY, fontsize=9, transform=ax.transAxes)
    _note(
        ax,
        data.get(
            "roadmap_note",
            "Current main blockers: model balance, state / pole identity,\n"
            "background-sensitive widths and stretched neutral reference.\n"
            "Literature agreement is a comparison, not a convergence test.",
        ),
    )


def _anchor_candidates(ax, data):
    _panel(
        ax,
        "A",
        "What the potential must reproduce",
        "CAS11 / DZ / 40: three independently checked pipelines",
    )
    rows = data["anchor_candidates"]
    x = np.asarray([row["R_bohr"] for row in rows])
    for key, color, marker, label in (
        ("position_ev", BLUE, "o", "Position above each neutral threshold"),
        ("full_width_ev", ORANGE, "s", "Full width"),
    ):
        values = [row[key] for row in rows]
        ax.scatter(x, values, c=color, marker=marker, s=65, label=label, zorder=5)
        for r, value in zip(x, values, strict=True):
            ax.annotate(
                f"{value:.3f}",
                (r, value),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                color=color,
                fontsize=10,
            )
    ax.set(
        xlim=(1.85, 2.56),
        ylim=(0, 4.35),
        xlabel="Internuclear distance R / bohr",
        ylabel="Native fit candidate / eV",
    )
    ax.set_xticks(x, ["1.900", "2.1323", "2.500"])
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    _note(
        ax,
        data.get(
            "anchor_note",
            "Reuse: anchor geometries, raw phases and import checks.\n"
            "Markers are candidates; state / pole identity remains open.\n"
            "Widths need background controls before fitting a potential.",
        ),
    )


def _channel_coverage(ax, data):
    _panel(
        ax,
        "E",
        "Which target channels are missing?",
        "Physical A2 projection; four links on two tracked routes",
    )
    rows = data["a2_coverage_links"]
    y = np.arange(len(rows))
    for offset, key, color, label in (
        (-0.23, "retained_five_projection_squared", RED, "Retained roots 1–5"),
        (0, "first_six_projection_squared", ORANGE, "Diagnostic roots 1–6"),
        (0.23, "first_twelve_projection_squared", BLUE, "Diagnostic roots 1–12"),
    ):
        values = [row[key] for row in rows]
        ax.barh(y + offset, values, height=0.20, color=color, label=label)
        for i, value in enumerate(values):
            if offset == -0.23:
                ax.text(0.015, i + offset, f"{value:.6f}", va="center", color=RED, fontsize=8)
    ax.set_yticks(y, [row["label"] for row in rows], fontsize=8.5)
    ax.set(
        xlim=(0, 1),
        xlabel=r"$\sum |\langle$left state | right root$\rangle|^2$ over selected roots",
    )
    ax.invert_yaxis()
    ax.set_ylim(4.5, -0.5)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    omitted = data["a2_singlet_midpoint_omitted_weights"]
    _note(
        ax,
        data.get(
            "channel_coverage_note",
            f"On the singlet 5 → 6 link: root7 = {omitted['7']:.3f}, root9 = {omitted['9']:.3f}.\n"
            "Next: fixed-orbital explicit-A2-channel scattering control.\n"
            "Projection is a coverage diagnostic, not a continuity verdict.",
        ),
    )


def _anchor_channel_sensitivity(ax, data):
    _panel(
        ax,
        "E",
        "Where do channels change scattering?",
        "Five → nine A2 per spin; phase sensitivity across the full grid",
    )
    rows = data["anchor_channel_sensitivity"]
    y = np.arange(len(rows))
    values = [row["maximum_phase_difference_rad_mod_pi"] for row in rows]
    budget = data["channel_phase_budget_rad_mod_pi"]
    ax.axvspan(0, budget, color=GREEN, alpha=0.05)
    ax.axvline(budget, color=INK, linestyle="--", lw=1.1, label="0.05-rad budget")
    ax.barh(y, values, color=[GREEN if row["passed"] else RED for row in rows], height=0.5)
    for i, value in enumerate(values):
        ax.text(value + 0.001, i, f"{value:.5f}", va="center", fontsize=10)
    ax.set_yticks(y, [row["label"] for row in rows], fontsize=9)
    ax.set(xlim=(0, 0.085), ylim=(2.6, -0.6), xlabel="Maximum phase change modulo π / rad")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    _note(ax, data["channel_sensitivity_note"])


def _identifiability(ax, data):
    _panel(
        ax,
        "B",
        "Can good phases hide wrong widths?",
        "Known ±10% widths; every unused phase audit within 0.001 rad",
    )
    rows = data["identifiability_scales"]
    y, left = np.arange(len(rows)), np.zeros(len(rows))
    for degree, color in ((3, GRAY), (7, ORANGE), (11, BLUE)):
        values = np.asarray([row["width_counterexamples_by_degree"][str(degree)] for row in rows])
        ax.barh(y, values, left=left, color=color, label=f"Degree {degree}", height=0.52)
        left += values
    for i, value in enumerate(left):
        ax.text(value + 0.8, i, f"{int(value)} / 54", va="center", fontsize=9)
    ax.set_yticks(y, [row["label"] for row in rows], fontsize=9)
    ax.set(
        xlim=(0, 66), ylim=(4.8, -0.8), xlabel="Known wrong-width counterexamples / tested cells"
    )
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    _note(ax, data["identifiability_note"])


def _physical_poles(ax, data):
    _panel(
        ax,
        "B",
        "Can physical outgoing poles be recovered?",
        "Real shell Hamiltonian; independent exact Jost / Lambert-W oracle",
    )
    rows = data["physical_shell_convergence"]
    x = np.arange(len(rows))
    for key, color, marker, label in (
        ("maximum_s_error", BLUE, "o", "Real-axis S error"),
        ("maximum_q_error", GREEN, "s", "Outgoing pole q error"),
    ):
        values = [r[key] for r in rows]
        ax.semilogy(x, values, color=color, marker=marker, lw=1.8, label=label)
        for i, value in enumerate(values):
            ax.annotate(
                f"{value:.2e}",
                (i, value),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=8.5,
                color=color,
            )
    ax.axhline(1e-7, color=GRAY, ls="--", lw=1, label="Finest S / q budget")
    ax.set_xticks(x, [f"{r['ode_rtol']:.0e}" for r in rows])
    ax.set(
        xlim=(-0.3, 2.3),
        ylim=(1e-13, 2e-3),
        xlabel="ODE relative tolerance (three declared levels)",
        ylabel="Maximum absolute error (dimensionless)",
    )
    ax.legend(loc="lower left", frameon=False, fontsize=8.5)
    _note(ax, data["physical_shell_note"])


def _tail_cutoff(ax, data):
    _panel(
        ax,
        "B",
        "Can a unitary boundary still be wrong?",
        "Exact two-wave dipole-tail model; omit coupling beyond matching radius",
    )
    rows = data["tail_cutoff_controls"]
    x = [r["radius_ratio"] for r in rows]
    y = [r["s_entry_error"] for r in rows]
    ax.loglog(x, y, color=BLUE, marker="o", lw=1.8, label="Error versus uncut exact S")
    for radius, error in zip(x, y, strict=True):
        ax.annotate(
            f"{error:.5f}",
            (radius, error),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=BLUE,
        )
    ax.set_xticks(x, [str(int(v)) for v in x])
    ax.set(
        xlim=(5, 200),
        ylim=(0.001, 0.09),
        xlabel="Matching radius B / inner radius a",
        ylabel="Maximum absolute S-entry error",
    )
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    _note(ax, data["tail_cutoff_note"])


def _combined_refinement(ax, data):
    _panel(
        ax,
        "B",
        data.get("combined_refinement_title", "Combined polar exterior: refinement rejects"),
        "Unequal thresholds + inverse-square coupling; outgoing boundary map Y",
    )
    x = np.arange(2)
    for method, color, marker, label in [
        ("ode", BLUE, "o", "Riccati ODE: rtol 1e-5 / 1e-8 / 1e-11"),
        (
            "spectral",
            ORANGE,
            "s",
            data.get("combined_spectral_label", "Spectral BVP: degrees 64 / 128 / 256"),
        ),
    ]:
        row = data["combined_refinement"][method]
        values = [row["coarse_to_middle"], row["middle_to_finest"]]
        ax.semilogy(x, values, color=color, marker=marker, lw=1.8, label=label)
        for i, value in enumerate(values):
            ax.annotate(
                f"{value:.2e}",
                (i, value),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=9,
                color=color,
            )
    ax.axhline(1e-7, color=GRAY, ls="--", lw=1, label="Final refinement budget: 1e-7")
    ax.set_xticks(x, ["Coarse → middle", "Middle → finest"])
    ax.set(xlim=(-0.25, 1.35), ylim=(1e-9, 0.8), ylabel="Maximum absolute change in Y entry")
    ax.legend(loc="upper right", frameon=False, fontsize=8.2)
    _note(ax, data["combined_refinement_note"])


def _native_exterior_convergence(ax, data):
    _panel(
        ax,
        "B",
        "Independent exterior: fourth-order control",
        "Smooth coupled toy intervals; independent direct ODE reference",
    )
    x = np.arange(3)
    for row, color, marker in zip(
        data["native_exterior_convergence"], [BLUE, ORANGE], ["o", "s"], strict=True
    ):
        values = row["maximum_k_errors"]
        label = "Two open channels" if row["thresholds"][1] < 1 else "One open / one closed"
        ax.semilogy(x, values, color=color, marker=marker, lw=1.8, label=label)
        for i, value in enumerate(values):
            ax.annotate(
                f"{value:.2e}",
                (i, value),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=9,
                color=color,
            )
    ax.set_xticks(x, ["0.2", "0.1", "0.05"])
    ax.set(xlim=(-0.25, 2.35), ylim=(5e-8, 1e-3), xlabel="Radial step / bohr")
    ax.set_ylabel("Maximum absolute K-entry error")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    _note(ax, data["native_exterior_note"])


def _native_exterior_pilot(ax, data):
    _panel(
        ax,
        "B",
        "Native exterior: reproduction rejects",
        "Full-channel equilibrium-B1 pilot; each error / its 1e-6 budget",
    )
    x = np.arange(3)
    for metric, color, marker, label in [
        ("native_k_error", BLUE, "o", "K entries"),
        ("native_s_error", ORANGE, "s", "S entries"),
        ("native_phase_error_rad_mod_pi", RED, "^", "Phase sum modulo π"),
    ]:
        values = [row[metric] / 1e-6 for row in data["native_exterior_pilot_finest"]]
        ax.semilogy(x, values, color=color, marker=marker, lw=1.6, label=label)
        for i, value in enumerate(values):
            ax.annotate(
                f"{value:.3g}×",
                (i, value),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=8.5,
                color=color,
            )
    ax.axhline(1, color=GRAY, ls="--", lw=1, label="Reproduction budget")
    ax.set_xticks(x, ["0.1", "2.5", "5.0"])
    ax.set(xlim=(-0.25, 2.35), ylim=(1e-3, 10), xlabel="Requested energy / eV")
    ax.set_ylabel("Finest error / declared tolerance (≤1 passes)")
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    _note(ax, data["native_exterior_pilot_note"])


def _snapshot_data(snapshot: Path) -> dict:
    """Resolve a full snapshot or a delta against a hash-pinned earlier release."""
    data = json.loads(snapshot.read_text())
    if "base_snapshot" not in data:
        return data
    base = snapshot.with_name(data["base_snapshot"]).read_bytes()
    if hashlib.sha256(base).hexdigest() != data["base_snapshot_sha256"]:
        raise ValueError("Historical base snapshot hash differs")
    result = _snapshot_data(snapshot.with_name(data["base_snapshot"]))
    result.update(data["updates"])
    result["source_file_sha256"].update(data["additional_source_file_sha256"])
    return result


def render(snapshot: Path, output: Path) -> None:
    """Write a PNG and text-preserving SVG from the dated plotted-data snapshot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = _snapshot_data(snapshot)
    for name, digest in data["source_file_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Snapshot provenance differs from current source: {name}")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": "#b6c0ca",
            "axes.facecolor": "white",
            "figure.facecolor": "#f7f9fb",
            "xtick.color": INK,
            "ytick.color": INK,
            "svg.fonttype": "none",
            "svg.hashsalt": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(20, 13.3))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.84, bottom=0.18, wspace=0.42, hspace=0.70)
    fig.text(
        0.04,
        0.966,
        "CO scattering: reliable calculations, an electronic model still to qualify",
        fontsize=23,
        fontweight="bold",
    )
    fig.text(
        0.04,
        0.934,
        data.get(
            "figure_subtitle",
            "What we have achieved • Why agreement alone is insufficient • "
            "What unlocks a predictive 2-D potential",
        ),
        fontsize=13,
        color=GRAY,
    )
    date = data["snapshot_utc"].split(".")[0].replace("T", " ")
    fig.text(
        0.04,
        0.906,
        f"Evidence snapshot: {date} UTC   |   Blue: our results   Green: established checks   "
        "Orange/red: sensitivity or open gates",
        fontsize=10,
        color=GRAY,
    )
    panels = (
        _anchor_candidates if "anchor_candidates" in data else _resonances,
        (
            _native_exterior_pilot
            if "native_exterior_pilot_finest" in data
            else _native_exterior_convergence
            if "native_exterior_convergence" in data
            else _combined_refinement
            if "combined_refinement" in data
            else _tail_cutoff
            if "tail_cutoff_controls" in data
            else _physical_poles
            if "physical_shell_convergence" in data
            else _identifiability
            if "identifiability_scales" in data
            else _diagnostics
        ),
        _neutral,
        _memory,
        (
            _anchor_channel_sensitivity
            if "anchor_channel_sensitivity" in data
            else _channel_coverage
            if "a2_coverage_links" in data
            else _tracking
        ),
        _roadmap,
    )
    for ax, draw in zip(axes.flat, panels, strict=True):
        draw(ax, data)
    fig.text(
        0.04,
        0.053,
        data.get(
            "source_footer",
            "REFERENCE LOCATORS   Laporta et al. 2012: preprint p.4 (unadjusted width) • "
            "Dora et al. 2016: p.6, Table 4; dipole p.4, Table 2 • "
            "Dora & Tennyson 2020: p.4, Table 2",
        ),
        fontsize=9,
        color=GRAY,
    )
    fig.text(
        0.04,
        0.030,
        "Atomic units internally; eV at presentation. C2v state counts are components. "
        "AO-following transport is not physical wavefunction overlap. "
        + data.get(
            "verification_footer",
            "12-root / latest B1 values await independent publication.",
        ),
        fontsize=9,
        color=GRAY,
    )
    fig.text(
        0.04,
        0.012,
        data.get(
            "key_footer",
            "KEY   CAS10/11/12: ten active electrons in 10/11/12 orbitals • "
            "DZ/TZ: double-/triple-zeta basis • "
            "39/60: SEP virtual orbitals • 40/50/41: retained C2v target components",
        ),
        fontsize=8.5,
        color=GRAY,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=180)
    svg = StringIO()
    fig.savefig(svg, format="svg", metadata={"Date": None})
    output.with_suffix(".svg").write_text(
        "\n".join(line.rstrip() for line in svg.getvalue().splitlines()) + "\n"
    )
    plt.close(fig)


def main() -> None:
    """Render the dated overview with its declared evidence snapshot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument(
        "--output", type=Path, default=ROOT / "docs/physics/figures/co-progress-20261017T1340"
    )
    args = parser.parse_args()
    render(args.snapshot, args.output)


if __name__ == "__main__":
    main()
