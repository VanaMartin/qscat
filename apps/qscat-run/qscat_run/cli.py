"""The `qscat-run` click command group: `validate`, `list`, `init`, `run`, `fetch`.

`run` resolves a config, runs it (`runner.run_experiment`), and writes its
artifacts (`artifacts.write_artifacts`) -- both `ti` and `td` methods (any
subset, including both at once) actually solve. `--dry-run` resolves +
prints the plan (grids, sizes, energy count) without solving anything.
"""

from __future__ import annotations

import warnings
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import click
import yaml

from qscat_run import presets
from qscat_run.artifact_store import ArtifactStoreError, fetch, load_pointer
from qscat_run.artifacts import write_artifacts
from qscat_run.config import ConfigError, load_config, validate_config
from qscat_run.runner import run_experiment

__all__ = ["main"]


def _resolve_method_grids(cfg: Any, method: str) -> str:
    """Build `method`'s grids and return a one-line description of them.

    `ti`/`td` share the tensor-grid resolver; the two 1-D-nuclear
    approximations do not. `lcp` and `nrm` both need the preset's paired
    two-ECS-angle electronic decks plus the fine nuclear deck, which the
    ti/td `TensorGrid` shape cannot express, so each has its own resolver --
    routing them through `resolve_grid` would make `validate`/`--dry-run`
    reject configs `run` supports.
    """
    if method in ("lcp", "nrm"):
        # NOTE the order: both resolvers return (nuclear, elec_a, elec_b),
        # nuclear FIRST.
        resolver = presets.resolve_lcp_grids if method == "lcp" else presets.resolve_nrm_grids
        g_nuc, g_ea, _g_eb = resolver(cfg)
        return f"nuclear={g_nuc.points.size} electronic={g_ea.points.size} (x2 ECS angles)"
    tg = presets.resolve_grid(cfg, method)
    return f"shape={tg.shape} size={tg.size}"


@click.group()
def main() -> None:
    """qscat-run -- a config-driven CLI for 2-D electron-diatomic model experiments."""


@main.command("validate")
@click.argument("config_path", metavar="CONFIG", type=click.Path(exists=True, dir_okay=False))
def validate_cmd(config_path: str) -> None:
    """Parse + validate CONFIG: schema, the (molecule, observable) validity
    matrix, and preset/grid resolution. No solve -- fast, for CI/edit loops.
    """
    cfg = load_config(config_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        validate_config(cfg)
    for w in caught:
        click.echo(f"WARNING: {w.message}", err=True)

    resolved = presets.resolve_defaults(cfg)
    for method in resolved.methods:
        try:
            _resolve_method_grids(resolved, method)
        except ConfigError:
            raise
        except Exception as exc:
            raise ConfigError(f"failed to build the {method!r} grid: {exc}") from exc

    click.echo(f"OK: {config_path} is valid (molecule={cfg.molecule}, methods={list(cfg.methods)})")


@main.command("list")
def list_cmd() -> None:
    """List known molecules, their preset variants, and the (molecule,
    observable) validity matrix."""
    click.echo("Molecules (preset variants -- valid observables):")
    for molecule in sorted(presets.MODELS):
        variants = sorted(presets.available_presets(molecule))
        kinds = sorted(presets.VALIDITY.get(molecule, frozenset()))
        warn_kinds = sorted(presets.WARN_OBSERVABLES.get(molecule, frozenset()))
        note = f"  (closed-in-range: {warn_kinds})" if warn_kinds else ""
        click.echo(f"  {molecule:<4} presets={variants} observables={kinds}{note}")


def _default_channel_count(kind: str, preset: presets.MoleculePreset) -> int:
    if kind == "ve":
        return max(1, min(preset.n_vib - 1, 3))
    if kind == "dr":
        return 3
    return 1


# Time-domain evolution defaults per molecule for the `init --methods td`
# scaffold. Not part of `presets.MoleculePreset` (the brief's field list has
# no dt/n_steps/order) -- N2/H2P mirror their validated decks (N2's
# `TD_WORKING_GRID`; H2P's eMoScat DR evolution: dt=10, order-3); NO/F2 have no validated TD
# experiment yet (see CLAUDE.md's diatomic note), so their numbers are a
# documented, unvalidated placeholder.
_TD_DEFAULTS: dict[str, dict[str, float | int]] = {
    "N2": {"dt": 1.0, "n_steps": 1500, "order": 3},
    "NO": {"dt": 1.0, "n_steps": 1000, "order": 3},
    "F2": {"dt": 1.0, "n_steps": 1000, "order": 3},
    "H2P": {"dt": 10.0, "n_steps": 2000, "order": 3},
}


def _scaffold_yaml(molecule: str, obs_kinds: list[str], methods: list[str]) -> str:
    variant = presets.DEFAULT_PRESET
    preset = presets.PRESETS.get(f"{molecule}:{variant}")
    if preset is None:
        variant = next(iter(sorted(presets.available_presets(molecule))), presets.DEFAULT_PRESET)
        preset = presets.PRESETS.get(f"{molecule}:{variant}")
    if preset is None:  # pragma: no cover -- every registered molecule has >=1 preset
        raise ConfigError(f"no preset available for molecule {molecule!r}")

    observables = [
        {"kind": kind, "channels": _default_channel_count(kind, preset)} for kind in obs_kinds
    ]
    cfg: dict[str, Any] = {
        "molecule": molecule,
        "methods": methods,
        "observables": observables,
        # Scaffolded in the same form the schema stores, so a generated config
        # and a resolved one look alike and can be diffed against each other.
        "energies": {
            "ranges": [
                {"start": r.start, "stop": r.stop, "step": r.step}
                for r in (preset.default_energies.ranges or ())
            ]
        },
        "grid": {"preset": preset.variant},
        "v_init": 0,
        "artifacts": {"cross_section": True},
        "backend": "auto",
        "output_dir": f"runs/{molecule.lower()}",
    }
    if "td" in methods:
        td_defaults = _TD_DEFAULTS.get(molecule, {"dt": 1.0, "n_steps": 1000, "order": 3})
        # Per-observable-kind mapping (`ve` electronic / `da`/`dr` nuclear --
        # see `presets.MoleculePreset`'s docstring): only emit an entry for a
        # kind this run actually requests AND the preset has a default for.
        kind_test_functions = {
            "ve": preset.ve_test_function,
            "da": preset.da_test_function,
            "dr": preset.dr_test_function,
        }
        test_function_map = {
            kind: {
                "r0_out": tf.r0_out,
                "p0_out": tf.p0_out,
                "sigma_out": tf.sigma_out,
            }
            for kind in obs_kinds
            if (tf := kind_test_functions.get(kind)) is not None
        }
        cfg["td"] = {
            **td_defaults,
            "extractors": ["flow", "delta", "tw"],
            "incident": {
                "r0": preset.default_incident.r0,
                "p0": preset.default_incident.p0,
                "sigma": preset.default_incident.sigma,
            },
            "test_function": test_function_map,
        }

    valid_kinds = sorted(presets.VALIDITY.get(molecule, frozenset()))
    variants = sorted(presets.available_presets(molecule))
    header = (
        f"# qscat-run starter config for {molecule}, scaffolded by `qscat-run init`.\n"
        f"# molecule: N2 | NO | F2 | H2P\n"
        f"# methods: any subset of [ti, td, lcp, nrm] (lcp = local-complex-potential DA,\n"
        f"#          F2/NO only; nrm = nonlocal resonance model VE + DA, N2/NO/F2)\n"
        f"# observables: a list of {{kind, channels}}; {molecule} supports kind in {valid_kinds}\n"
        f"# grid.preset: one of {variants} (or an explicit {{electronic, nuclear}} grid)\n"
        "# See docs/superpowers/specs/2026-08-01-qscat-run-cli-design.md for the full schema.\n"
    )
    body = yaml.safe_dump(cfg, sort_keys=False)
    return header + body


@main.command("init")
@click.argument("molecule")
@click.option(
    "--observables",
    default=None,
    help=(
        "Comma-separated observable kinds, e.g. 've,da' (default: the molecule's first valid kind)."
    ),
)
@click.option(
    "--methods",
    default="ti",
    help="Comma-separated methods, e.g. 'ti,td' or 'ti,lcp,nrm' (default: 'ti').",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to write the scaffolded YAML config.",
)
def init_cmd(molecule: str, observables: str | None, methods: str, output_path: str) -> None:
    """Scaffold a fully-commented starter config for MOLECULE, seeded from
    its preset (N2 | NO | F2 | H2P)."""
    molecule = molecule.upper()
    if molecule not in presets.MODELS:
        raise ConfigError(f"unknown molecule {molecule!r}; choose one of {sorted(presets.MODELS)}")

    method_list = [m.strip() for m in methods.split(",") if m.strip()] or ["ti"]
    if observables:
        obs_kinds = [o.strip() for o in observables.split(",") if o.strip()]
    else:
        valid_kinds = sorted(presets.VALIDITY.get(molecule, frozenset()))
        obs_kinds = valid_kinds[:1]

    text = _scaffold_yaml(molecule, obs_kinds, method_list)
    Path(output_path).write_text(text)
    click.echo(f"wrote {output_path}")


@main.command("run")
@click.argument("config_path", metavar="CONFIG", type=click.Path(exists=True, dir_okay=False))
@click.option("--output", "output_dir", default=None, help="Override the config's output_dir.")
@click.option(
    "--backend",
    type=click.Choice(["auto", "mumps", "scipy"]),
    default=None,
    help="Override the config's backend.",
)
@click.option("--dry-run", is_flag=True, help="Resolve + print the plan without solving.")
def run_cmd(config_path: str, output_dir: str | None, backend: str | None, dry_run: bool) -> None:
    """Run CONFIG end-to-end: parse -> validate -> resolve -> solve -> write
    artifacts. `--dry-run` resolves and prints the plan without solving or
    writing anything.
    """
    cfg = load_config(config_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        validate_config(cfg)
    for w in caught:
        click.echo(f"WARNING: {w.message}", err=True)

    if backend is not None:
        cfg = replace(cfg, backend=backend)
    if output_dir is not None:
        cfg = replace(cfg, output_dir=output_dir)

    resolved = presets.resolve_defaults(cfg)

    if dry_run:
        click.echo(f"molecule: {resolved.molecule}")
        click.echo(f"methods: {list(resolved.methods)}")
        click.echo("observables:")
        for obs in resolved.observables:
            click.echo(f"  - kind={obs.kind} channels={obs.channels}")
        n_energies = len(resolved.energies.as_array()) if resolved.energies is not None else 0
        click.echo(f"energies: {n_energies}")
        for method in resolved.methods:
            try:
                described = _resolve_method_grids(resolved, method)
            except ConfigError:
                raise
            except Exception as exc:
                raise ConfigError(f"failed to build the {method!r} grid: {exc}") from exc
            click.echo(f"grid[{method}]: {described}")
        click.echo(f"output_dir: {resolved.output_dir}")
        click.echo(f"backend: {resolved.backend}")
        return

    result = run_experiment(resolved)
    out_dir = resolved.output_dir
    timestamp = datetime.now(UTC).isoformat()
    write_artifacts(result, resolved, out_dir, timestamp=timestamp)
    click.echo(f"wrote artifacts to {out_dir}")


@main.command("fetch")
@click.argument("directory", metavar="DIR", type=click.Path(exists=True, file_okay=False), nargs=-1)
@click.option("--list", "list_only", is_flag=True, help="Print the URLs without downloading.")
def fetch_cmd(directory: tuple[str, ...], list_only: bool) -> None:
    """Download the published artifacts DIR points at.

    Large run outputs are not committed: a sweep is reproducible from its
    `config.resolved.yaml`, but reproducing costs minutes to hours, so the
    bytes live in public object storage and DIR carries a small
    `artifacts.json` naming them. Reads are anonymous -- no account, no
    credentials.

    Every file is verified against the digest recorded at publication; one
    already present and correct is skipped, so re-running is free.
    """
    if not directory:
        raise click.UsageError("give at least one run directory")
    for d in directory:
        try:
            pointer = load_pointer(d)
            if list_only:
                click.echo(f"{d}  (from {pointer.git_sha[:7]})")
                for name in sorted(pointer.artifacts):
                    click.echo(f"  {pointer.url_for(name)}")
                continue
            written = fetch(d)
        except ArtifactStoreError as exc:
            raise click.ClickException(str(exc)) from exc
        if written:
            click.echo(f"{d}: fetched {len(written)} file(s)")
        else:
            click.echo(f"{d}: already complete")
