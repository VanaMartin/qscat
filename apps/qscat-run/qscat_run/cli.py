"""The `qscat-run` click command group: `validate`, `list`, `init`, `run`.

`run` is a stub in this task (Task 1) -- the shared-work runner (`runner.py`)
and artifact writers (`artifacts.py`) land in Task 2/3; it raises a
`click.ClickException` pointing at `validate` in the meantime.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import click
import yaml

from qscat_run import presets
from qscat_run.config import ConfigError, load_config, validate_config

__all__ = ["main"]


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
            presets.resolve_grid(resolved, method)
        except ConfigError:
            raise
        except Exception as exc:  # noqa: BLE001 -- surface as an actionable ConfigError
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
# no dt/n_steps/order) -- N2/H2P mirror their validated decks
# (`TD_WORKING_GRID`, `validation.h2plus.td_dr`); NO/F2 have no validated TD
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
        "energies": {
            "min": preset.default_energies.min,
            "max": preset.default_energies.max,
            "step": preset.default_energies.step,
        },
        "grid": {"preset": preset.variant},
        "v_init": 0,
        "artifacts": {"cross_section": True},
        "backend": "auto",
        "output_dir": f"runs/{molecule.lower()}",
    }
    if "td" in methods:
        td_defaults = _TD_DEFAULTS.get(molecule, {"dt": 1.0, "n_steps": 1000, "order": 3})
        cfg["td"] = {
            **td_defaults,
            "extractors": ["flow", "delta", "tw"],
            "incident": {
                "r0": preset.default_incident.r0,
                "p0": preset.default_incident.p0,
                "sigma": preset.default_incident.sigma,
            },
            "test_function": {
                "r0_out": preset.default_test_function.r0_out,
                "p0_out": preset.default_test_function.p0_out,
                "sigma_out": preset.default_test_function.sigma_out,
            },
        }

    valid_kinds = sorted(presets.VALIDITY.get(molecule, frozenset()))
    variants = sorted(presets.available_presets(molecule))
    header = (
        f"# qscat-run starter config for {molecule}, scaffolded by `qscat-run init`.\n"
        f"# molecule: N2 | NO | F2 | H2P\n"
        f"# methods: any subset of [ti, td]\n"
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
        "Comma-separated observable kinds, e.g. 've,da' "
        "(default: the molecule's first valid kind)."
    ),
)
@click.option(
    "--methods",
    default="ti",
    help="Comma-separated methods, e.g. 'ti,td' (default: 'ti').",
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
    """Run CONFIG end-to-end: parse -> resolve -> solve -> write artifacts.

    Not yet implemented -- the shared-work runner (`runner.py`) and artifact
    writers (`artifacts.py`) are a later task; use `qscat-run validate` to
    check a config in the meantime.
    """
    raise click.ClickException(
        "run is implemented in a later task (the runner/artifacts writers); "
        "use `qscat-run validate CONFIG` to check a config for now."
    )
