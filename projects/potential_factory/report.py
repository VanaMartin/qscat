"""FitReport: the model, per-tier residuals, verdicts, and provenance."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

__all__ = ["FitReport", "TierResult", "Tolerances", "ecs_bounded"]

Status = Literal["met", "not met", "not attempted"]


@dataclass(frozen=True)
class TierResult:
    name: str
    status: Status
    rms: float
    max: float
    detail: str


@dataclass(frozen=True)
class Tolerances:
    """PLACEHOLDERS (Hartree / relative) until a measured sensitivity
    budget replaces them."""

    v0_rms: float = 2e-4
    omega_e_rel: float = 0.01
    e_res_rms: float = 1e-3
    gamma_rel: float = 0.10
    gamma_floor: float = 2e-3
    coupling_log_rms: float = 0.2


@dataclass
class FitReport:
    target_name: str
    parameters: dict[str, float]
    tiers: list[TierResult]
    ecs_bounds_deg: dict[str, float]
    crossing_R: float | None
    da_threshold_sign: int | None
    provenance: dict[str, dict[str, str]] = field(default_factory=dict)

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def from_json(cls, path: str | Path) -> FitReport:
        d = json.loads(Path(path).read_text())
        d["tiers"] = [TierResult(**t) for t in d["tiers"]]
        return cls(**d)


def ecs_bounded(model, pair, R_tail, *, nuclear_deg: float | None = None) -> dict[str, float]:
    """Growth of `|V|` on the ECS tails relative to the real region; raises if > 10x.

    The returned dict mixes two different kinds of number, and says which is
    which:

    `electronic_max_deg` (45) and `nuclear_max_deg` (90) are ANALYTIC bounds
    of the ansatz, not measurements. Under `r -> r e^{i theta}` the Gaussian
    well `e^{-alpha r^2}` picks up `r^2 = |r|^2 e^{2 i theta}`, so its
    exponent is `-alpha |r|^2 (cos 2theta + i sin 2theta)` and it decays iff
    `cos 2theta > 0`, i.e. `theta < 45 deg`. The EMO neutral curve's
    exponentials `e^{-beta (R - R_e)}` pick up only `e^{i theta}` and stay
    bounded iff `cos theta > 0`, i.e. `theta < 90 deg`. Both hold for every
    member of this ansatz family, independently of any grid.

    `probed_electronic_deg` and `probed_nuclear_deg` are the angles this call
    actually EVALUATED at -- the ECS tail angle of `pair.grid_a` and the
    direction of `R_tail` -- and `tail_growth` is the one measured quantity:
    `max|V|` over those two tails divided by `max|V|` over the real region.
    Note only `grid_a` is probed; a `pair`'s second grid sits at a different
    angle by construction, so a `probed_electronic_deg` well under 45 deg is
    not on its own a statement about the pair.

    `nuclear_deg` may be passed when the caller built `R_tail` from a known
    angle; otherwise it is derived from the tail's own direction, which needs
    at least two points.
    """
    r = pair.grid_a.points
    real = r.imag == 0.0
    tail = np.asarray(R_tail, dtype=np.complex128)
    v_r = np.abs(model.surface(r, model.R_e))
    v_R = np.abs(model.v0(tail))
    ref = max(float(np.max(v_r[real][r[real].real > 0.3])), 1e-12)
    growth = max(float(np.max(v_r[~real])), float(np.max(v_R))) / ref
    if growth > 10.0:
        raise ValueError(f"potential grows {growth:.1f}x on the ECS tail; not absorbing")
    # A GridSpec is validated to carry at most one distinct nonzero tail
    # angle, so the maximum over its elements IS the ECS angle (0 for a grid
    # with no complex tail at all).
    probed_e = max(float(el.angle_deg) for el in pair.grid_a.spec.elements)
    if nuclear_deg is None:
        if tail.size < 2:
            raise ValueError("nuclear_deg must be given when R_tail has fewer than two points")
        nuclear_deg = float(np.rad2deg(np.angle(tail[-1] - tail[0])))
    return {
        "electronic_max_deg": 45.0,
        "nuclear_max_deg": 90.0,
        "probed_electronic_deg": probed_e,
        "probed_nuclear_deg": float(nuclear_deg),
        "tail_growth": growth,
    }
