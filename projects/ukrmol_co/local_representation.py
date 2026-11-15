"""Evaluate the frozen local CO curves; no fit, extrapolation or engine calls.

See docs/physics/co-effective-resonance-target.md for scope and conventions.
Downloads are explicit (README); archived source is never extracted or executed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import tarfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.typing import ArrayLike, NDArray

DELIVERY_URL = (
    "https://data.qscat.org/co-local-representation-delivery-20261115T2020/"
    "delivery-capsule.tar.437ea9e58c35.gz"
)
DELIVERY_SHA256 = "437ea9e58c351b2a940bd2f57be05f2c09d34c3bb0c906735bd273b419c53555"
DOMAIN_BOHR = (1.9, 2.5)
_CURVE_MEMBERS = (
    "representation.json",
    "inputs/frozen-c2-model.json",
    "inputs/accepted-neutral-result.json",
)
_EVIDENCE_MEMBERS = (
    *_CURVE_MEMBERS,
    "ten-qualified-pairs.json",
    "inputs/lower-frozen-score.json",
    "inputs/upper-frozen-score.json",
)


def _real_finite(values: ArrayLike) -> NDArray[np.float64]:
    raw = np.asarray(values)
    if raw.dtype.kind not in "fiu" or not np.all(np.isfinite(raw)):
        raise ValueError("Expected finite real numeric values, not complex/bool/string values")
    return np.array(raw, dtype=float, copy=True)


def _radius(radius_bohr: ArrayLike) -> NDArray[np.float64]:
    radius = _real_finite(radius_bohr)
    if np.any((radius < DOMAIN_BOHR[0]) | (radius > DOMAIN_BOHR[1])):
        raise ValueError("R must lie in [1.9, 2.5] bohr; no extrapolation or clipping")
    return radius


def _cubic(knots: NDArray, coefficients: NDArray, radius: NDArray) -> NDArray:
    # Only the interval index is capped, so the right endpoint uses the last piece.
    index = np.minimum(np.searchsorted(knots, radius, side="right") - 1, len(knots) - 2)
    dx = radius - knots[index]
    value = np.zeros_like(radius)
    for row in coefficients:
        value = value * dx + row[index]
    return value


def _scalar_or_array(value: NDArray) -> float | complex | NDArray:
    return value.item() if np.ndim(value) == 0 else value


@dataclass(frozen=True, eq=False)
class LocalRepresentation:
    """Finite real-bohr CO curves, in Hartree; not a global scattering model.

    Construction from arrays is useful for synthetic controls, not scientific
    acceptance. Use load_representation for the authenticated accepted delivery.
    """

    neutral_coefficients: ArrayLike
    knots_bohr: ArrayLike
    er_coefficients: ArrayLike
    log_gamma_coefficients: ArrayLike

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            value = _real_finite(getattr(self, name))
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        knots = self.knots_bohr
        if (
            self.neutral_coefficients.shape != (7,)
            or knots.ndim != 1
            or len(knots) < 4
            or tuple(knots[[0, -1]]) != DOMAIN_BOHR
            or np.any(np.diff(knots) <= 0)
            or self.er_coefficients.shape != (4, len(knots) - 1)
            or self.log_gamma_coefficients.shape != (4, len(knots) - 1)
        ):
            raise ValueError("Expected degree-6 neutral and ordered cubic pieces over [1.9, 2.5]")

    def v0(self, radius_bohr: ArrayLike) -> float | NDArray:
        """5Z neutral relative to its own molecular total at R=2.1323 bohr."""
        radius = _radius(radius_bohr)
        return _scalar_or_array(chebval((radius - 2.2) / 0.3, self.neutral_coefficients))

    def er(self, radius_bohr: ArrayLike) -> float | NDArray:
        """Native incident resonance position (Hartree)."""
        radius = _radius(radius_bohr)
        return _scalar_or_array(_cubic(self.knots_bohr, self.er_coefficients, radius))

    def full_gamma(self, radius_bohr: ArrayLike) -> float | NDArray:
        """Positive full width (Hartree), not the half width."""
        radius = _radius(radius_bohr)
        with np.errstate(over="ignore", under="ignore"):
            width = np.exp(_cubic(self.knots_bohr, self.log_gamma_coefficients, radius))
        if not np.all(np.isfinite(width)) or np.any(width <= 0):
            raise ValueError("Full width is not finite and positive")
        return _scalar_or_array(width)

    def incident_pole(self, radius_bohr: ArrayLike) -> complex | NDArray:
        """Incident pole z=Er-i*full_Gamma/2, in Hartree."""
        return self.er(radius_bohr) - 0.5j * self.full_gamma(radius_bohr)

    def vres_hybrid(self, radius_bohr: ArrayLike) -> float | NDArray:
        """Mixed-method V0_5Z+Er_native; NOT a correlated anion total."""
        return self.v0(radius_bohr) + self.er(radius_bohr)


def _read_json_members(path: str | Path, sha256: str, names: tuple | list) -> dict:
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != sha256:
        raise ValueError(f"Capsule SHA256 mismatch: expected {sha256}")
    # Parse only authenticated bytes and named regular files, without extractall.
    result = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        members = archive.getmembers()
        for name in names:
            matches = [m for m in members if m.name == name]
            if len(matches) != 1 or not matches[0].isfile():
                raise ValueError(f"Capsule needs exactly one regular member: {name}")
            with archive.extractfile(matches[0]) as stream:
                result[name] = json.load(stream)
    return result


def _from_payloads(payloads: dict) -> LocalRepresentation:
    acceptance = payloads["representation.json"]
    model = payloads["inputs/frozen-c2-model.json"]
    neutral = payloads["inputs/accepted-neutral-result.json"]
    if (
        acceptance["accepted_for_scoped_local_use"] is not True
        or acceptance["no_extrapolation"] is not True
        or acceptance["domain_bohr"] != list(DOMAIN_BOHR)
        or model["accepted_representation"] is not False
        or model["family"] != "CubicSpline_not_a_knot_Er_logFullGamma"
        or neutral["status"] != "pass"
        or neutral["reference_R_bohr"] != 2.1323
    ):
        raise ValueError("Expected scoped acceptance with the unchanged historical candidate")
    er = model["models"]["Er"]
    gamma = model["models"]["log_full_Gamma"]
    if er["knots_bohr"] != gamma["knots_bohr"] or er["knots_bohr"] != model["development_R_bohr"]:
        raise ValueError("Frozen resonance pieces must share the development knots")
    return LocalRepresentation(
        neutral["fits"]["5z"]["coefficients_hartree"],
        er["knots_bohr"],
        er["coefficients_descending_local_power"],
        gamma["coefficients_descending_local_power"],
    )


def load_representation(capsule: str | Path) -> LocalRepresentation:
    """Authenticate the pinned external capsule, then load saved coefficients."""
    return _from_payloads(_read_json_members(capsule, DELIVERY_SHA256, _CURVE_MEMBERS))


def load_plot_evidence(capsule: str | Path) -> dict:
    """Read only accepted local curves/pairs/scores, not protected S or pole data."""
    return _read_json_members(capsule, DELIVERY_SHA256, _EVIDENCE_MEMBERS)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--r-bohr", type=float, nargs="+", help="Real radii in [1.9, 2.5] bohr")
    parser.add_argument("--plot-dir", type=Path, help="New external directory for saved-only plots")
    args = parser.parse_args()
    curve = load_representation(args.capsule)
    if args.r_bohr:
        radius = np.asarray(args.r_bohr)
        print(
            json.dumps(
                {
                    "R_bohr": radius.tolist(),
                    "energy_unit": "Hartree",
                    "neutral_zero": "Own molecular total at R=2.1323 bohr, not minimum/fragments",
                    "V0_5Z": curve.v0(radius).tolist(),
                    "Er_native": curve.er(radius).tolist(),
                    "full_Gamma": curve.full_gamma(radius).tolist(),
                    "Vres_hybrid_5Z_plus_native_NOT_correlated_anion": curve.vres_hybrid(
                        radius
                    ).tolist(),
                    "scope": "Local finite model; no uniform/electronic/global accuracy claim",
                },
                indent=2,
            )
        )
    if args.plot_dir:
        from projects.ukrmol_co.local_representation_plot import plot_saved_curves

        for path in plot_saved_curves(curve, load_plot_evidence(args.capsule), args.plot_dir):
            print(path)


if __name__ == "__main__":
    main()
