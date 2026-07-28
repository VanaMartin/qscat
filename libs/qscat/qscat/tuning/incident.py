"""Incident-state + test-function placement for the FEM-DVR-ECS
discretisation tuner (design spec: "the TW analysis").

The grid must CONTAIN and RESOLVE not just the potential's eigenstates but
the incident state and the test functions (the outgoing-channel
projections / flux-extraction points) -- these can dominate the required
extent and resolution (see
`docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md`'s
"Incident state + test-function placement" section):

- **TI route**: the incident is the channel function (Bessel/Coulomb) and
  the test functions are the exit-channel projections, both set directly
  by the energy range and already covered by the channel-representation
  probes (`qscat.tuning.probe_channel_representation`). There is no
  position/width to place -- `IncidentSpec` below is a no-op for this
  route beyond that existing probe.
- **TD route (Tannor-Weeks)**: the incident is a GAUSSIAN WAVEPACKET
  `g(r) = exp(-(r-position)^2 / (2 sigma^2)) exp(i impulse r)` (the
  `qscat.core.wavepacket.gaussian_coeffs` envelope, up to its own
  normalization prefactor) placed far out, with an `observation` boundary
  where the outgoing test function is read. `IncidentSpec` models this
  case: `required_extent` says how far the real region must reach to
  contain it, and `tw_analysis` auto-places it for a target energy range.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["IncidentSpec", "required_extent", "tw_analysis"]

# How many sigma of Gaussian tail the real region must contain beyond the
# wavepacket's centre `position`. At 5 sigma, exp(-5**2/2) ~ 4e-6 of the
# peak amplitude -- negligible, so the packet (and its tail) is genuinely
# INSIDE the real region rather than clipped at the boundary.
_TAIL_SIGMAS = 5.0

# --- TW analysis (auto-placing the TD Gaussian for an energy range) ---

# Safety factor (< 1) applied to the marginal energy-spread bound `sigma <=
# impulse / (E_max - E_min)` derived in `tw_analysis`'s docstring, so the
# realized spread comfortably exceeds the required half-range rather than
# just touching it.
_SIGMA_SAFETY_FACTOR = 0.6

# How many sigma of standoff to place the wavepacket beyond the radius
# where the model's electron-molecule interaction has decayed away -- gives
# the packet room to start as a genuinely free (non-interacting) Gaussian.
_STANDOFF_SIGMAS = 6.0

# Where `v_int(r, R)` counts as "decayed away", relative to its r=0 value.
_INTERACTION_DECAY_RTOL = 1e-3

# Domain + sample count for `_interaction_extent`'s scan -- generous enough
# to bracket the Gaussian-in-r interaction cutoff of every registered model
# (N2/NO/F2's `alpha_c` in [0.4, 3.0], decaying well inside 10 bohr).
_INTERACTION_SCAN_MAX = 30.0
_INTERACTION_SCAN_N = 3000


@dataclass(frozen=True)
class IncidentSpec:
    """TD Gaussian-wavepacket placement + observation boundary.

    `position` (bohr): the wavepacket centre `r0` at `t=0`.
    `impulse` (bohr^-1): the mean momentum `p0`; negative launches the
        packet INWARD (toward the interaction region), matching the
        `n2_2d_td_cross_section` convention (`wp_in`'s `p0 < 0` -- see
        `projects/n2_2d_td_cross_section/convergence.py`).
    `sigma` (bohr): the Gaussian width.
    `observation`: the flux-extraction / test-function boundary (bohr);
        `None` if the caller has not set one (`required_extent` then falls
        back to just the wavepacket's own tail).

    For the TI route this dataclass is not used -- see the module
    docstring.
    """

    position: float
    impulse: float
    sigma: float
    observation: float | None = None

    def required_extent(self) -> float:
        """How far the REAL region must reach to CONTAIN this wavepacket
        and its observation boundary: `max(position + N*sigma,
        observation)`, `N = _TAIL_SIGMAS`.

        A METHOD (not a free function) so `propose_grid`'s existing
        `getattr(incident, "required_extent", lambda: 0.0)()` call works
        against an `IncidentSpec` unchanged.
        """
        tail = self.position + _TAIL_SIGMAS * self.sigma
        observation = self.observation if self.observation is not None else 0.0
        return max(tail, observation)

    def incident_energy(self) -> float:
        """The wavepacket's mean kinetic energy, `impulse**2 / 2` (mass 1,
        the electronic-coordinate convention every `IncidentSpec` is
        defined against).

        `propose_grid` reads this to make sure the RESOLUTION (not just
        the extent) of the a-priori mesh covers the incident's local
        wavenumber -- see its docstring for why `required_extent` alone
        is not enough whenever a hand-built `IncidentSpec` implies an
        energy above `energy_range`'s own `E_max`.
        """
        return self.impulse**2 / 2.0


def required_extent(spec: IncidentSpec) -> float:
    """Free-function form of `IncidentSpec.required_extent`, for callers
    who prefer a function over a bound method."""
    return spec.required_extent()


def _interaction_extent(model: ResonanceModel) -> float:
    """The electronic radius `r` beyond which `model.v_int(r, R)` has
    decayed to `_INTERACTION_DECAY_RTOL` of its `r=0` value.

    Uses only `v_int` (part of the `ResonanceModel` protocol) at an
    ARBITRARY fixed `R`: for the shared Morse+sigmoid+Gaussian-in-r form,
    `v_int(r, R) = -lambda(R) * exp(-alpha_c * r**2)`, so the ratio
    `|v_int(r, R)| / |v_int(0, R)| = exp(-alpha_c * r**2)` does not depend
    on `R` at all -- any fixed `R` gives the same extent. That keeps this
    estimate genuinely model-independent (works for any `ResonanceModel`,
    not just `DiatomicResonanceModel`), as long as `v_int` decays
    monotonically in `r` away from a maximum at `r=0` (true of every
    registered model).
    """
    r = np.linspace(0.0, _INTERACTION_SCAN_MAX, _INTERACTION_SCAN_N)
    v = np.abs(np.asarray(model.v_int(r, 1.0), dtype=np.complex128))
    v0 = float(v[0])
    if v0 <= 0.0:
        return 0.0
    below = np.nonzero(v / v0 < _INTERACTION_DECAY_RTOL)[0]
    if below.size == 0:
        return float(r[-1])
    return float(r[int(below[0])])


def tw_analysis(model: ResonanceModel, energy_range: tuple[float, float]) -> IncidentSpec:
    """Auto-place the TD Gaussian wavepacket for `energy_range = (E_min,
    E_max)` -- the design spec's "TW analysis", BEST-EFFORT and not
    independently calibrated against the eMoScat decks: Task 8 calibrated
    the mesh's de-Broglie phase constant `C` (`qscat.tuning.mesh.
    _PHASE_COEFF`), not this placement heuristic, which remains a documented
    follow-on (the design spec lists the TW auto-tune as best-effort).

    Physics: `g(r) = exp(-(r-position)^2 / (2 sigma^2)) exp(i impulse r)`
    has mean energy `impulse**2 / 2` and, near the mean, an energy spread
    `delta_E ~ impulse * delta_p` with `delta_p ~ 1 / (2 sigma)` -- i.e.
    `delta_E ~ impulse / (2 sigma)`. Inverting `[E_min, E_max] ->
    (impulse, sigma, position)`:

    - `impulse`: `p = sqrt(2 * E_centre)`, `E_centre = (E_min+E_max)/2`, so
      the packet's mean energy sits at the range's centre; launched
      INWARD (`impulse = -p`), matching `n2_2d_td_cross_section`'s
      `wp_in` convention.
    - `sigma`: the spread must COVER the range, `delta_E >= (E_max -
      E_min) / 2 <=> sigma <= p / (E_max - E_min)`; `sigma` is set to
      `_SIGMA_SAFETY_FACTOR` times that marginal bound so the realized
      spread comfortably brackets `[E_min, E_max]` rather than just
      touching its edges (a smaller `sigma` only widens the momentum/
      energy spread further, so this is the conservative direction).
    - `position`: `_interaction_extent(model) + _STANDOFF_SIGMAS * sigma`
      -- far enough out that the packet starts genuinely outside the
      interaction region, with room to spare.
    - `observation`: set equal to `position` (a sane symmetric default --
      the outgoing test function is read at the same radius the packet
      launches from; construct `IncidentSpec` directly for a different
      boundary).
    """
    e_min, e_max = energy_range
    if e_max <= e_min:
        raise ValueError(f"energy_range must have E_max > E_min, got {energy_range}")

    e_centre = 0.5 * (e_min + e_max)
    p = math.sqrt(2.0 * e_centre)

    sigma = _SIGMA_SAFETY_FACTOR * p / (e_max - e_min)

    position = _interaction_extent(model) + _STANDOFF_SIGMAS * sigma

    return IncidentSpec(position=position, impulse=-p, sigma=sigma, observation=position)
