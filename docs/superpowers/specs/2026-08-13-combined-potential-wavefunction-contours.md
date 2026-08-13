# Combined potential + wavefunction contours — Design Spec

**Date:** 2026-08-13
**Status:** Approved design (pending review) — extends the `qscat.viz` contour overlay
**Lifecycle:** visualisation add-on on top of `qscat.viz` (no new physics).

## Goal

Draw **both** overlays on one domain-coloured image: the solid-white `|ψ|`
contours (existing) *plus* a second, **dotted** set of **potential** contours
placed at **physically relevant energies**. The result reads directly — the
dotted lines are the classical turning surfaces "where the wavefunction is
*allowed* at energy E", the solid lines are "where it *actually is*" — exposing
turning points, tunnelling into the forbidden region, and resonance trapping.

## Physical principle for the levels

A state's total energy is `E_tot = ε_v + E_collision`. The equipotential
`V(r,R) = E_tot` is the classical turning surface. So the "relevant" potential
levels are the total PES sampled at the energies in play: the vibrational
thresholds `ε_v` and/or `ε_{v_init} + E` over the collision-energy range being
shown. Contouring the **total potential surface** (Hartree) at those levels is
what makes the overlay physically meaningful; the potential field the caller
supplies must be in the same energy units as the levels.

## API — dedicated potential overlay (chosen)

`plot_wavefunction_2d` keeps its `|ψ|` contour params and gains a parallel,
independent potential overlay drawn when `potential_levels` is given:

```python
plot_wavefunction_2d(
    ..., contours=True,                 # solid-white |psi| overlay (existing)
    # --- new dedicated potential overlay ---
    potential=<nodal array | callable>, # same-grid field (project_values) or analytic V(r,R)
    potential_levels: Sequence[float] | Literal["auto"] | None = None,
    potential_style: str = ":",         # dotted
    potential_color: str = "0.75",      # subtle grey (overridable)
    potential_alpha: float = 0.7,
    potential_linewidth: float = 0.6,
    potential_labels: bool = True,      # inline clabel of the energy on each line
    potential_label_fmt: str = "%.3f",
    # for potential_levels="auto":
    eps: npt.ArrayLike | None = None,   # vibrational energies
    v_init: int = 0,
    energies: npt.ArrayLike | None = None,  # collision energies -> E_tot levels
)
```

- The potential overlay draws iff `potential is not None and potential_levels is
  not None`. It is independent of `contours`, so any of {|ψ| only, potential
  only, both} is expressible.
- `potential_levels` is an explicit list, or `"auto"` (uses the helper below with
  `eps`/`v_init`/`energies`), clipped to the shown potential range.
- Reuses the existing same-grid `potential` handling (nodal array via
  `EquidistantProjector.project_values`, or a callable) — unchanged.
- Relationship to the existing `contour_field="potential"`: that stays for a
  potential-only *primary* figure; the new params are the way to COMBINE. Both
  route through one internal `_draw_potential_contours` helper (no duplicate
  logic). If this proves redundant in practice we can deprecate
  `contour_field="potential"` later (pre-1.0, cheap).

## Level-selection helper

`qscat.viz.energy_contour_levels(*, eps=None, v_init=0, energies=None,
include_thresholds=True, e_range=None, max_levels=12, min_spacing=None) ->
list[float]`:

- Collects `{ε_v}` (if `include_thresholds`) ∪ `{ε_{v_init} + E for E in
  energies}`.
- Clips to `e_range` (defaults to the shown potential's finite range when called
  from the plot), thins by `min_spacing`, caps at `max_levels` (keep the most
  spread-out subset), returns sorted levels.
- Pure array logic (model-independent): `eps`/`energies` are plain arrays, so it
  lives in `qscat.viz` and works for any molecule. `ScatteringProblem.eps` feeds
  it directly.

## Styling & readability

- Potential: `linestyle=":"`, thin, subtle grey, `alpha≈0.7` — visually distinct
  from the solid-white `|ψ|` lines.
- `potential_labels=True` → `ax.clabel(cset, fmt=potential_label_fmt, fontsize=…)`
  prints the energy (e.g. `0.100`) inline on each dotted line, so the turning
  surface's energy is self-documenting (no separate legend).
- Contour orientation shares the imshow-inverted axes (same handling as the
  existing overlay) so both sets align with the colour.

## Validation / tests

- **Visual**: re-render the N₂ Ψ₊ demo with `contours=True` + a total-PES
  `potential` at `potential_levels="auto"` (from the problem's `eps` and the
  shown collision energies) — confirm dotted energy-labelled turning lines
  overlaid on the solid `|ψ|` contours, aligned to the colour.
- **Unit** (`test_viz.py`):
  - `energy_contour_levels`: thresholds + `ε_{v_init}+E` selected; `e_range`
    clip, `min_spacing` thinning, `max_levels` cap; sorted output.
  - Both overlays present: with `contours=True` and `potential`/`potential_levels`
    set, TWO `QuadContourSet`s are added; the potential one has
    `linestyle=":"`, the requested colour/alpha; `clabel` texts exist when
    `potential_labels=True`.
  - potential overlay independent of `contours` (potential-only draws one set).

## Deliverables

- `energy_contour_levels` in `qscat.viz` (+ export) and its tests.
- `plot_wavefunction_2d` extended with the potential-overlay params + a
  `_draw_potential_contours` helper shared with `contour_field="potential"`.
- Updated N₂ demo figure (both overlays) + CHANGELOG.

## Verification

`uv run pytest libs/qscat/tests/test_viz.py` green; ruff + mypy clean; the demo
shows dotted, energy-labelled turning contours over the solid `|ψ|` contours,
both aligned to the domain-coloured field.
