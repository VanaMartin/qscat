# Wavefunction contour overlay for `qscat.viz` — Design Spec

**Date:** 2026-08-02
**Status:** Approved design (pending review) — extends `qscat.viz.plot_wavefunction_2d`
**Lifecycle:** visualisation add-on on top of the validated `qscat.viz` projector.

## Goal

Overlay contour lines on the domain-coloured 2-D wavefunction image — **thin
white lines at 0.6 opacity** — with the contoured quantity selectable between the
wavefunction **magnitude** and the **potential**, mirroring eMoScat's
`display_wf.py` (`ax.contour(X, Y, |Z|², colors="w", alpha=0.5)`). Off by default
so existing callers are unchanged; carried into the future `animate()`.

## API (additive to `plot_wavefunction_2d`)

```python
plot_wavefunction_2d(
    projector, state, *, mag, path=None, inverse=False, title=None,
    xlabel=..., ylabel=..., ax=None,
    # --- new, all optional ---
    contours: bool | int | Sequence[float] = False,        # off | default | N levels | explicit
    contour_field: Literal["magnitude", "potential"] = "magnitude",
    potential: Callable[[np.ndarray, np.ndarray], np.ndarray] | np.ndarray | None = None,
    contour_color: str = "white",
    contour_alpha: float = 0.6,
    contour_linewidth: float = 0.6,   # "thin"
) -> Any
```

- `contour_field="magnitude"` — contour `|field|` from the SAME projected state
  the colouring uses (no second projection).
- `contour_field="potential"` — contour a potential supplied via the optional
  `potential=` argument as a real field **on the same discretisation (tensor)
  grid** as the state, shape `(n0_nodes, n1_nodes)` or flat. It is projected
  through the SAME `EquidistantProjector`, so contours land on the identical
  sampling grid as the colour. Selecting `"potential"` without `potential=`
  raises `ValueError`.

  Convention note: a wavefunction `state` is √w-scaled (so `project` divides by
  √w to recover ψ), whereas a potential is nodal *values* `V(x_i)` whose
  interpolant is `Σ V(x_i) L_i(x)` (no √w). So the projector gains a
  `project_values(field)` method (= `project(√w ⊙ field)`) that interpolates a
  nodal-value field on the same grid; the potential contour uses it.

## Levels — derived from `mag` (the chosen scheme)

- `contours=True`:
  - **magnitude**: levels tied to the colour scale — `k · mag / 5` for
    `k = 1 … 20` (i.e. up to `4·mag`, eMoScat's `wlevels` range, but on the
    magnitude scale since we contour `|ψ|`, not `|ψ|²`). Contours and colours
    then share one scale, and — crucially for animation — the levels are FIXED
    across frames (no auto-ranging "breathing").
  - **potential**: `mag` is a wavefunction scale, so it does not define potential
    levels; default to 10 levels spanning the finite range of `V` on the grid.
    (Documented asymmetry; `contours=int`/`Sequence` override.)
- `contours=<int N>`: N levels (magnitude: `k·mag/5`, `k=1…N`; potential: N over range).
- `contours=<Sequence>`: exactly those levels.

## Two implementation details (handled, flagged for honesty)

1. **Orientation.** The image is drawn with `imshow(origin="upper")` and a
   flipped-`r` extent; `ax.contour` works in data coordinates. The overlay must
   use coordinate arrays / `origin` consistent with the image so lines land on
   the features (the standard imshow+contour alignment gotcha). Concretely:
   `ax.contour(projector.axis1, projector.axis0, Z, levels=..., origin="upper",
   extent=<same as imshow>)`, `Z` shape `(n0, n1)`.
2. **Animation.** A `QuadContourSet` cannot be `set_data`'d — each frame removes
   the previous contour set and draws a new one (a standard matplotlib pattern;
   slightly slower than pure-`imshow` blitting). The `animate()` spec will note
   this; here we just make `plot_wavefunction_2d` return/track the contour set so
   an animation loop can remove it.

## Styling

`colors=contour_color` (default white), `alpha=contour_alpha` (0.6),
`linewidths=contour_linewidth` (0.6). Applies to whichever `contour_field` is
chosen. All overridable (e.g. grey potential contours).

## Validation

- **Visual**: re-render the N₂ Ψ₊ demo with `contours=True` (magnitude), and a
  second with `contour_field="potential"` (the model PEC) — confirm thin white
  lines at 0.6 α overlay correctly (alignment check).
- **Test** (`test_viz.py`): with `contours=True`, a `QuadContourSet` is added to
  the axes with the requested `color`/`alpha`/`linewidth`; magnitude levels equal
  `k·mag/5`; `contour_field="potential"` without `potential=` raises `ValueError`;
  `contours=False` adds no contour set (existing behaviour unchanged).

## Deliverables

- `plot_wavefunction_2d` extended with the params above + a small
  `_contour_levels(contours, contour_field, mag, Z)` helper.
- Tests + an updated N₂ demo figure with contours.
- CHANGELOG entry.

## Verification

`uv run pytest libs/qscat/tests/test_viz.py` green; ruff + mypy clean; the
demo renders thin white 0.6-α contours aligned to the coloured field.
