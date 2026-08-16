# Gallery Enabling Features Implementation Plan (Plan B)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the two capabilities the showcase gallery needs before any figure can be produced — region-split domain colouring in `qscat.viz`, and a `reference:` block in the `qscat-run` config so a published dataset can be overlaid on a computed cross section.

**Architecture:** Both are small, self-contained additions to existing modules. `qscat.viz.coloring` gains a pure-numpy helper that computes a per-point magnitude scale from region boundaries, and `mag` becomes array-accepting throughout the viz call chain. `qscat_run` gains a `reference.py` loader that reads a data file **by path** — it never imports `validation`, so the enforced layering holds.

**Tech Stack:** numpy, matplotlib (viz), PyYAML + click (qscat-run config), pytest.

**Source spec:** `docs/superpowers/specs/2026-08-16-documentation-showcase-design.md`, sections B1 and B2.

## Scope: why this is its own plan

The spec described one "Plan B" covering both these features *and* the gallery itself (eight pages, production runs on `sadaharu`, figure dedup, README rework). Splitting off the two features is a deliberate deviation, for two reasons:

1. These are code with tests and merge on their own. The gallery is content production that consumes them.
2. A gallery plan written *now* would be full of placeholders. Page copy depends on what the figures actually show, and the cost tables depend on measured wall-clock from real runs. Neither is knowable until these features exist and a trial run has happened.

The gallery becomes its own plan, written after this merges and a first run is on disk.

## Global Constraints

- Atomic units throughout (Hartree, bohr). Do not introduce unit conversions outside `qscat.units`.
- The public API is exactly each submodule's `__all__` (ADR 0004). **Any new public name must be added to its `__all__` AND to its page under `docs/api/`** — `tests/test_api_docs_coverage.py` fails otherwise. Write autosummary entries as **bare unqualified names** under a `.. currentmodule::`; the gate does not match fully-qualified dotted entries.
- The docs must still build under the CI command: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`.
- `qscat_run` **must not import `validation`** — enforced by `apps/qscat-run/tests/test_no_validation_import.py`.
- `qscat.viz.coloring` must stay pure numpy (no matplotlib import), so it is testable without the `plot` extra.
- Existing callers passing a scalar `mag` must keep working unchanged. This is a widening, not a breaking change.
- Run tests with `uv run --no-sync pytest`.
- Preserve the repo's commit trailers on every commit:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS`

## Verified facts this plan is built on

Measured on 2026-08-16 against the working tree. Do not re-derive; do verify anything you change.

1. `qscat/viz/coloring.py` exposes `complex_to_hsv(z, mag=1.0, *, inverse=False)`, `complex_to_rgb(...)`, `hsv_to_rgb(hsv)`. `__all__ = ["complex_to_hsv", "complex_to_rgb", "hsv_to_rgb"]`. It imports only numpy.
2. The brightness path is `r = np.abs(c) / mag`; then (non-inverse) `saturation = 1/max(1, r)`, `value = min(1, r)`. So a point whose magnitude is `1e-6 * mag` gets `value = 1e-6` — black. That is the defect being fixed.
3. `mag` threads through as a **scalar float** in three more places: `plot_wavefunction_2d(..., mag: float, ...)` (`plot.py:32`), `WavefunctionArtist.__init__(..., mag: float, ...)` (`artist.py:108`), and `animate_wavefunction(..., mag: float, ...)` (`animate.py:76`).
4. `artist.py:34` `contour_levels(n, contour_field, mag: float, z)` derives magnitude contours as `k * mag / 5.0` for `k` in `1..n`. Those levels **must be scalars** — this is the one place an array `mag` cannot simply flow through.
5. `validation/n2/data/CSVE.V00.J00` is whitespace-delimited, **400 rows × 32 columns**: column 0 is energy in Hartree (strictly increasing), columns 1..31 are σ in bohr² for `v' = 0..30` (`v'=0` elastic). Loaded elsewhere with a plain `np.loadtxt`.
6. `ExperimentConfig` (`apps/qscat-run/qscat_run/config.py`) is a frozen dataclass with fields `molecule, methods, observables, output_dir, energies, grid, v_init, td, artifacts, backend`. `load_config` parses YAML then constructs it; `validate_config` raises `ConfigError` (a `click.ClickException`) with an actionable message for the first problem found.
7. `artifacts.py` writes cross sections through `_write_cross_section_csv(path, energies, series)`, `_write_cross_section_npz(...)`, and `_write_cross_section_png(...)`, where `series` is `dict[str, ndarray]` and **every array is indexed by the same `energies` vector**.

## A spec correction, made deliberately

The spec said reference series would be "written as extra columns in `cross_section.csv` under a `ref:` key prefix". **That is wrong and this plan does not do it.**

Fact 7 is why: every column in `cross_section.csv` shares one energy axis. A reference dataset has its own sampling — Houfek's has 400 energies (fact 5), while a config typically requests a few dozen. Putting reference values in that table would require interpolating published data onto our grid and presenting the result as the reference. That misrepresents someone else's measurement, in the one figure whose entire purpose is honest external comparison.

Instead: the reference keeps its own energy axis, is written to its **own** `reference.csv` and its own arrays in the npz, and is overlaid on the PNG plotted against its own energies.

---

## File Structure

**Create:**
- `apps/qscat-run/qscat_run/reference.py` — reference-dataset loading, keyed by named format.
- `apps/qscat-run/tests/test_reference.py` — loader + config-integration tests.
- `apps/qscat-run/examples/n2-ve-vs-houfek.yaml` — a config exercising `reference:`.

**Modify:**
- `libs/qscat/qscat/viz/coloring.py` — array-accepting `mag`; new `region_magnitudes`.
- `libs/qscat/qscat/viz/__init__.py` — export `region_magnitudes`.
- `libs/qscat/qscat/viz/artist.py` — accept array `mag`; scalar reference for contour levels.
- `libs/qscat/qscat/viz/plot.py`, `animate.py` — widen the `mag` type.
- `libs/qscat/tests/test_viz_coloring.py` — new tests (create if absent).
- `docs/api/viz.md` — add `region_magnitudes` to the autosummary.
- `apps/qscat-run/qscat_run/config.py` — `ReferenceSpec`, `_load_reference`, `ExperimentConfig.reference`, validation.
- `apps/qscat-run/qscat_run/artifacts.py` — write and overlay reference series.
- `apps/qscat-run/README.md` — document `reference:` in the observables/config tables.

---

### Task 1: Region-split magnitude scaling in `qscat.viz`

**Files:**
- Modify: `libs/qscat/qscat/viz/coloring.py`, `libs/qscat/qscat/viz/__init__.py`
- Test: `libs/qscat/tests/test_viz_coloring.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `region_magnitudes(magnitude, *, axis, boundaries, percentile=99.5, floor=1e-300) -> NDArray[np.float64]` returning an array the same shape as `magnitude`; and `complex_to_hsv`/`complex_to_rgb` accepting `mag: float | NDArray[np.float64]`. Task 2 widens the callers that pass `mag` onward.

- [ ] **Step 1: Write the failing tests**

Create `libs/qscat/tests/test_viz_coloring.py` (if the file exists, append these):

```python
"""Region-split magnitude scaling for domain colouring.

The defect being fixed: `complex_to_hsv` normalises the whole field by one
scalar `mag`, so in a time-dependent wavefunction -- where the incident packet
outweighs the resonant and outgoing amplitude by orders of magnitude -- every
interesting feature renders black. See docs/physics/ for the physics; this
module is pure numpy.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.viz import complex_to_hsv, complex_to_rgb, region_magnitudes


def _two_region_field() -> np.ndarray:
    """Left half amplitude 1.0, right half amplitude 1e-6 -- the incident-packet
    vs outgoing-tail contrast that motivates the whole feature."""
    z = np.ones((8, 10), dtype=np.complex128)
    z[:, 5:] = 1e-6
    return z


def test_scalar_mag_is_unchanged() -> None:
    # Regression guard: the widening must not alter existing scalar behaviour.
    z = np.array([[1 + 1j, 0.5], [0.0, 2.0]], dtype=np.complex128)
    hsv = complex_to_hsv(z, 2.0)
    assert hsv.shape == z.shape + (3,)
    assert np.all((hsv >= 0.0) & (hsv <= 1.0))
    # |z|=2 at mag=2 saturates value to exactly 1.0
    assert hsv[1, 1, 2] == pytest.approx(1.0)


def test_scalar_mag_leaves_the_weak_region_black() -> None:
    # This is the BUG, asserted so the fix has something to improve on.
    z = _two_region_field()
    hsv = complex_to_hsv(z, 1.0)
    assert hsv[0, 0, 2] == pytest.approx(1.0)      # strong region: full brightness
    assert hsv[0, 9, 2] < 1e-5                      # weak region: invisible


def test_region_magnitudes_rescales_each_region_independently() -> None:
    z = _two_region_field()
    mag = region_magnitudes(np.abs(z), axis=1, boundaries=[5])
    assert mag.shape == z.shape
    hsv = complex_to_hsv(z, mag)
    # Both regions now reach full brightness -- the point of the feature.
    assert hsv[0, 0, 2] == pytest.approx(1.0)
    assert hsv[0, 9, 2] == pytest.approx(1.0)


def test_region_magnitudes_uses_a_percentile_not_the_max() -> None:
    # One hot pixel must not re-flatten its region. With max-scaling the bulk
    # would sit at 1/1000; with a percentile it stays visible.
    m = np.ones((1, 100))
    m[0, 0] = 1000.0
    mag = region_magnitudes(m, axis=1, boundaries=[], percentile=90.0)
    assert mag[0, 50] == pytest.approx(1.0, rel=1e-6)


def test_region_magnitudes_is_positive_even_for_an_all_zero_region() -> None:
    # An empty region must not produce a 0 scale and a divide-by-zero downstream.
    m = np.zeros((4, 6))
    mag = region_magnitudes(m, axis=1, boundaries=[3])
    assert np.all(mag > 0.0)
    assert np.all(np.isfinite(complex_to_rgb(m.astype(np.complex128), mag)))


def test_region_magnitudes_rejects_bad_boundaries() -> None:
    m = np.ones((4, 6))
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[0])        # empty leading region
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[6])        # empty trailing region
    with pytest.raises(ValueError):
        region_magnitudes(m, axis=1, boundaries=[4, 2])     # not increasing


def test_complex_to_hsv_rejects_a_non_broadcastable_mag() -> None:
    z = np.ones((4, 6), dtype=np.complex128)
    with pytest.raises(ValueError):
        complex_to_hsv(z, np.ones((3, 3)))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest libs/qscat/tests/test_viz_coloring.py -v`

Expected: FAIL — `ImportError`/`AttributeError` on `region_magnitudes`, since it does not exist yet. `test_scalar_mag_is_unchanged` and `test_scalar_mag_leaves_the_weak_region_black` may already pass; that is fine, they are regression guards.

- [ ] **Step 3: Implement `region_magnitudes` and widen `mag`**

In `libs/qscat/qscat/viz/coloring.py`, add to `__all__`: `"region_magnitudes"`.

Add the function:

```python
def region_magnitudes(
    magnitude: npt.ArrayLike,
    *,
    axis: int,
    boundaries: Sequence[int],
    percentile: float = 99.5,
    floor: float = 1e-300,
) -> npt.NDArray[np.float64]:
    """Per-point brightness scale that normalises each region separately.

    A single global ``mag`` ties the whole field to its largest feature, so a
    region whose amplitude is orders of magnitude smaller renders black. This
    splits ``axis`` at ``boundaries`` and gives every point the scale of its own
    region, so each region is visible on its own terms.

    The scale is a percentile rather than the maximum, so one outlying point
    cannot flatten the region it sits in.

    Parameters
    ----------
    magnitude : array_like
        Non-negative field, e.g. ``np.abs(psi)``.
    axis : int
        Axis the split runs along.
    boundaries : sequence of int
        Strictly increasing split indices along ``axis``, each in
        ``1 .. n-1``. An empty sequence means one region (the whole field).
    percentile : float
        Percentile of each region's magnitudes mapped to full brightness.
    floor : float
        Lower clamp, so an all-zero region yields a positive scale rather than
        a divide-by-zero downstream.

    Returns
    -------
    ndarray
        Same shape as ``magnitude``; every point carries its region's scale.
        Pass it straight to ``complex_to_rgb`` as ``mag``.
    """
    m = np.asarray(magnitude, dtype=np.float64)
    n = m.shape[axis]
    bounds = [int(b) for b in boundaries]
    if any(b < 1 or b > n - 1 for b in bounds):
        raise ValueError(
            f"boundaries must lie in 1..{n - 1} along axis {axis} (length {n}); got {bounds}"
        )
    if any(b >= c for b, c in zip(bounds, bounds[1:], strict=True)):
        raise ValueError(f"boundaries must be strictly increasing; got {bounds}")

    edges = [0, *bounds, n]
    out = np.empty_like(m)
    for lo, hi in zip(edges, edges[1:], strict=True):
        sl: list[slice] = [slice(None)] * m.ndim
        sl[axis] = slice(lo, hi)
        block = m[tuple(sl)]
        scale = float(np.percentile(block, percentile)) if block.size else 0.0
        out[tuple(sl)] = max(scale, floor)
    return out
```

Add `from collections.abc import Sequence` to the imports.

Then widen the two entry points. In both `complex_to_hsv` and `complex_to_rgb`, change the annotation `mag: float = 1.0` to `mag: float | npt.NDArray[np.float64] = 1.0`, and in `complex_to_hsv` replace `r = np.abs(c) / mag` with:

```python
    scale = np.asarray(mag, dtype=np.float64)
    if scale.ndim and scale.shape != c.shape:
        try:
            scale = np.broadcast_to(scale, c.shape)
        except ValueError as exc:
            raise ValueError(
                f"mag of shape {scale.shape} is not broadcastable to the field shape {c.shape}"
            ) from exc
    r = np.abs(c) / scale
```

Update the module docstring's description of `Mag` to say it may be a scalar or a per-point array, and mention `region_magnitudes`. Leave the existing print-mode TODO alone.

- [ ] **Step 4: Export the new name**

In `libs/qscat/qscat/viz/__init__.py`, add `region_magnitudes` to the `from .coloring import ...` line and to `__all__`, keeping the existing ordering style.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_viz_coloring.py -v`

Expected: all PASS.

Then run the existing viz suite to confirm nothing regressed:
`uv run --no-sync pytest libs/qscat/tests -k viz -v`

Expected: no failures.

- [ ] **Step 6: Document the new public name**

Add `region_magnitudes` to the `autosummary` block in `docs/api/viz.md`, as a **bare name** (the list sits under `.. currentmodule:: qscat.viz`). Place it next to the other colouring helpers.

Run: `uv run --no-sync pytest tests/test_api_docs_coverage.py -q`
Expected: 13 passed. This test fails if a new public name has no docs entry — that is its job.

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: `build succeeded.` with no warnings.

- [ ] **Step 7: Commit**

```bash
git add libs/qscat/qscat/viz/coloring.py libs/qscat/qscat/viz/__init__.py \
        libs/qscat/tests/test_viz_coloring.py docs/api/viz.md
git commit -m "feat(viz): per-region magnitude scaling for domain colouring

complex_to_hsv normalised the whole field by one scalar mag, so a region
orders of magnitude weaker than the brightest feature rendered black. In a
time-dependent wavefunction that is exactly the interesting part: the
resonant and outgoing amplitude next to the incident packet.

region_magnitudes splits an axis and gives each region its own scale, taken
at a percentile so a single hot pixel cannot flatten its region. mag is now
scalar-or-array everywhere; scalar callers are unaffected.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

### Task 2: Thread array `mag` through the plotting call chain

**Files:**
- Modify: `libs/qscat/qscat/viz/artist.py`, `libs/qscat/qscat/viz/plot.py`, `libs/qscat/qscat/viz/animate.py`
- Test: `libs/qscat/tests/test_viz_coloring.py` (append)

**Interfaces:**
- Consumes: `region_magnitudes` and array-accepting `complex_to_hsv`/`complex_to_rgb` from Task 1.
- Produces: `plot_wavefunction_2d`, `WavefunctionArtist`, and `animate_wavefunction` all accept `mag: float | NDArray[np.float64]`. No new public names.

**The one non-mechanical part:** `contour_levels` (`artist.py:34`) computes magnitude contour levels as `k * mag / 5.0`. Contour levels must be scalars, so an array `mag` cannot flow into it. Resolve it by deriving a scalar reference once in `WavefunctionArtist.__init__` — `float(np.max(np.asarray(mag)))` — and passing that to `contour_levels`, leaving the array for the image itself. Document that choice in a comment: contours describe the field as a whole, so they key off the brightest region's scale.

- [ ] **Step 1: Write the failing test**

Append to `libs/qscat/tests/test_viz_coloring.py`:

```python
def test_artist_accepts_an_array_mag() -> None:
    # The call chain must carry a per-point scale all the way to the image,
    # while contour levels (which must be scalars) key off the largest region.
    plt = pytest.importorskip("matplotlib.pyplot")
    from qscat.viz import WavefunctionArtist, region_magnitudes

    z = np.ones((8, 10), dtype=np.complex128)
    z[:, 5:] = 1e-6
    mag = region_magnitudes(np.abs(z), axis=1, boundaries=[5])

    fig, ax = plt.subplots()
    try:
        artist = WavefunctionArtist(
            _IdentityProjector(z.shape), ax=ax, mag=mag, contours=3
        )
        changed = artist.update(z.ravel())
        assert changed, "update() returned no artists"
        rgb = artist._image.get_array()
        # Weak region is visible, not black -- the whole point.
        assert float(np.max(rgb[0, 9])) > 0.5
    finally:
        plt.close(fig)
```

`_IdentityProjector` is a minimal stand-in for `EquidistantProjector`. Define it once near the top of the test module:

```python
class _IdentityProjector:
    """Minimal EquidistantProjector stand-in: reshapes a flat state to a grid."""

    def __init__(self, shape: tuple[int, int]) -> None:
        self.shape = shape
        self.x = np.arange(shape[1], dtype=np.float64)
        self.y = np.arange(shape[0], dtype=np.float64)

    def __call__(self, state: np.ndarray) -> np.ndarray:
        return np.asarray(state).reshape(self.shape)
```

**Before writing this test, read `libs/qscat/qscat/viz/artist.py` and match `_IdentityProjector` to the attributes and call form `WavefunctionArtist` actually uses.** If the real projector's interface differs (attribute names, a `project()` method rather than `__call__`), adapt the stand-in to the real one rather than changing `WavefunctionArtist` to suit the test. If constructing `WavefunctionArtist` needs additional required arguments, supply them; do not add defaults to production code to make the test shorter.

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run --no-sync pytest libs/qscat/tests/test_viz_coloring.py::test_artist_accepts_an_array_mag -v`

Expected: FAIL — either a type/broadcast error from `contour_levels` receiving an array, or a numpy "truth value of an array is ambiguous" error.

- [ ] **Step 3: Widen the three signatures**

In `artist.py`: change `WavefunctionArtist.__init__`'s `mag: float` to `mag: float | npt.NDArray[np.float64]`. Keep `self.mag = mag` (the image path wants the array) and add alongside it:

```python
        # Contour levels must be scalars, so they key off the brightest
        # region's scale; the per-point array still drives the image.
        self._contour_mag = float(np.max(np.asarray(mag, dtype=np.float64)))
```

Replace the `self.mag` argument at the `contour_levels(...)` call site with `self._contour_mag`. Leave `contour_levels`' own signature as `mag: float`.

In `plot.py` and `animate.py`: change `mag: float` to `mag: float | npt.NDArray[np.float64]` in the signatures, and update each docstring's `mag` entry to say a per-point array from `region_magnitudes` is accepted.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run --no-sync pytest libs/qscat/tests/test_viz_coloring.py -v`
Expected: all PASS.

Run: `uv run --no-sync pytest libs/qscat/tests -k viz -v`
Expected: no failures — every existing scalar-`mag` caller still works.

- [ ] **Step 5: Type-check**

Run: `uv run mypy libs/qscat`
Expected: no new errors. The library is type-clean today, so any error here is yours.

- [ ] **Step 6: Commit**

```bash
git add libs/qscat/qscat/viz/artist.py libs/qscat/qscat/viz/plot.py \
        libs/qscat/qscat/viz/animate.py libs/qscat/tests/test_viz_coloring.py
git commit -m "feat(viz): carry a per-point mag through artist, plot and animate

Completes the region-scaling path: the array reaches the domain-coloured
image, while contour levels -- which must be scalars -- key off the largest
region's scale, since contours describe the field as a whole.

Scalar mag callers are untouched.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

### Task 3: A `reference:` block in the qscat-run config

**Files:**
- Create: `apps/qscat-run/qscat_run/reference.py`, `apps/qscat-run/tests/test_reference.py`, `apps/qscat-run/examples/n2-ve-vs-houfek.yaml`
- Modify: `apps/qscat-run/qscat_run/config.py`, `apps/qscat-run/qscat_run/artifacts.py`, `apps/qscat-run/README.md`

**Interfaces:**
- Consumes: nothing from Tasks 1-2.
- Produces: `ReferenceSpec(path: str, format: str, label: str | None, channels: tuple[int, ...] | None)`; `ExperimentConfig.reference: tuple[ReferenceSpec, ...] = ()`; and `qscat_run.reference.load_reference(spec: ReferenceSpec, base_dir: Path) -> dict[str, tuple[NDArray, NDArray]]` mapping a series key to `(energies, sigma)`.

**Read the spec correction above before starting.** The reference keeps its own energy axis. It is NOT a column in `cross_section.csv`.

- [ ] **Step 1: Write the failing tests**

Create `apps/qscat-run/tests/test_reference.py`:

```python
"""The `reference:` config block: overlay a published dataset on a computed
cross section.

The loader reads a data file BY PATH. It must never import `validation` --
`test_no_validation_import.py` enforces that layering, and the whole point of
naming the file in the config is that qscat_run stays independent of it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from qscat_run.config import ConfigError, load_config, validate_config
from qscat_run.reference import load_reference

REPO_ROOT = Path(__file__).resolve().parents[3]
HOUFEK = REPO_ROOT / "validation" / "n2" / "data" / "CSVE.V00.J00"


def _write(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "cfg.yaml"
    p.write_text(body)
    return p


BASE = """
molecule: N2
methods: [ti]
observables: [{kind: ve, channels: 2}]
output_dir: runs/x
"""


def test_config_without_reference_has_an_empty_tuple(tmp_path: Path) -> None:
    cfg = load_config(_write(tmp_path, BASE))
    assert cfg.reference == ()


def test_reference_block_is_parsed(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            BASE
            + f"""
reference:
  - path: {HOUFEK}
    format: houfek
    label: "Houfek (2006)"
    channels: [0, 1]
""",
        )
    )
    assert len(cfg.reference) == 1
    ref = cfg.reference[0]
    assert ref.format == "houfek"
    assert ref.label == "Houfek (2006)"
    assert ref.channels == (0, 1)


def test_unknown_format_is_an_actionable_config_error(tmp_path: Path) -> None:
    cfg = load_config(
        _write(tmp_path, BASE + f"\nreference:\n  - path: {HOUFEK}\n    format: nope\n")
    )
    with pytest.raises(ConfigError, match="nope"):
        validate_config(cfg)


def test_missing_file_fails_at_validate_time_not_at_plot_time(tmp_path: Path) -> None:
    cfg = load_config(
        _write(tmp_path, BASE + "\nreference:\n  - path: no/such/file.dat\n    format: houfek\n")
    )
    with pytest.raises(ConfigError, match="no/such/file.dat"):
        validate_config(cfg)


def test_relative_path_resolves_against_the_config_file(tmp_path: Path) -> None:
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "r.dat").write_text("0.1 1.0 2.0\n0.2 3.0 4.0\n")
    cfg = load_config(
        _write(tmp_path, BASE + "\nreference:\n  - path: data/r.dat\n    format: houfek\n")
    )
    validate_config(cfg)  # must not raise
    series = load_reference(cfg.reference[0], tmp_path)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1"}
    energy, sigma = series["ref:ve:ch0"]
    assert energy.tolist() == [0.1, 0.2]
    assert sigma.tolist() == [1.0, 3.0]


def test_houfek_loader_reads_the_committed_dataset() -> None:
    from qscat_run.reference import ReferenceSpec

    spec = ReferenceSpec(path=str(HOUFEK), format="houfek", label=None, channels=(0, 1, 2))
    series = load_reference(spec, REPO_ROOT)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1", "ref:ve:ch2"}
    energy, sigma = series["ref:ve:ch0"]
    assert energy.shape == (400,)
    assert sigma.shape == (400,)
    assert np.all(np.diff(energy) > 0.0)
    assert np.all(sigma >= 0.0)


def test_channels_omitted_loads_every_column(tmp_path: Path) -> None:
    (tmp_path / "r.dat").write_text("0.1 1.0 2.0 3.0\n0.2 4.0 5.0 6.0\n")
    from qscat_run.reference import ReferenceSpec

    spec = ReferenceSpec(path="r.dat", format="houfek", label=None, channels=None)
    series = load_reference(spec, tmp_path)
    assert set(series) == {"ref:ve:ch0", "ref:ve:ch1", "ref:ve:ch2"}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_reference.py -v`
Expected: FAIL — `ModuleNotFoundError: qscat_run.reference`.

- [ ] **Step 3: Implement the loader**

Create `apps/qscat-run/qscat_run/reference.py`:

```python
"""Published reference datasets overlaid on a computed cross section.

Reads a data file BY PATH, named in the config. This module deliberately does
NOT import `validation` -- `qscat_run` must stay independent of it (see
`tests/test_no_validation_import.py`), and naming the file in the config is
what keeps that true while still letting a run cite committed data.

A reference keeps its OWN energy axis. It is never interpolated onto the run's
energies: doing so would fabricate values and present them as someone else's
measurement, in the one figure whose purpose is honest external comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

__all__ = ["ReferenceSpec", "REFERENCE_FORMATS", "load_reference", "resolve_path"]

Series = dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]]


@dataclass(frozen=True)
class ReferenceSpec:
    """One `reference:` entry: a data file, its format, and how to label it."""

    path: str
    format: str
    label: str | None = None
    channels: tuple[int, ...] | None = None


def resolve_path(spec: ReferenceSpec, base_dir: Path) -> Path:
    """Absolute path for `spec.path`, resolved against `base_dir` when relative."""
    p = Path(spec.path)
    return p if p.is_absolute() else (Path(base_dir) / p)


def _load_columns(path: Path, channels: tuple[int, ...] | None) -> Series:
    """Whitespace-delimited table: column 0 is energy (Hartree), columns 1..N
    are sigma (bohr^2) for successive final channels."""
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] < 2:
        raise ValueError(
            f"{path}: expected a 2-D table with an energy column and at least one "
            f"cross-section column, got shape {raw.shape}"
        )
    energy = raw[:, 0].astype(np.float64)
    n_channels = raw.shape[1] - 1
    wanted = tuple(range(n_channels)) if channels is None else channels
    bad = [c for c in wanted if c < 0 or c >= n_channels]
    if bad:
        raise ValueError(
            f"{path}: requested channel(s) {bad} but the file has {n_channels} "
            f"(valid 0..{n_channels - 1})"
        )
    return {f"ref:ve:ch{c}": (energy, raw[:, c + 1].astype(np.float64)) for c in wanted}


REFERENCE_FORMATS = {"houfek": _load_columns}


def load_reference(spec: ReferenceSpec, base_dir: Path) -> Series:
    """Load one reference dataset as `{series_key: (energies, sigma)}`."""
    loader = REFERENCE_FORMATS.get(spec.format)
    if loader is None:
        raise ValueError(
            f"unknown reference format {spec.format!r}; "
            f"choose one of {sorted(REFERENCE_FORMATS)}"
        )
    return loader(resolve_path(spec, base_dir), spec.channels)
```

- [ ] **Step 4: Wire it into the config**

In `apps/qscat-run/qscat_run/config.py`:

- Import `ReferenceSpec` from `qscat_run.reference` (a leaf module — no cycle).
- Add a loader beside the other `_load_*` helpers:

```python
def _load_reference(raw: list[Any] | None) -> tuple[ReferenceSpec, ...]:
    if not raw:
        return ()
    out: list[ReferenceSpec] = []
    for i, item in enumerate(raw):
        if not isinstance(item, dict) or "path" not in item:
            raise ConfigError(f"reference[{i}] must be a mapping with a 'path' key")
        chans = item.get("channels")
        out.append(
            ReferenceSpec(
                path=str(item["path"]),
                format=str(item.get("format", "houfek")),
                label=None if item.get("label") is None else str(item["label"]),
                channels=None if chans is None else tuple(int(c) for c in chans),
            )
        )
    return tuple(out)
```

- Add `reference: tuple[ReferenceSpec, ...] = ()` to `ExperimentConfig`, after `artifacts`.
- In `load_config`, pass `reference=_load_reference(raw.get("reference"))`.
- **`load_config` must record the config file's directory** so relative reference paths resolve against it. Add a field `config_dir: str | None = None` to `ExperimentConfig` and set it to `str(Path(path).resolve().parent)` in `load_config`. Keep it last so existing positional construction in tests is unaffected.
- In `validate_config`, after the existing checks, add:

```python
    from qscat_run import reference as _reference

    base = Path(cfg.config_dir) if cfg.config_dir else Path.cwd()
    for i, ref in enumerate(cfg.reference):
        if ref.format not in _reference.REFERENCE_FORMATS:
            raise ConfigError(
                f"reference[{i}]: unknown format {ref.format!r}; "
                f"choose one of {sorted(_reference.REFERENCE_FORMATS)}"
            )
        resolved = _reference.resolve_path(ref, base)
        if not resolved.is_file():
            raise ConfigError(
                f"reference[{i}]: no such file {ref.path!r} (looked at {resolved})"
            )
```

Also extend `validate_config`'s docstring check-order sentence to mention the reference checks.

- [ ] **Step 5: Run the config tests**

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_reference.py apps/qscat-run/tests/test_config.py -v`
Expected: all PASS.

Run: `uv run --no-sync pytest apps/qscat-run/tests/test_no_validation_import.py -v`
Expected: PASS — the layering rule still holds.

- [ ] **Step 6: Write the reference into the artifacts**

In `apps/qscat-run/qscat_run/artifacts.py`:

- Add a writer beside the existing ones:

```python
def _write_reference(
    out_dir: Path, series: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]]
) -> None:
    """Reference datasets keep their OWN energy axis, so they get their own
    files rather than columns in `cross_section.csv` (whose rows are the run's
    energies). Interpolating published data onto our grid would fabricate
    values and present them as the reference."""
    if not series:
        return
    with (out_dir / "reference.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "energy", "sigma"])
        for key, (energy, sigma) in series.items():
            for e, s in zip(energy, sigma, strict=True):
                w.writerow([key, float(e), float(s)])
    flat: dict[str, npt.NDArray[np.float64]] = {}
    for key, (energy, sigma) in series.items():
        flat[f"{key}:energy"] = energy
        flat[f"{key}:sigma"] = sigma
    np.savez(out_dir / "reference.npz", **flat)  # type: ignore[arg-type]
```

- Give `_write_cross_section_png` an extra keyword argument
  `reference: dict[str, tuple[NDArray, NDArray]] | None = None`, and after the
  existing solver-series loop add:

```python
    for key, (r_energy, r_sigma) in (reference or {}).items():
        masked = np.where(r_sigma > 0.0, r_sigma, np.nan)
        ax.plot(r_energy, masked, "--", linewidth=1.0, alpha=0.8, label=key)
```

  Dashed and thinner so the reference is visually distinct from computed curves.

- In `write_artifacts`, load each `cfg.reference` entry (resolving against
  `cfg.config_dir`, falling back to the current directory), merge the series
  dicts, pass them to `_write_cross_section_png`, and call `_write_reference`.

- [ ] **Step 7: Add an example config**

Create `apps/qscat-run/examples/n2-ve-vs-houfek.yaml`:

```yaml
# N2 vibrational excitation vs Karel Houfek's independent published data.
#
# The reference file is COMMITTED (validation/n2/data/CSVE.V00.J00), so this
# comparison is reproducible from a bare clone. qscat-run reads it by path and
# does not import `validation` -- the layering rule holds.
#
# This is a PRODUCTION deck: `grid: {preset: emoscat}` is the real,
# convergence-tested grid, not the tiny fast grids the other examples use.
# Run it under Docker: docker/run.sh apps/qscat-run/examples/n2-ve-vs-houfek.yaml runs/n2-ve
molecule: N2
methods: [ti]
observables:
  - {kind: ve, channels: [0, 1, 2]}
energies: {min: 0.04, max: 0.18, step: 0.001}
grid: {preset: emoscat}
v_init: 0
reference:
  - path: ../../../validation/n2/data/CSVE.V00.J00
    format: houfek
    label: "Houfek (2006) CSVE.V00.J00"
    channels: [0, 1, 2]
artifacts: {cross_section: true}
backend: auto
output_dir: runs/n2-ve-vs-houfek
```

Verify the relative path resolves from the example's own directory:
`uv run qscat-run validate apps/qscat-run/examples/n2-ve-vs-houfek.yaml`
Expected: no output (valid). If the path is wrong, `validate_config` says so with the path it looked at — fix the YAML, not the resolver.

- [ ] **Step 8: Run the qscat-run suite**

Run: `uv run --no-sync pytest apps/qscat-run -q`
Expected: no failures. `test_examples.py` globs the examples directory, so the new config is schema-gated automatically.

- [ ] **Step 9: Document it**

In `apps/qscat-run/README.md`, add a row to the "Observables → config knob → artifact" table:

| **published reference overlay** | `reference: [{path: ..., format: houfek, channels: [...]}]` | `reference.{csv,npz}` + dashed overlay on `cross_section.png` (keys `ref:...`) |

And add a short paragraph under it noting that a reference keeps its own energy axis — it is not interpolated onto the run's energies, and so is written separately rather than as extra `cross_section.csv` columns.

- [ ] **Step 10: Commit**

```bash
git add apps/qscat-run/qscat_run/reference.py apps/qscat-run/qscat_run/config.py \
        apps/qscat-run/qscat_run/artifacts.py apps/qscat-run/tests/test_reference.py \
        apps/qscat-run/examples/n2-ve-vs-houfek.yaml apps/qscat-run/README.md
git commit -m "feat(qscat-run): overlay a published reference dataset from config

The flagship gallery figure compares N2 vibrational excitation against
Houfek's independent published data. That overlay previously existed only
inside validation/n2/ti_curve.py, so the most important figure was the one
NOT produced from a config.

The loader reads the data file by path, so qscat_run still does not import
`validation` and the layering rule holds.

A reference keeps its own energy axis and is written to reference.{csv,npz}
rather than as columns in cross_section.csv, whose rows are the run's
energies. Interpolating published data onto our grid would fabricate values
and present them as someone else's measurement, in exactly the figure whose
purpose is honest external comparison.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

## Definition of done

- `region_magnitudes` is public, documented under `docs/api/viz.md`, and the two-region test proves a 1e-6 region renders visibly instead of black.
- Every existing scalar-`mag` caller still works; `uv run mypy libs/qscat` is clean.
- `reference:` parses, validates with actionable errors for an unknown format and a missing file, and resolves relative paths against the config file.
- `apps/qscat-run/tests/test_no_validation_import.py` still passes.
- `uv run qscat-run validate apps/qscat-run/examples/n2-ve-vs-houfek.yaml` is clean.
- `uv run --no-sync pytest tests/test_api_docs_coverage.py` passes (a new public name without docs fails it).
- `uv run sphinx-build -b html -W --keep-going docs docs/_build/html` succeeds.
- `uv run --no-sync pytest -n 8 -m "not slow" -q` is green.

## Explicitly not in this plan

The gallery itself: the eight pages, the production runs on `sadaharu`, the provenance layout, the figure dedup, and the README rework — including the stale `validation/diatomic/da_curves.py` reference at `README.md:104`, which stays broken. Those go in the gallery plan, written once this merges and a first real run exists to quote numbers from.

Also out: the print-mode brightness inversion still deferred in `coloring.py`, and tuning the split radius or percentile for any specific molecule — that is a per-figure choice the gallery plan makes with real data in front of it.
