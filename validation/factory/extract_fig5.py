"""Vector-extract the VE cross sections of Alt & Houfek, PRA 103, 032829 (2021),
Fig. 5 (p. 032829-6): the paper's own NONLOCAL (NRM) and LOCAL (LCP) curves
for 0 -> v', v' = 0..5, spin-orbit resolved -- theory curves, for the
theory-vs-theory overlay of the factory model's spin-orbit-resolved exact 2-D
cross section. Nothing experimental is on this page.

`uv run --with pymupdf python -m validation.factory.extract_fig5`

Same route as `extract_fig2.py`: the figure is embedded as vector paths (the
NRM curve a filled polyline outline in blue, the LCP a dashed one in green,
2-11k vertices per panel). Each panel's axes are calibrated from its tick
RECTANGLES (the labels are glyph outlines, not text; their values are read
off the rendered page once and fixed in `PANELS`). The centreline is the
UPPER ENVELOPE of the outline per x-bin, shifted down by half the stroke
width: a comb of meV-wide peaks is where the Fig. 2 median-of-edges
centreline fails (an x-bin on a peak holds both walls), while the envelope
keeps every peak's height to half a stroke width (0.25 pt, i.e. 0.3-0.7 % of
each panel's range). Energy resolution is one bin, `BIN_PT` = 0.15 pt:
0.9-1.4 meV depending on the panel.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from validation.factory.extract_fig2 import PDF, _points

__all__ = ["PANELS", "extract", "main"]

OUT = Path(__file__).parent / "data" / "o2"
PAGE = 5  # zero-based; printed page 032829-6
BIN_PT = 0.15
HALF_STROKE_PT = 0.25
BLUE = (0.0, 0.45, 0.74)  # NRM
GREEN = (0.47, 0.67, 0.19)  # LCP
# The legend key samples, relative to each panel's frame (measured on the
# page): x from `x1 - 36.5` to `x1 - 9` pt, the NRM line 7.25 pt and the LCP
# dashes 16.2 pt below the top edge, each stroke ~0.5 pt tall.
KEY_X_PT = (36.5, 9.0)
KEY_ROWS_PT = (7.25, 16.2)
KEY_ROW_HALF_PT = 1.0

# Per panel: v', the frame (x0, x1, y0, y1) in PDF points, and the tick VALUES
# (ascending in energy / cross section) matching the tick rectangles found on
# the frame's bottom / left edges. Read off the rendered page (2026-08-25).
PANELS: tuple[dict[str, object], ...] = (
    {
        "v": 0,
        "frame": (133.26, 290.54, 88.27, 209.7),
        "x": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        "y": [100.0 * k for k in range(11)],
    },
    {
        "v": 1,
        "frame": (350.35, 507.72, 88.27, 209.7),
        "x": [0.4, 0.6, 0.8, 1.0, 1.2],
        "y": [5.0 * k for k in range(8)],
    },
    {
        "v": 2,
        "frame": (133.26, 290.63, 264.97, 386.4),
        "x": [0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
        "y": [0.5 * k for k in range(11)],
    },
    {
        "v": 3,
        "frame": (350.35, 507.72, 264.97, 386.4),
        "x": [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
        "y": [0.1 * k for k in range(11)],
    },
    {
        "v": 4,
        "frame": (133.26, 290.63, 441.62, 563.09),
        "x": [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
        "y": [0.02 * k for k in range(9)],
    },
    {
        "v": 5,
        "frame": (355.69, 507.72, 441.62, 563.09),
        "x": [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2],
        "y": [0.005 * k for k in range(5)],
    },
)


def _tick_rects(draw) -> np.ndarray:
    rects = []
    for d in draw:
        if d.get("fill") is None or tuple(round(c, 2) for c in d["fill"]) != (0.0, 0.0, 0.0):
            continue
        for it in d["items"]:
            if it[0] == "re":
                r = it[1]
                rects.append((r.x0, r.y0, r.x1, r.y1))
    return np.array(rects)


def _merge(v: np.ndarray) -> np.ndarray:
    a = np.array(sorted(np.round(v, 2)))
    if a.size == 0:
        return a
    return np.array([g.mean() for g in np.split(a, np.flatnonzero(np.diff(a) > 0.6) + 1)])


def _calibrate(rects: np.ndarray, panel: dict[str, object]) -> tuple[float, float, float, float]:
    """Linear maps `E = ax * x + bx`, `sigma = ay * y + by` from the panel's tick
    rectangles on its bottom (x) and left (y) frame edges."""
    x0, x1, y0, y1 = panel["frame"]  # type: ignore[misc]
    w = rects[:, 2] - rects[:, 0]
    h = rects[:, 3] - rects[:, 1]
    mx = (w < 0.6) & (h > 1.5) & (h < 8) & (np.abs(rects[:, 3] - y1) < 0.6)
    mx &= (rects[:, 0] > x0 - 1) & (rects[:, 2] < x1 + 1)
    my = (h < 0.6) & (w > 1.5) & (w < 8) & (np.abs(rects[:, 0] - x0) < 0.6)
    my &= (rects[:, 1] > y0 - 1) & (rects[:, 3] < y1 + 1)
    xt = _merge((rects[mx, 0] + rects[mx, 2]) / 2)
    yt = _merge((rects[my, 1] + rects[my, 3]) / 2)
    xv = np.asarray(panel["x"], dtype=float)
    yv = np.asarray(panel["y"], dtype=float)[::-1]  # PDF y grows downward
    if xt.size != xv.size or yt.size != yv.size:
        raise RuntimeError(
            f"panel v'={panel['v']}: {xt.size} x-ticks / {yt.size} y-ticks found, "
            f"{xv.size} / {yv.size} values given"
        )
    ax, bx = np.polyfit(xt, xv, 1)
    ay, by = np.polyfit(yt, yv, 1)
    return float(ax), float(bx), float(ay), float(by)


def _envelope(P: np.ndarray, ax: float, bx: float, ay: float, by: float) -> np.ndarray:
    """Upper envelope of the outline per x-bin, minus half the stroke, in (E, sigma)."""
    edges = np.arange(P[:, 0].min(), P[:, 0].max() + BIN_PT, BIN_PT)
    k = np.digitize(P[:, 0], edges)
    rows = []
    for kk in np.unique(k):
        m = k == kk
        rows.append((P[m, 0].mean(), P[m, 1].min() + HALF_STROKE_PT))
    pts = np.array(rows)
    return np.column_stack([ax * pts[:, 0] + bx, ay * pts[:, 1] + by])


def _curve_paths(draw, colour: tuple[float, float, float], frame) -> tuple[list, list]:
    """`(curves, swatches)` of one colour inside `frame`: the curve paths (hundreds
    of vertices) and the legend key samples (a few items). gnuplot strokes the
    key sample with the SAME path as the curve, so the NRM key line is part of
    the NRM curve's own drawing: the swatch boxes are used to mask it out."""
    x0, x1, y0, y1 = frame
    curves, swatches = [], []
    for d in draw:
        f = d.get("fill")
        if f is None or tuple(round(c, 2) for c in f) != colour:
            continue
        r = d["rect"]
        inside = r.x0 >= x0 - 2 and r.x1 <= x1 + 2 and r.y0 >= y0 - 8 and r.y1 <= y1 + 2
        if not inside:
            continue
        (curves if len(d["items"]) > 100 else swatches).append(d)
    return curves, swatches


def _mask_legend(P: np.ndarray, frame) -> np.ndarray:
    """Drop the outline points of the legend key samples. gnuplot strokes each
    key sample with the SAME path as its curve (on five of the six panels no
    separate swatch path exists), at a fixed place relative to the frame:
    `KEY_X_PT` from the right edge, the NRM line `KEY_ROWS_PT[0]` and the LCP
    dashes `KEY_ROWS_PT[1]` below the top edge (measured on the page). Only
    two thin bands are removed, so a real peak crossing them keeps its top
    to within a band's half-height (0.0002 of a panel's range)."""
    _x0, x1, y0, _y1 = frame
    in_x = (P[:, 0] >= x1 - KEY_X_PT[0]) & (P[:, 0] <= x1 - KEY_X_PT[1])
    in_rows = np.zeros(P.shape[0], dtype=bool)
    for row in KEY_ROWS_PT:
        in_rows |= np.abs(P[:, 1] - (y0 + row)) < KEY_ROW_HALF_PT
    return P[~(in_x & in_rows)]


FIGURES = (("fig5", PAGE, "Fig. 5 (p. 032829-6)", PANELS, (("nrm", BLUE), ("lcp", GREEN))),)


def extract(pdf: Path = PDF, out: Path = OUT) -> dict[str, np.ndarray]:
    import pymupdf

    doc = pymupdf.open(pdf)
    out.mkdir(parents=True, exist_ok=True)
    result: dict[str, np.ndarray] = {}
    for tag, page_no, where, panels, curves in FIGURES:
        draw = doc[page_no].get_drawings()
        rects = _tick_rects(draw)
        header = (
            f"extracted from Alt & Houfek, PRA 103, 032829 (2021) {where} vector "
            "paths by validation/factory/extract_fig5.py; upper-envelope centreline, "
            f"resolution one {BIN_PT}-pt bin, heights to half a stroke width (0.25 pt)"
        )
        for panel in panels:
            cal = _calibrate(rects, panel)
            v = panel["v"]
            for name, colour in curves:
                paths, _ = _curve_paths(draw, colour, panel["frame"])
                if not paths:
                    raise RuntimeError(f"{tag} panel v'={v}: no {name} curve found")
                P = _mask_legend(np.vstack([_points(d) for d in paths]), panel["frame"])
                curve = _envelope(P, *cal)
                curve = curve[np.argsort(curve[:, 0])]
                key = f"{tag}_ve_0{v}_{name}"
                result[key] = curve
                np.savetxt(
                    out / f"{key}.csv",
                    curve,
                    delimiter=",",
                    header=f"{header}\nE_eV,sigma_a0^2 (0 -> {v}, {name.upper()})",
                    comments="# ",
                )
    return result


def main() -> None:
    curves = extract()
    for k, c in curves.items():
        span = f"[{c[0, 0]:.3f}, {c[-1, 0]:.3f}] eV"
        print(f"{k}: {c.shape[0]} points, E in {span}, max {c[:, 1].max():.4g}")


if __name__ == "__main__":
    main()
