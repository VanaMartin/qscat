"""Extract the O2 curves of Alt & Houfek, Phys. Rev. A 103, 032829 (2021), Fig. 2.

The published figure (p. 032829-3) is embedded as vector graphics: every curve
is a thin FILLED polyline outline of the original stroked line, dashes and
dots are sub-paths of one filled path, and the tick labels are glyph outlines. So the
curves are recovered EXACTLY rather than digitised by eye: the axes are
calibrated from the tick marks (ten on `R`, 1.5..6.0 bohr; eight on `E`,
-6..+1 eV), and each curve's centreline is the median of its outline points
per 0.6-pt bin in `x` (0.016 bohr). The outline half-width is ~0.3 pt, so the
vertical precision is about 0.02 eV; that is the uncertainty floor the O2
target carries.

Curves taken (colour as filled in the PDF):
    v0           the 3Sigma_g^- neutral curve (black, full)
    v_ion_bound  the 2Pi_g anion curve where it is bound, R >= 2.289 (blue, full)
    e_res_dashed the real part of the resonance energy, R < 2.289 (green, dashed)
    gamma_x2     the local width Gamma(R) x 2 (black, dash-dotted), 0 beyond 2.41

Run: `python -m validation.factory.extract_fig2` with the paper's PDF at
`reference/literature/alt-houfek-2021-pra103-032829.pdf` (gitignored, see
reference/literature/README.md). Writes `validation/factory/data/o2/*.csv`.
Needs `pymupdf` (`uv run --with pymupdf ...`).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

PDF = Path("reference/literature/alt-houfek-2021-pra103-032829.pdf")
OUT = Path("validation/factory/data/o2")
PAGE = 2  # zero-based; printed page 032829-3
PLOT = (340.0, 535.0, 545.0, 675.0)  # PDF points, the region holding Fig. 2
CURVES = {"v0": 306, "v_ion_bound": 308, "e_res_dashed": 309, "gamma_x2": 307}
BIN_PT = 0.6
PRECISION_EV = 0.02


def _points(d, step: float = 0.2) -> np.ndarray:
    """Points ON the outline, sampled every `step` pt along each segment.

    The outlines are polylines (the EPS converter flattened every stroke to
    line segments; V_0's outline has 96 of them, ~one vertex per 0.09 bohr per
    edge, which is the figure's own resolution). Binning the bare vertices
    under-samples the two edges where the curve bends; sampling along the
    segments lets each bin's median average the edges properly.
    """
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            a, b = np.array([it[1].x, it[1].y]), np.array([it[2].x, it[2].y])
            n = max(2, int(np.ceil(np.linalg.norm(b - a) / step)) + 1)
            for t in np.linspace(0.0, 1.0, n):
                q = a + t * (b - a)
                pts.append((float(q[0]), float(q[1])))
        elif it[0] == "c":
            p0, p1, p2, p3 = (np.array([q.x, q.y]) for q in it[1:5])
            for t in np.linspace(0.0, 1.0, 9):
                b = (
                    (1 - t) ** 3 * p0
                    + 3 * (1 - t) ** 2 * t * p1
                    + 3 * (1 - t) * t**2 * p2
                    + t**3 * p3
                )
                pts.append((float(b[0]), float(b[1])))
        elif it[0] == "re":
            r = it[1]
            pts += [(r.x0, r.y0), (r.x1, r.y1)]
    return np.array(pts)


def _ticks(draw, plot) -> tuple[np.ndarray, np.ndarray]:
    import pymupdf

    region = pymupdf.Rect(*plot)
    xt, yt = [], []
    for d in draw:
        if not d["rect"].intersects(region):
            continue
        col = d.get("fill") or d.get("color")
        if col is None or tuple(round(c, 2) for c in col) != (0.0, 0.0, 0.0):
            continue
        r = d["rect"]
        if r.width < 1.5 and 2 < r.height < 8 and r.y1 > 650:
            xt.append((r.x0 + r.x1) / 2)
        if r.height < 1.5 and 2 < r.width < 8 and r.x0 < 372:
            yt.append((r.y0 + r.y1) / 2)

    def merge(v: list[float]) -> np.ndarray:  # each tick is two path pieces 0.1 pt apart
        a = np.array(sorted(v))
        return np.array([g.mean() for g in np.split(a, np.flatnonzero(np.diff(a) > 0.5) + 1)])

    return merge(xt), merge(yt)


def _bin(P: np.ndarray, axis: int, binw: float) -> np.ndarray:
    """Median of the outline points per bin along `axis`; returns (x, y) centreline points."""
    edges = np.arange(P[:, axis].min(), P[:, axis].max() + binw, binw)
    k = np.digitize(P[:, axis], edges)
    rows = []
    for kk in np.unique(k):
        m = k == kk
        if m.sum() < 2:
            continue
        rows.append(
            (np.median(P[m, 0]), np.median(P[m, 1]))
            if axis == 1
            else (P[m, 0].mean(), np.median(P[m, 1]))
        )
    return np.array(rows)


def _centreline(P: np.ndarray, ax: float, bx: float, ay: float, by: float) -> np.ndarray:
    """Centreline in (R, E).

    x-binned where the curve is flat; where it is steep (|dy/dx| > 1 in PDF
    points) the outline's two edges are far apart in y within one x-bin, so
    those runs are re-binned in y instead -- but only from the outline points
    inside that run's own x-window, otherwise a y-bin at a given energy would
    also collect the OTHER wall of the well and report their midpoint.
    """
    cx = _bin(P, 0, BIN_PT)
    slope = np.gradient(cx[:, 1], cx[:, 0]) if cx.shape[0] > 2 else np.zeros(cx.shape[0])
    steep = np.abs(slope) > 1.0
    pieces = [cx[~steep]]
    # contiguous steep runs
    idx = np.flatnonzero(steep)
    if idx.size:
        runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
        for run in runs:
            x_lo, x_hi = cx[run[0], 0] - BIN_PT, cx[run[-1], 0] + BIN_PT
            sel = P[(P[:, 0] >= x_lo) & (P[:, 0] <= x_hi)]
            if sel.shape[0] >= 4:
                pieces.append(_bin(sel, 1, BIN_PT))
    pts = np.vstack(pieces)
    pts = pts[np.argsort(pts[:, 0])]
    # merge points closer than 0.02 pt in x (x-binned and y-binned pieces overlap
    # at a run's ends); a curve table must be strictly increasing in R
    groups = np.split(np.arange(pts.shape[0]), np.flatnonzero(np.diff(pts[:, 0]) > 0.02) + 1)
    pts = np.array([pts[g].mean(axis=0) for g in groups])
    return np.column_stack([ax * pts[:, 0] + bx, ay * pts[:, 1] + by])


def extract(pdf: Path = PDF, out: Path = OUT) -> dict[str, np.ndarray]:
    import pymupdf

    page = pymupdf.open(pdf)[PAGE]
    draw = page.get_drawings()
    xt, yt = _ticks(draw, PLOT)
    if xt.size != 10 or yt.size != 8:
        raise RuntimeError(f"expected 10 x ticks and 8 y ticks, found {xt.size} and {yt.size}")
    ax, bx = np.polyfit(xt, 1.5 + 0.5 * np.arange(10), 1)
    ay, by = np.polyfit(yt, (-6.0 + np.arange(8))[::-1], 1)
    out.mkdir(parents=True, exist_ok=True)
    curves: dict[str, np.ndarray] = {}
    for name, idx in CURVES.items():
        P = _points(draw[idx])
        arr = _centreline(P, ax, bx, ay, by)
        curves[name] = arr
        np.savetxt(
            out / f"{name}.csv",
            arr,
            delimiter=",",
            header=(
                "extracted from Alt & Houfek, PRA 103, 032829 (2021) Fig. 2 (p. 032829-3)"
                f" vector paths by validation/factory/extract_fig2.py;"
                f" vertical precision ~{PRECISION_EV} eV\n"
                "R_bohr,E_eV"
            ),
            comments="# ",
        )
    return curves


if __name__ == "__main__":
    c = extract()
    for name, arr in c.items():
        j = int(np.argmin(arr[:, 1]))
        print(
            f"{name}: {arr.shape[0]} pts, R[{arr[0, 0]:.3f}, {arr[-1, 0]:.3f}], "
            f"min E = {arr[j, 1]:.3f} eV at R = {arr[j, 0]:.3f}"
        )
