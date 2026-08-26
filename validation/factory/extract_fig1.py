"""Vector-extract the spin-orbit splitting of the O2^- 2Pi_g curve from Alt &
Houfek, PRA 103, 032829 (2021), Fig. 1 (p. 032829-3): the Gaussian FIT the
paper itself used to build the 2Pi_{1/2} / 2Pi_{3/2} curves (symmetric at
+-Delta_SO(R)/2 around 2Pi_g, Sec. III A).

`uv run --with pymupdf python -m validation.factory.extract_fig1`

Same route as `extract_fig2.py`: the fit curve is a filled polyline outline
(12k vertices), the axes are calibrated from the tick rectangles, and the
centreline is the median of the outline's two edges per x-bin -- the curve
is smooth and single-valued, so no steep-run handling is needed. Only the
fit curve is taken; the MOLPRO points it was fitted to are not.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from validation.factory.extract_fig2 import PDF, _bin, _points

__all__ = ["extract", "main"]

OUT = Path(__file__).parent / "data" / "o2" / "so_split.csv"
PAGE = 2  # zero-based; printed page 032829-3
FRAME = (362.31, 532.94, 75.14, 191.06)  # PDF points; R 1..16 bohr, split 12..20 meV
X_TICKS = [2.0 * k for k in range(1, 9)]  # bohr, at the 8 bottom tick rectangles
Y_TICKS = [12.0 + k for k in range(9)]  # meV, at the 9 left tick rectangles
BLUE = (0.0, 0.45, 0.74)
BIN_PT = 0.4


def _merge(v: np.ndarray) -> np.ndarray:
    a = np.array(sorted(np.round(v, 2)))
    if a.size == 0:
        return a
    return np.array([g.mean() for g in np.split(a, np.flatnonzero(np.diff(a) > 0.6) + 1)])


def _calibrate(draw) -> tuple[float, float, float, float]:
    x0, x1, y0, y1 = FRAME
    rects = []
    for d in draw:
        if d.get("fill") != (0.0, 0.0, 0.0) or d["rect"].y1 > y1 + 20:
            continue
        for it in d["items"]:
            if it[0] == "re":
                r = it[1]
                rects.append((r.x0, r.y0, r.x1, r.y1))
    R = np.array(rects)
    w, h = R[:, 2] - R[:, 0], R[:, 3] - R[:, 1]
    mx = (w < 0.6) & (h > 1.5) & (h < 8) & (np.abs(R[:, 3] - y1) < 1.5)
    my = (h < 0.6) & (w > 1.5) & (w < 8) & (np.abs(R[:, 0] - x0) < 1.5)
    xt = _merge((R[mx, 0] + R[mx, 2]) / 2)
    yt = _merge((R[my, 1] + R[my, 3]) / 2)
    if xt.size != len(X_TICKS) or yt.size != len(Y_TICKS):
        raise RuntimeError(f"found {xt.size} x-ticks / {yt.size} y-ticks")
    ax, bx = np.polyfit(xt, X_TICKS, 1)
    ay, by = np.polyfit(yt, Y_TICKS[::-1], 1)  # PDF y grows downward
    return float(ax), float(bx), float(ay), float(by)


def extract(pdf: Path = PDF, out: Path = OUT) -> np.ndarray:
    import pymupdf

    draw = pymupdf.open(pdf)[PAGE].get_drawings()
    ax, bx, ay, by = _calibrate(draw)
    x0, x1, y0, y1 = FRAME
    curves = [
        d
        for d in draw
        if d.get("fill") is not None
        and tuple(round(c, 2) for c in d["fill"]) == BLUE
        and len(d["items"]) > 100
        and d["rect"].y1 <= y1 + 2
    ]
    if len(curves) != 1:
        raise RuntimeError(f"expected one fit curve, found {len(curves)}")
    P = _points(curves[0])
    # the legend key line sits at the top right; drop points above the frame's
    # first 12 pt right of R = 10 bohr, where the curve is flat at 12.2 meV
    key = (P[:, 0] > x0 + 0.55 * (x1 - x0)) & (P[:, 1] < y0 + 14.0)
    P = P[~key]
    c = _bin(P, 0, BIN_PT)
    table = np.column_stack([ax * c[:, 0] + bx, ay * c[:, 1] + by])
    table = table[np.argsort(table[:, 0])]
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        out,
        table,
        delimiter=",",
        header=(
            "extracted from Alt & Houfek, PRA 103, 032829 (2021) Fig. 1 (p. 032829-3) vector "
            "paths by validation/factory/extract_fig1.py -- the paper's Gaussian fit of the "
            "spin-orbit splitting of the O2^- 2Pi_g curve; precision ~0.05 meV\n"
            "R_bohr,split_meV"
        ),
        comments="# ",
    )
    return table


def main() -> None:
    t = extract()
    print(
        f"so_split: {t.shape[0]} points, R in [{t[0, 0]:.2f}, {t[-1, 0]:.2f}] bohr, "
        f"split {t[:, 1].min():.2f}..{t[:, 1].max():.2f} meV, at 2.28 bohr "
        f"{np.interp(2.28, t[:, 0], t[:, 1]):.2f} meV"
    )


if __name__ == "__main__":
    main()
