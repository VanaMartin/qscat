"""The committed Fig. 5 extraction: the paper's NRM/LCP VE curves are what the
extractor reproduces from the PDF (when it is present), and what the tables
claim is sane -- every panel's maximum sits inside its own axis range and the
elastic peak is the ~940 a0^2 the page shows."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data" / "o2"
PDF = Path(__file__).resolve().parents[2] / "reference/literature/alt-houfek-2021-pra103-032829.pdf"

# (v', the NRM curve's tallest peak as printed on the panel, panel y range).
# Before the legend key samples were masked, every "maximum" but the first
# was the key LINE itself at 0.94 of the range -- a lesson kept here.
_PANELS = (
    (0, 939.0, 1000.0),
    (1, 23.7, 35.0),
    (2, 3.75, 5.0),
    (3, 0.562, 1.0),
    (4, 0.0752, 0.16),
    (5, 0.0172, 0.02),
)


@pytest.mark.parametrize("v,peak,ymax", _PANELS)
def test_committed_curves_are_sane(v, peak, ymax):
    for name in ("nrm", "lcp"):
        d = np.loadtxt(DATA / f"fig5_ve_0{v}_{name}.csv", delimiter=",")
        assert np.all(np.diff(d[:, 0]) > 0)  # strictly increasing energy
        assert d[:, 1].max() <= ymax * 1.02 and d[:, 1].min() >= -0.01 * ymax
    nrm = np.loadtxt(DATA / f"fig5_ve_0{v}_nrm.csv", delimiter=",")
    assert abs(nrm[:, 1].max() - peak) < 0.02 * ymax


@pytest.mark.skipif(not PDF.exists(), reason="the source PDF is gitignored")
def test_extractor_reproduces_the_committed_tables(tmp_path):
    pytest.importorskip("pymupdf")
    from validation.factory.extract_fig5 import extract

    curves = extract(out=tmp_path)
    for key, c in curves.items():
        committed = np.loadtxt(DATA / f"{key}.csv", delimiter=",")
        np.testing.assert_allclose(c, committed, rtol=0, atol=1e-9)
