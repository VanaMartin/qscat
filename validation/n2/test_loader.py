import numpy as np

from validation.n2 import loader


def test_shape_and_grid():
    d = loader.load()
    assert d.energy.shape == (400,)
    assert d.sigma.shape == (400, 31)
    assert d.energy[0] == 5e-4
    assert abs(d.energy[-1] - 0.2) < 1e-12
    # strictly increasing
    assert np.all(np.diff(d.energy) > 0)


def test_nonnegative_and_elastic_column():
    d = loader.load()
    assert np.all(d.sigma >= 0.0)
    # elastic (v=0->0) is column index 0, grows into the resonance region
    assert d.sigma[-1, 0] > d.sigma[0, 0]


def test_threshold_ordering():
    # channel v=0->(j+1) opens at an energy >= where v=0->j opens (higher channels open later)
    d = loader.load()

    def first_open(j):
        nz = np.nonzero(d.sigma[:, j] > 0)[0]
        return d.energy[nz[0]] if nz.size else np.inf

    opens = [first_open(j) for j in range(1, 31)]  # skip elastic
    finite = [o for o in opens if np.isfinite(o)]
    assert finite == sorted(finite)


def test_integrity_checks_all_pass():
    results = loader.integrity_checks()
    assert results, "expected integrity checks"
    assert all(ok for _name, ok, _detail in results), results
