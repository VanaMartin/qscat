# Getting started

## Install

`qscat` is **repo-only**: it is not published to PyPI (and will not be until
the qscat citation article is out). Install it from a clone:

```bash
git clone https://github.com/VanaMartin/qscat
cd qscat
uv sync --all-packages
```

`qscat` imports with only numpy/scipy/mpmath. `plot` (matplotlib, for the
figure helpers) and `mumps` are optional extras, pulled in by name:

```bash
uv sync --all-packages --extra plot
uv sync --all-packages --extra mumps
```

The MUMPS backend builds against a *system* MUMPS, so that extra only works
where one is present — see the [package
README](https://github.com/VanaMartin/qscat/blob/main/libs/qscat/README.md).

## Your first cross section

```python
import numpy as np
from qscat.core import ScatteringProblem
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2

grid = TensorGrid([
    electronic_grid(r_max=16.0, order=7, n_complex=5),
    nuclear_grid(r_max=22.0, quadrature=10, n_complex=5),
])
prob = ScatteringProblem(grid=grid, model=N2, n_vib=4, v_init=0)

sigma = prob.ve_cross_section(vprimes=[0, 1, 2], E=np.array([0.10, 0.15, 0.20]))
print(sigma)  # (3, 3) bohr², [E, v']
```

`ScatteringProblem` is the recommended entry point. It bundles the grid, model,
and vibrational basis once and exposes every observable as a method:

- `prob.ve_cross_section(vprimes, E)` — vibrational excitation
- `prob.da_cross_section(E)` — dissociative attachment
- `prob.dr_cross_section(E)` — dissociative recombination (ions, e.g. H₂⁺)
- `prob.td_ve_cross_section(vprimes, E, ...)` / `prob.td_da_cross_section(E, ...)`
  — the time-dependent (wavepacket-propagation) routes

## Choosing a grid

Grids are per-potential FEM-DVR-ECS tensor products. For a first pass, use the
`qscat.model` molecules (`N2`, `NO`, `F2`, `H2P`, the fitted `O2`) with grids sized like the
example above; for production, the `qscat.tuning` module derives a minimal grid
at a target precision from the potential and energy range.

## The config CLI

For reproducible, dockerized runs, the `qscat-run` CLI drives whole experiments
from a YAML config (TI/TD, multiple observables, artifacts + manifest). See the
top-level `README.md`.
