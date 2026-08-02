# qscat

**Validated quantum-scattering numerics for electron–diatomic collisions.**

`qscat` is a CPU-first, laptop-runnable Python library for computing
electron–molecule scattering observables with the **FEM-DVR** (finite-element
discrete-variable representation) method and **exterior complex scaling (ECS)**.
It provides the numerical building blocks — sparse tensor Hamiltonians, ECS
grids, direct sparse solvers, time propagators, special functions — and the
model-independent scattering engine on top of them: vibrational-excitation (VE),
dissociative-attachment (DA), and dissociative-recombination (DR) cross sections,
by both time-independent (driven Lippmann–Schwinger) and time-dependent
(wavepacket-propagation) routes. Everything is in **atomic units**.

Every method is checked against an analytic benchmark, a conservation law, a
convergence study, or an independent reference implementation before it ships.

## Installation

```bash
pip install qscat                 # core (numpy, scipy, mpmath only)
pip install "qscat[plot]"         # + matplotlib, for the figure helpers
```

`qscat` imports with only numpy/scipy/mpmath; matplotlib is an optional extra
used lazily by `qscat.core.plot_cross_sections`.

**Optional MUMPS backend.** `qscat.linalg.SparseLU` can use a complex-symmetric
MUMPS solver (much faster and lighter than SuperLU on the large ECS matrices).
It builds against a *system* MUMPS, so `pip install "qscat[mumps]"` only works
where one is present; the supported route is a prebuilt MUMPS (e.g. conda-forge).
Without it, the SciPy SuperLU backend is used automatically.

## Quickstart

An N₂ vibrational-excitation cross section, σ₀→v′(E), on a small grid:

```python
import numpy as np
from qscat.core import ScatteringProblem
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.model import N2

# Tensor-product FEM-DVR-ECS grid: electronic (r) × nuclear (R)
grid = TensorGrid([
    electronic_grid(r_max=16.0, order=7, n_complex=5),
    nuclear_grid(r_max=22.0, quadrature=10, n_complex=5),
])

# One problem = grid + model + how many vibrational states (solved once, reused)
prob = ScatteringProblem(grid=grid, model=N2, n_vib=4, v_init=0)

# σ for v=0 → v'=0,1,2 at three collision energies (Hartree)
sigma = prob.ve_cross_section(vprimes=[0, 1, 2], E=np.array([0.10, 0.15, 0.20]))
print(sigma)  # shape (3, 3) bohr², [E, v']
```

`ScatteringProblem` is the recommended entry point; it also exposes
`.da_cross_section(...)`, `.dr_cross_section(...)`, and the time-dependent
`.td_ve_cross_section(...)` / `.td_da_cross_section(...)`. The underlying
functional solvers (`qscat.core.ve_cross_section`, …) remain available as the
low-level layer.

See the `docs/` tree in the repository for the theory notes and more examples
(TI vs TD, DA/DR, the discretisation tuner, and the `qscat-run` config CLI).

## What's inside

| Module | Purpose |
|---|---|
| `qscat.units` | atomic-unit constants |
| `qscat.linalg` | sparse `kron_sum`, `SparseLU` (SuperLU/MUMPS), ECS `c_product` |
| `qscat.special` | Riccati–Bessel / Coulomb radial functions |
| `qscat.dvr` | FEM-DVR-ECS grids, kinetic assembly, N-D tensor Hamiltonians |
| `qscat.ecs` | exterior-complex-scaling map + resonance-pole finder |
| `qscat.evolution` | Crank–Nicolson and diagonal-Padé time propagators |
| `qscat.core` | model-independent VE/DA/DR engine (TI + TD) |
| `qscat.model` | molecule models (N₂, NO, F₂, H₂⁺) + the `ResonanceModel` protocol |
| `qscat.tuning` | automatic FEM-DVR-ECS grid tuner |

## Citing qscat

**If you use qscat directly in research, citation is required.** Please cite the
software (see [`CITATION.cff`](https://github.com/VanaMartin/qscat/blob/main/CITATION.cff))
and the relevant papers below.

*Method / software paper (in preparation — the paper to cite for the code
itself; update this entry with the final reference on publication):*

- M. Vana et al., *"QSCAT: validated FEM-DVR-ECS numerics for electron–molecule
  scattering,"* in preparation (2026). **[reference to be finalized]**

*Prior published work this code builds on / continues (fill in the exact
citations you want required):*

- **[Add citation — prior article #1: title, authors, journal, year, DOI]**
- **[Add citation — prior article #2: title, authors, journal, year, DOI]**

> Maintainers: keep this list and `CITATION.cff` in sync. The forthcoming
> method paper is the canonical reference for the code; the prior articles are
> the physics/method lineage a direct user is expected to cite.

## License

BSD-3-Clause © 2026 Martin Vana. See [`LICENSE`](LICENSE).
