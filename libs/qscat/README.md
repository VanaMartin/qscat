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
| `qscat.viz` | wavefunction rendering and animation (needs the `plot` extra) |

## Relation to existing work

The method is not new: qscat implements Rescigno & McCurdy's FEM-DVR with
exterior complex scaling (Phys. Rev. A **62**, 032706, 2000) applied to the
exactly-solvable two-dimensional resonant-collision model of Houfek, Rescigno &
McCurdy (Phys. Rev. A **73**, 032721, 2006). What qscat adds is a released,
validated implementation of the scattering observables built on it.

The grid layer overlaps one other published package,
[`quantumgrid`](https://pypi.org/project/quantumgrid/) (McCurdy, Streeter &
Barbalinardo, MIT) — a teaching-oriented FEM-DVR-ECS package for the
one-dimensional Schrödinger equation, with no scattering observables. The
established *ab initio* electron–molecule suites (UKRmol+, ePolyScat, FERM3D,
Quantemol-EC) solve a different problem by a different method and are complements
rather than alternatives. A survey of what is and is not already published as
code — including what appears to have no released counterpart, and the limits of
that claim — is in
[`docs/related-work.md`](https://github.com/VanaMartin/qscat/blob/main/docs/related-work.md).

## Citing qscat

If you use qscat in research, please cite the software (see
[`CITATION.cff`](https://github.com/VanaMartin/qscat/blob/main/CITATION.cff))
together with the papers below.

**The method this code implements:**

- M. Váňa and K. Houfek, *Time-dependent formulation of the two-dimensional
  model of resonant electron collisions with diatomic molecules and
  interpretation of the vibrational excitation cross sections*,
  Phys. Rev. A **95**, 022714 (2017).
  [doi:10.1103/PhysRevA.95.022714](https://doi.org/10.1103/PhysRevA.95.022714)

**The model itself:**

- K. Houfek, T. N. Rescigno and C. W. McCurdy, *Numerically solvable model for
  resonant collisions of electrons with diatomic molecules*,
  Phys. Rev. A **73**, 032721 (2006).
  [doi:10.1103/PhysRevA.73.032721](https://doi.org/10.1103/PhysRevA.73.032721)

**The numerical methods**, if you use the grid or propagator directly:

- T. N. Rescigno and C. W. McCurdy, *Numerical grid methods for
  quantum-mechanical scattering problems*, Phys. Rev. A **62**, 032706 (2000)
  — the FEM-DVR grid with exterior complex scaling.
  [doi:10.1103/PhysRevA.62.032706](https://doi.org/10.1103/PhysRevA.62.032706)
- D. J. Tannor and D. E. Weeks, *Wave packet correlation function formulation
  of scattering theory*, J. Chem. Phys. **98**, 3884 (1993) — the
  correlation-function energy extraction behind every time-dependent cross
  section. [doi:10.1063/1.464016](https://doi.org/10.1063/1.464016)
- W. van Dijk and F. M. Toyama, *Accurate numerical solutions of the
  time-dependent Schrödinger equation*, Phys. Rev. E **75**, 036707 (2007) —
  the order-N diagonal-Padé propagator.
  [doi:10.1103/PhysRevE.75.036707](https://doi.org/10.1103/PhysRevE.75.036707)

## License

BSD-3-Clause © 2026 Martin Vana. See
[`LICENSE`](https://github.com/VanaMartin/qscat/blob/main/LICENSE).
