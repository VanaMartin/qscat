# Related work

What qscat shares with code that is already published, and what appears to have
no released counterpart. Surveyed **2026-08-16** across PyPI, GitHub, Zenodo,
and the published literature; the survey's method and its blind spots are stated
at the end, so a reader can judge how much the "no counterpart found" claims are
worth.

## Summary

| Layer | Status |
|---|---|
| FEM-DVR grid, Gauss–Lobatto quadrature, bridge functions, ECS contour, Crank–Nicolson propagation | **Also implemented elsewhere** — `quantumgrid` (Python, MIT), plus several partial implementations in Fortran/Julia |
| Driven Lippmann–Schwinger solve, VE/DA/DR cross sections, LCP approximation, resonance levels, TD extractors, grid tuner | **No released implementation found** in any language |
| Electron–molecule scattering as a research goal | Well served by *ab initio* suites (UKRmol+, ePolyScat, FERM3D) that solve a **different problem** with a different method |

The method qscat implements is not new — it is Rescigno & McCurdy's FEM-DVR-ECS
(2000) applied to Houfek, Rescigno & McCurdy's two-dimensional resonant-collision
model (2006). What appears to be new is that it is *released as a library*.

## The closest relative: quantumGrid

[`quantumgrid`](https://pypi.org/project/quantumgrid/) —
<https://github.com/cwmccurdy/quantumGrid> (fork: `zstreeter/quantumGrid`).

| | |
|---|---|
| Authors | C. W. McCurdy, Z. Streeter, G. Barbalinardo (UC Davis / LBNL) |
| License | MIT |
| Latest release | 2.0.30, April 2021 (PyPI); repository last updated December 2024 |
| Stated purpose | "Exterior Complex Scaling Finite-Element Method Discrete Variable Representation grid for general physics problems" — written for a UC Davis graduate course in time-dependent quantum mechanics and released as open source |

This is the same method lineage as qscat: McCurdy is a co-author of both the
FEM-DVR-ECS method paper and the two-dimensional model qscat solves.

**Shared with qscat** (per its documentation and examples): Gauss–Lobatto
FEM-DVR grid with bridge functions between elements; kinetic-energy matrix
assembly in that basis; ECS applied as a complex final element; a Crank–Nicolson
propagator; propagation on two coupled potential curves; a two-dimensional grid
used for a two-*electron* problem (`Hamiltonian_2D`, its `Two-electron_ECS`
example); built-in Morse and Coulomb potentials.

**Not in it:** any scattering observable. There is no driven-equation or
Lippmann–Schwinger solve, no cross section of any kind, no local-complex-potential
approximation, no resonance-pole search or complex quasi-bound levels, no
propagator beyond Crank–Nicolson, no time-dependent energy-extraction transform,
no sparse or out-of-core linear algebra, no grid-tuning, and no Coulomb (ionic)
channels. Its two-dimensional example couples two electronic coordinates, not an
electronic coordinate to a nuclear one, so it does not express the resonant
electron–diatomic model at all.

**Honest overlap:** `qscat.dvr` and `qscat.ecs` re-implement numerics that
`quantumgrid` already publishes. Both are independent implementations of the same
published method, which makes `quantumgrid` a legitimate external cross-check for
qscat's one-dimensional grid — a use qscat does not currently make of it. Anyone
who needs *only* a FEM-DVR-ECS grid for a one-dimensional problem, and
especially anyone teaching the method, should look at `quantumgrid` first.

## Other implementations of the same numerical ingredients

None of these compute electron–molecule scattering observables; they are grid,
basis, or single-coordinate solvers.

| Code | Language | What it is |
|---|---|---|
| [COLOSS](https://github.com/jinleiphys/COLOSS) — Lei, Liu & Ren, *Comput. Phys. Commun.* (2025), [arXiv:2407.16425](https://arxiv.org/abs/2407.16425) | Fortran | Complex-scaled optical + Coulomb scattering solver for **nuclear** reactions; complex scaling to turn oscillatory boundary conditions into decaying ones |
| [SEECS](https://github.com/banana-bred/SEECS) | Fortran | One-dimensional rovibrational Schrödinger equation with exterior complex scaling |
| [FEDVR.jl](https://github.com/jagot/FEDVR.jl) | Julia | FE-DVR basis package (archived) |
| [LBNL-AMO-MCTDHF](https://github.com/LBNL-AMO-MCTDHF/V1) | Fortran | MCTDHF electron dynamics from the same group; uses FEM-DVR internally, targets a different physical problem |
| [Discvar](https://github.com/KenHino/Discvar), [`dvr_py`](https://github.com/richford/dvr_py), and ~15 similar repositories | Python | DVR basis/Hamiltonian construction; no ECS, no scattering |
| [tRecX](https://arxiv.org/abs/2403.11918) | C++ | General time-dependent Schrödinger solver using infinite-range ECS for strong-field/attosecond problems |

## Electron–molecule scattering codes with a different method and scope

These solve the *ab initio* electron–molecule problem for real molecules —
typically the fixed-nuclei electronic scattering problem, from which nuclear
dynamics is treated separately or not at all. qscat solves a reduced,
exactly-solvable model in which the electronic and nuclear coordinates are
treated on the same footing, without the Born–Oppenheimer or local-complex-potential
approximations. They are complements, not alternatives.

| Code | Method | Availability |
|---|---|---|
| [UKRmol+](https://arxiv.org/abs/1908.03018) | Molecular R-matrix | Open source (Zenodo) |
| [ePolyScat](https://github.com/rrlucchese/ePolyScat) | Single-centre expansion / Schwinger variational | Open source; also hosted on the AMOS Gateway |
| [FERM3D](https://arxiv.org/abs/physics/0607062) | Finite-element R-matrix | Published code |
| Quantemol-EC | Expert system driving UKRmol+ | Commercial |
| [PyOpenCAP](https://github.com/gayverjr/opencap) | Complex absorbing potential on quantum-chemistry data | Open source; extracts resonance position and width, not cross sections |
| [jitr](https://pypi.org/project/jitr/) | Calculable R-matrix on a Lagrange–Legendre mesh | Open source; nuclear physics |

## The resonance-model literature releases data, not solvers

The line of work qscat continues — the nonlocal resonance model and the
exactly-solvable two-dimensional model (Domcke; Čížek, Horáček; Houfek, Rescigno
& McCurdy) — has published extensively without releasing code. In particular no
public repository, Zenodo deposit, or supplementary archive was found for:

- K. Houfek, T. N. Rescigno, C. W. McCurdy, *Numerically solvable model for
  resonant collisions of electrons with diatomic molecules*, Phys. Rev. A **73**,
  032721 (2006), <https://doi.org/10.1103/PhysRevA.73.032721> — the model qscat
  solves;
- M. Váňa, K. Houfek, Phys. Rev. A **95**, 022714 (2017),
  <https://doi.org/10.1103/PhysRevA.95.022714> — the time-dependent formulation
  qscat's TD route implements.

What this community distributes is *cross-section data* (LXCat, the Belgrade
ACol/VAMDC database) rather than the solvers that produced it. The C++/CUDA code
behind the second reference, from which qscat's physics was re-derived, was never
released either — qscat is its first public form.

## What appears to have no released counterpart

Stated as survey results, not as priority claims:

- **The two-dimensional resonant electron–diatomic model as a library.** No
  released code was found that solves the electronic × nuclear driven problem and
  returns vibrational-excitation, dissociative-attachment, or
  dissociative-recombination cross sections.
- **Exact solution and approximation side by side.** qscat computes the exact
  two-dimensional result and the local-complex-potential approximation of the
  same quantity from one model definition, which is what makes the approximation
  testable rather than merely usable.
- **Three time-dependent extractors driven by one propagation.** The
  Tannor–Weeks, fixed-point (Dirac) and flux extractors share a single propagated
  wavepacket, so their disagreement measures discretisation error rather than
  differing dynamics. No released implementation of any of the three for this
  class of problem was found.
- **An automatic FEM-DVR-ECS discretisation tuner.** Grid parameters in this
  field are published as hand-tuned tables; no code was found that derives them
  from the potential and energy range.
- **A complex-symmetric MUMPS backend with symbolic-analysis reuse across an
  energy sweep**, exposed behind a solver-agnostic API with SuperLU as the
  fallback and differential oracle.

## How this survey was done, and what it cannot see

Queries run on 2026-08-16: PyPI project and simple-index lookups; GitHub
repository search for `FEM-DVR`, `exterior complex scaling`, `discrete variable
representation`, `electron molecule scattering`, `dissociative attachment`,
`resonant electron scattering`, `vibrational excitation cross section`; GitHub
code search for `"exterior complex scaling"` and `FEM-DVR`; Zenodo software
search; and literature search for released implementations of the nonlocal
resonance model, the local-complex-potential model, and Tannor–Weeks cross-section
extraction.

Notable negatives: the GitHub repository search returned **no** repositories
whose name or description matches `FEM-DVR`, `dissociative attachment`,
`resonant electron scattering`, or `vibrational excitation cross section`; the
code search for `"exterior complex scaling"` returned essentially only qscat and
`quantumGrid`. No Debian, Fedora, or conda-forge package exists for any code in
this space.

**Limits.** This covers what is publicly indexed. It cannot see unreleased group
codes (the Charles University and LBNL Fortran/C++ codes in this field are, as
far as could be determined, unreleased), supplementary material behind paywalls,
institutional repositories, or code shared privately on request. Treat "no
released implementation found" as a statement about the public record on the
survey date, not about what exists.

## On the name

`qscat` here is *quantum scattering*, and the project's domain is `qscat.org` —
<https://qscat.org> is the project's landing page and `data.qscat.org` serves the
published artifacts. The acronym is used elsewhere for unrelated things — most
visibly the
[QGIS Shoreline Change Analysis Tool](https://github.com/qscat/qscat)
(<https://doi.org/10.1016/j.envsoft.2024.106263>) and, historically, as shorthand
for NASA's QuikSCAT scatterometer mission. There is no functional overlap and no
packaging collision: this project is repo-only and publishes nothing to PyPI,
the PyPI distribution name `qscat` was unregistered as of 2026-08-16, and the
QGIS tool ships through the QGIS plugin repository rather than PyPI. Readers
searching literature or code for "QSCAT" should expect all three.
