# qscat

**Validated quantum-scattering numerics for electron–diatomic collisions.**

qscat is a CPU-first Python library for electron–molecule scattering with the
FEM-DVR method and exterior complex scaling (ECS): vibrational-excitation (VE),
dissociative-attachment (DA), and dissociative-recombination (DR) cross sections,
by time-independent and time-dependent routes. Everything is in atomic units, and
every method is validated against an analytic benchmark, a conservation law, a
convergence study, or an independent reference.

```{toctree}
:maxdepth: 2
:caption: Contents

getting-started
api
```

## Quick links

- **Install & first cross section:** {doc}`getting-started`
- **API reference:** {doc}`api`
- **Theory notes:** the `docs/physics/` directory in the repository
- **Design decisions:** the `docs/adr/` directory (ADRs)

## Citing

If you use qscat in research, citation is required — see `CITATION.cff` in the
repository and the "Citing qscat" section of the package README.

## Indices

- {ref}`genindex`
- {ref}`modindex`
