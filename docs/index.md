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
api/index
related-work
```

```{toctree}
:maxdepth: 2
:caption: Theory — technical

physics/README
physics/validation-harnesses
physics/discretisation
physics/solvers
physics/engine
physics/time-dependent
physics/dissociation
physics/resonances
physics/open-directions
```

```{toctree}
:maxdepth: 1
:caption: Theory — molecules

molecules/n2
molecules/no-f2
molecules/o2
molecules/h2plus
```

## Quick links

- **Install & first cross section:** {doc}`getting-started`
- **API reference:** {doc}`api/index`
- **How this relates to existing codes:** {doc}`related-work`
- **Theory notes:** {doc}`physics/README`
- **Design decisions:** the `docs/adr/` directory (ADRs)

## Citing

If you use qscat in research, citation is required — see `CITATION.cff` in the
repository and the "Citing qscat" section of the package README.

## Indices

- {ref}`genindex`
- {ref}`modindex`
