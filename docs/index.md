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
:maxdepth: 1
:caption: Theory

physics/README
physics/femdvr-ecs
physics/nd-tensor-hamiltonian
physics/discretisation-tuning
physics/mumps-sparse-backend
physics/ti-energy-sweep-reuse
physics/shift-invert-eigensolver
physics/qscat-core-scattering
physics/n2-resonance
physics/n2-cross-section
physics/n2-2d-cross-section
physics/n2-td-cross-section
physics/n2-2d-td-cross-section
physics/td-extractors
physics/td-da
physics/diatomic-ve-cross-sections
physics/nonlocal-resonance-model
physics/h2plus-dr
physics/lcp-resonance-levels
physics/exact-2d-resonances
physics/h2plus-resonance-states
physics/angular-coupled-channels
physics/optimization-targets
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
