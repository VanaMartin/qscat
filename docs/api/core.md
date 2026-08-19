# qscat.core

The model-independent electron–diatomic scattering engine. `ScatteringProblem`
is the recommended entry point: it bundles the grid, model, and vibrational
basis once and exposes every observable as a method. The functional solvers
below are the low-level layer those methods call.

This module never imports `qscat.model` at runtime — it depends only on the
`ResonanceModel` protocol, so a new molecule needs no change here.

## The problem object

```{eval-rst}
.. autoclass:: qscat.core.ScatteringProblem
   :members:
```

## Cross sections (time-independent)

```{eval-rst}
.. autofunction:: qscat.core.ve_cross_section
.. autofunction:: qscat.core.da_cross_section
.. autofunction:: qscat.core.dr_cross_section
```

## Cross sections (time-dependent)

```{eval-rst}
.. autofunction:: qscat.core.td_ve_cross_section
.. autofunction:: qscat.core.td_ve_cross_sections_all
.. autofunction:: qscat.core.td_da_cross_section
.. autofunction:: qscat.core.td_da_cross_sections_all
```

## Grids

```{eval-rst}
.. autofunction:: qscat.core.electronic_grid
.. autofunction:: qscat.core.nuclear_grid
.. autofunction:: qscat.core.fem_grid_exp_tail
.. autofunction:: qscat.core.segmented_grid
.. autofunction:: qscat.core.ecs_angle_family
.. autofunction:: qscat.core.assert_shared_real_nodes
```

## Channels

```{eval-rst}
.. autofunction:: qscat.core.channel_vector
.. autofunction:: qscat.core.anion_electronic_states
.. autofunction:: qscat.core.v_dr_diag
.. autofunction:: qscat.core.outgoing_channel
.. autofunction:: qscat.core.outgoing_channel_nuclear
```

## Vibrational structure

```{eval-rst}
.. autofunction:: qscat.core.vibrational_states
.. autoclass:: qscat.core.VibrationalBasis
   :members:
```

## Wavepacket and correlation

```{eval-rst}
.. autofunction:: qscat.core.gaussian_coeffs
.. autofunction:: qscat.core.initial_state
.. autofunction:: qscat.core.eta_incident
.. autofunction:: qscat.core.eta_outgoing
.. autofunction:: qscat.core.hankel_point_value
.. autofunction:: qscat.core.outgoing_surface_wave
.. autofunction:: qscat.core.propagate
.. autofunction:: qscat.core.sigma_from_correlations
```

## Time-dependent energy extractors

All three share one propagate-once protocol, so a single propagation can drive
every extractor. See `docs/physics/td-extractors.md`.

```{eval-rst}
.. autoclass:: qscat.core.Extractor
   :members:
.. autoclass:: qscat.core.TannorWeeks
   :members:
.. autoclass:: qscat.core.Dirac
   :members:
.. autoclass:: qscat.core.Flux
   :members:
```

## The LCP approximation

The local-complex-potential reduction and the Born–Oppenheimer resonance
levels built on it. These are the *approximation* under test against the exact
solvers above; see `docs/physics/diatomic-ve-cross-sections.md` and
`docs/physics/lcp-resonance-levels.md`.

```{eval-rst}
.. autofunction:: qscat.core.local_complex_potential
.. autofunction:: qscat.core.lcp_da_cross_section
.. autofunction:: qscat.core.resonance_levels
.. autofunction:: qscat.core.lcp_resonance_levels
.. autoclass:: qscat.core.ResonanceLevels
   :members:
```

## Exact resonance states

Poles of the full 2-D S-matrix, with no Born–Oppenheimer separation and no
local approximation — the objects `resonance_levels` above approximates.
Identification is by ECS angle stability in *both* coordinates; see
`docs/physics/exact-2d-resonances.md`.

```{eval-rst}
.. autofunction:: qscat.core.exact_resonance_states
.. autoclass:: qscat.core.ExactResonanceStates
   :members:
```

## Born–Oppenheimer reference states

The picture an exact pole departs from: electronic eigen-curves over a nuclear
grid, and the product states `φ_j(r; R)·χ_v(R)` built on them.
`electronic_curves` covers the bound (ion / Rydberg) case and `resonance_curve`
the resonant (neutral / anion) one, so one comparator serves both. See
`docs/physics/h2plus-resonance-states.md`.

```{eval-rst}
.. autofunction:: qscat.core.electronic_curves
.. autofunction:: qscat.core.resonance_curve
.. autofunction:: qscat.core.bo_basis
.. autofunction:: qscat.core.bo_basis_from_levels
.. autofunction:: qscat.core.n_eff
.. autofunction:: qscat.core.admissible_levels
.. autofunction:: qscat.core.basis_covers
.. autoclass:: qscat.core.ElectronicCurves
   :members:
.. autoclass:: qscat.core.BoState
   :members:
.. autoclass:: qscat.core.BoBasis
   :members:
```

## Assignment and verification

Is a pole a resonance at all, and which quasi-bound state is it? Angle
stability is necessary and not sufficient; overlap against a BO basis is the
check that catches a rotated-continuum eigenvalue sitting in a stable corner.

```{eval-rst}
.. autofunction:: qscat.core.overlap
.. autofunction:: qscat.core.real_weight
.. autofunction:: qscat.core.pair_by_overlap
.. autofunction:: qscat.core.pair_one_to_one
.. autofunction:: qscat.core.peak_positions
.. autofunction:: qscat.core.peak_alignment
.. autoclass:: qscat.core.OverlapPair
   :members:
.. autoclass:: qscat.core.PeakAlignment
   :members:
```

## Nonlocal resonance model

`qscat.core.nrm` -- the Feshbach-projection nonlocal resonance model (NRM) for
dissociative attachment (Houfek, Rescigno & McCurdy, Phys. Rev. A 77, 012710
(2008)). It sits between the LCP approximation above and the exact 2-D
solvers, and is kept as a subpackage that `qscat.core` does not import at
runtime (see the hard boundary above) -- import it explicitly, e.g.
`from qscat.core.nrm import nrm_da_cross_section`. See
`docs/physics/nonlocal-resonance-model.md`.

```{eval-rst}
.. autoclass:: qscat.core.nrm.DiscreteState
   :members:
.. autoclass:: qscat.core.nrm.AsymptoticDiscreteState
   :members:
.. autoclass:: qscat.core.nrm.PhysicalDiscreteState
   :members:
.. autofunction:: qscat.core.nrm.nrm_ingredients
.. autoclass:: qscat.core.nrm.NrmIngredients
   :members:
.. autofunction:: qscat.core.nrm.nrm_da_cross_section
.. autofunction:: qscat.core.nrm.nrm_ve_cross_section
```

## Plotting

```{eval-rst}
.. autofunction:: qscat.core.plot_cross_sections
.. autofunction:: qscat.core.plot_resonance_levels
.. autoclass:: qscat.core.ComparisonPanel
   :members:
.. autofunction:: qscat.core.plot_route_comparison
```
