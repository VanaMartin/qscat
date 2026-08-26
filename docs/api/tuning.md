# qscat.tuning

The automatic discretisation tuner: derive a minimal FEM-DVR-ECS grid at a
target precision from the potential and energy range, instead of hand-picking
element lengths. See `docs/physics/discretisation-tuning.md`, including the
documented limits of the 1-D convergence probes.

```{eval-rst}
.. currentmodule:: qscat.tuning

.. autosummary::
   :nosignatures:

   IncidentSpec
   PotentialProfile
   ProbeResult
   analyze_potential
   equidistribution_elements
   grid_cost
   interaction_region
   max_stable_angle
   optimal_real_mesh
   order_for_wavenumber
   probe_channel_representation
   probe_electronic
   probe_nuclear
   propose_grid
   refine
   refine_elements_in_window
   refine_to_2d_convergence
   required_extent
   resonance_curve_arrays
   tensor_cost
   tune_ecs_tail
   tw_analysis
```

```{eval-rst}
.. automodule:: qscat.tuning
   :members:
   :imported-members:
```
