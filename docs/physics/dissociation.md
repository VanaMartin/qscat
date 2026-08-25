# Dissociation and approximations

The channels where the molecule comes apart, and the two reduced models —
local and nonlocal — measured against the exact solver.

```{toctree}
:maxdepth: 1

diatomic-ve-cross-sections
nonlocal-resonance-model
h2plus-dr
```

- {doc}`diatomic-ve-cross-sections` — NO and F₂, and the local-complex-potential
  approximation measured against the exact oracle.
- {doc}`nonlocal-resonance-model` — the rung above the LCP: a nonlocal,
  energy-dependent kernel that keeps the energy dependence and the nonlocality
  the LCP throws away. On **vibrational excitation** it reproduces the exact
  oracle to better than 0.7 % on both N₂ and F₂, elastic and first-inelastic
  alike, and the note argues that is physics rather than luck — an R-independent
  discrete state carries no derivative couplings, so the model is formally exact
  and the residual is discretization error. On **dissociative attachment** it
  reproduces the oracle on F₂ and on NO; the five-to-eight-order NO "collapse"
  this page used to report was the exact 2-D oracle's own volume-T-matrix defect,
  charged to the model, and is retracted.
- {doc}`h2plus-dr` — dissociative recombination for an ionic target.
