# The scattering engine

The model-independent solver and the time-independent routes to a cross
section.

```{toctree}
:maxdepth: 1

qscat-core-scattering
n2-resonance
n2-cross-section
n2-2d-cross-section
nrm-vibrational-excitation
```

- {doc}`qscat-core-scattering` — the model-independent engine and the
  model/engine split it enforces.
- {doc}`n2-resonance` — locating the resonance pole.
- {doc}`n2-cross-section` — the one-dimensional time-independent route.
- {doc}`n2-2d-cross-section` — the exact two-dimensional driven solve, gated
  against independent published data.
- {doc}`nrm-vibrational-excitation` — the nonlocal resonance model's
  vibrational-excitation route: two-potential background + resonant
  T-matrices on the shared nonlocal kernel, reproducing the exact solver
  to better than 0.7 % on N₂ and F₂.
