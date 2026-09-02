# Resonances and levels

Quasi-bound states: the Born–Oppenheimer approximation to them, the exact
two-dimensional poles, and how to tell a resonance from an artefact.

```{toctree}
:maxdepth: 1

lcp-resonance-levels
exact-2d-resonances
h2plus-resonance-states
potential-factory
coupled-partial-waves
```

- {doc}`lcp-resonance-levels` — Born–Oppenheimer quasi-bound levels in the
  complex curve.
- {doc}`exact-2d-resonances` — the same levels without the approximation: poles
  of the full 2-D S-matrix, and what the Born–Oppenheimer error actually
  measures on N₂.
- {doc}`h2plus-resonance-states` — the same comparison on H₂⁺, against a σ_DR
  sweep: the Born–Oppenheimer error sorted by regime, and the four "resonances"
  that turned out not to be.
- {doc}`potential-factory` — fitting a richer model surface to a tiered target
  curve: round-tripped against N₂/NO/F₂'s own published parameters, then O₂
  from Alt & Houfek's published curves to its spin–orbit-resolved VE cross
  section on the paper's own nonlocal-model comb ({doc}`../molecules/o2`).
- {doc}`coupled-partial-waves` — the parked
  {doc}`angular-coupled-channels <angular-coupled-channels>` direction,
  delivered: does NO's single-partial-wave shape resonance stay a single pole
  and a good approximation once it is allowed to couple to neighbouring
  partial waves? Yes to both. Only `l = 1` hosts a resonance at all — O⁻ has
  one bound orbital, 2p — which explains the single pole rather than merely
  observing it, and the truncation costs 2–7 % on the angle-integrated VE
  cross section against a reference converged to 0.3–0.5 %. That observable is
  the wrong one for the question anyway, since it sums over the exit partial
  waves the anisotropy produces; the differential cross section has not been
  computed.
