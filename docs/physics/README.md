# Theory notes

One note per method: the derivation, the equations, the unit conventions
(atomic units throughout), the validation evidence, and the literature it
comes from. These are the working notes behind the implementation, so they
record limitations and negative results as well as what works.

The notes are grouped into seven sections in the sidebar:

- {doc}`discretisation` — turning a coordinate into a matrix, and choosing the
  grid.
- {doc}`solvers` — the sparse factorizations and eigensolves underneath.
- {doc}`engine` — the model-independent solver and the time-independent routes.
- {doc}`time-dependent` — the wavepacket route to the same answers.
- {doc}`dissociation` — the channels where the molecule comes apart, and the
  reduced models measured against the exact solver.
- {doc}`resonances` — quasi-bound states, exactly and approximately.
- {doc}`open-directions` — designed but parked.

{doc}`validation-harnesses` cuts across all of them: what each harness in
`validation/` gates, and how to run it.

For a per-molecule view — what has been computed on N₂, NO, F₂ and H₂⁺, and
which notes report it — start from the Molecules section of the sidebar.
