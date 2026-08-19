# Discretisation

How a continuous radial coordinate becomes a finite matrix, and how the grid
is chosen rather than guessed.

```{toctree}
:maxdepth: 1

femdvr-ecs
nd-tensor-hamiltonian
discretisation-tuning
```

- {doc}`femdvr-ecs` — the FEM-DVR grid with an exterior-complex-scaled tail,
  and the four analytic benchmarks that pin it down.
- {doc}`nd-tensor-hamiltonian` — the N-dimensional sparse tensor Hamiltonian.
- {doc}`discretisation-tuning` — deriving a grid from the potential instead of
  hand-picking element lengths, and where the 1-D probes are not sufficient.
