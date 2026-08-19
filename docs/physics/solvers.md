# Linear algebra and solvers

The sparse factorizations and eigensolves the two-dimensional problems run
on, and the reuse tricks that make a sweep affordable.

```{toctree}
:maxdepth: 1

mumps-sparse-backend
ti-energy-sweep-reuse
shift-invert-eigensolver
```

- {doc}`mumps-sparse-backend` — the complex-symmetric MUMPS backend.
- {doc}`ti-energy-sweep-reuse` — reusing the symbolic factorization across an
  energy sweep.
- {doc}`shift-invert-eigensolver` — the eigenpairs nearest a complex shift, for
  resonances that sit in the interior of the spectrum. Validated in 1-D only.
