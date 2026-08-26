# qscat.linalg

Dimension-general sparse linear algebra: Kronecker sums over arbitrary D, a
cached sparse LU with a SuperLU and a complex-symmetric MUMPS backend, and the
bilinear (non-conjugated) ECS inner product. See
`docs/physics/mumps-sparse-backend.md`.

```{eval-rst}
.. currentmodule:: qscat.linalg

.. autosummary::
   :nosignatures:

   kron_sum
   c_product
   SparseLU
   Ordering
   ShiftInvertEigs
   default_backend
   set_default_backend
   get_default_backend
```

```{eval-rst}
.. automodule:: qscat.linalg
   :members:
   :imported-members:
```
