# qscat.dvr

The FEM-DVR grid with an exterior-complex-scaled tail: grid construction,
kinetic-energy assembly (dense and sparse), and the N-dimensional tensor
layer. See `docs/physics/femdvr-ecs.md` and
`docs/physics/nd-tensor-hamiltonian.md`.

```{eval-rst}
.. currentmodule:: qscat.dvr

.. autosummary::
   :nosignatures:

   ElementSpec
   GridSpec
   FemDvrEcsGrid
   kinetic
   kinetic_sparse
   dvr_first_derivative_at_node
   dvr_interpolation_matrix
   hamiltonian
   eigen
   gll_nodes_weights
   diff_matrix
   TensorGrid
   kinetic_nd
   potential_nd
   hamiltonian_nd
```

```{eval-rst}
.. automodule:: qscat.dvr
   :members:
   :imported-members:
```
