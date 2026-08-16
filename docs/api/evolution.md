# qscat.evolution

Time propagators for `d/dt psi = -i H psi` with complex, possibly
non-Hermitian `H`. The order-N diagonal-Padé stepper generalizes
Crank–Nicolson (order 1); order 3 is what makes the time-dependent cross
sections converge.

```{eval-rst}
.. currentmodule:: qscat.evolution

.. autosummary::
   :nosignatures:

   make_cn_stepper
   make_pade_stepper
   make_sparse_cn_stepper
   pade_roots
```

```{eval-rst}
.. automodule:: qscat.evolution
   :members:
   :imported-members:
```
