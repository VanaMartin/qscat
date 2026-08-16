# qscat.model

Everything tied to a specific molecule: the `ResonanceModel` protocol that
`qscat.core` depends on, the shared neutral and ionic potential forms, and the
molecule registry. Adding a molecule is a registry entry plus validation,
never solver code.

```{eval-rst}
.. currentmodule:: qscat.model

.. autosummary::
   :nosignatures:

   ResonanceModel
   DiatomicResonanceModel
   IonicResonanceModel
   N2
   NO
   F2
   H2P
```

```{eval-rst}
.. automodule:: qscat.model
   :members:
   :imported-members:
```
