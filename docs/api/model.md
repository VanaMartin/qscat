# qscat.model

Everything tied to a specific molecule: the `ResonanceModel` protocol that
`qscat.core` depends on, the shared neutral and ionic potential forms, the
potential factory's flexible form (`FlexibleDiatomicModel`, with `SmoothR`
sigmoid or long-range-correct `TailR` coefficient functions, embedding the
published models exactly via `from_diatomic`), and the molecule registry —
`O2` being the first entry that is a fit rather than a published parameter
set. Adding a molecule is a registry entry plus validation, never solver code.

```{eval-rst}
.. currentmodule:: qscat.model

.. autosummary::
   :nosignatures:

   ResonanceModel
   DiatomicResonanceModel
   IonicResonanceModel
   FlexibleDiatomicModel
   SmoothR
   TailR
   from_diatomic
   N2
   NO
   F2
   H2P
   O2
```

```{eval-rst}
.. automodule:: qscat.model
   :members:
   :imported-members:
```
