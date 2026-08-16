# qscat.units and qscat.exceptions

## qscat.units

Atomic units are used throughout qscat: energies in Hartree, lengths in Bohr.
These conversions exist so unit handling is never scattered through method
code — there is exactly one place that knows the Hartree/eV factor.

```{eval-rst}
.. currentmodule:: qscat.units

.. autosummary::
   :nosignatures:

   HARTREE_TO_EV
   EV_TO_HARTREE
   hartree_to_ev
   ev_to_hartree

.. autodata:: qscat.units.HARTREE_TO_EV
.. autodata:: qscat.units.EV_TO_HARTREE
.. autofunction:: qscat.units.hartree_to_ev
.. autofunction:: qscat.units.ev_to_hartree
```

## qscat.exceptions

Every recoverable error qscat raises is a subclass of `QscatError`. Each also
subclasses the built-in it replaces, so catching the built-in remains valid.
Generic argument validation may still raise a plain `ValueError`/`TypeError`.

```{eval-rst}
.. currentmodule:: qscat.exceptions

.. autosummary::
   :nosignatures:

   QscatError
   GridError
   ModelError
   BackendError
   ConvergenceError

.. autoexception:: qscat.exceptions.QscatError
   :members:
.. autoexception:: qscat.exceptions.GridError
   :members:
.. autoexception:: qscat.exceptions.ModelError
   :members:
.. autoexception:: qscat.exceptions.BackendError
   :members:
.. autoexception:: qscat.exceptions.ConvergenceError
   :members:
```
