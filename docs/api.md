# API reference

Generated from the source by autodoc. The public surface of each module is its
`__all__` (see ADR 0004 for the stability policy).

## Exceptions

```{eval-rst}
.. automodule:: qscat.exceptions
   :members:
```

## qscat.core

The model-independent scattering engine. `ScatteringProblem` is the recommended
high-level entry point; the functional solvers below are the low-level layer.

```{eval-rst}
.. autoclass:: qscat.core.ScatteringProblem
   :members:

.. autofunction:: qscat.core.ve_cross_section
.. autofunction:: qscat.core.da_cross_section
.. autofunction:: qscat.core.dr_cross_section
.. autofunction:: qscat.core.td_ve_cross_section
.. autofunction:: qscat.core.td_da_cross_section
.. autofunction:: qscat.core.vibrational_states
```

## qscat.model

```{eval-rst}
.. automodule:: qscat.model
   :members:
   :imported-members:
```

## qscat.dvr

```{eval-rst}
.. automodule:: qscat.dvr
   :members:
   :imported-members:
```

## qscat.linalg

```{eval-rst}
.. automodule:: qscat.linalg
   :members:
   :imported-members:
```

## qscat.ecs

```{eval-rst}
.. automodule:: qscat.ecs
   :members:
   :imported-members:
```

## qscat.evolution

```{eval-rst}
.. automodule:: qscat.evolution
   :members:
   :imported-members:
```

## qscat.special

```{eval-rst}
.. automodule:: qscat.special
   :members:
   :imported-members:
```

## qscat.tuning

```{eval-rst}
.. automodule:: qscat.tuning
   :members:
```
