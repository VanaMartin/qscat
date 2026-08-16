# API reference

Generated from the source by autodoc, one page per subpackage.

## What is public

The public API is exactly the names exported in each submodule's `__all__`.
Anything with a leading underscore, and any name not in an `__all__`, is
private and may change without notice. While the version is `0.y.z`, minor
releases may contain breaking changes to the public API; they are called out in
`CHANGELOG.md`. See ADR 0004 (`docs/adr/0004-public-api-stability-policy.md`)
for the full policy, including the *provisional* marker and the post-1.0
deprecation rule.

Every recoverable error is a subclass of `qscat.exceptions.QscatError`, so
`except QscatError` is a valid catch-all.

## The layers

| layer | modules | what lives here |
|---|---|---|
| Engine | {doc}`core`, {doc}`model` | The model-independent scattering engine, and the per-molecule potentials it consumes. |
| Numerics | {doc}`dvr`, {doc}`ecs`, {doc}`linalg`, {doc}`evolution`, {doc}`special` | FEM-DVR grids, the complex-scaling map, sparse linear algebra, time propagators, radial/Coulomb special functions. |
| Tooling | {doc}`tuning`, {doc}`viz` | Deriving a grid from a potential; rendering wavefunctions. |
| Base | {doc}`base` | Unit conversions and the exception hierarchy. |

`qscat.core` never imports `qscat.model` at runtime — it depends only on the
`ResonanceModel` protocol. Adding a molecule is a registry entry, never solver
code.

```{toctree}
:maxdepth: 1

core
model
dvr
ecs
linalg
evolution
special
tuning
viz
base
```
