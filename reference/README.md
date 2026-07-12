# reference/ — read-only porting oracles

These trees are the prior QSCAT codebases, kept ONLY as sources of algorithms and
math for clean reimplementation. **Do not edit, build, or import them as dependencies.**

- `eMoScat/` — C++/CUDA electron–molecule scattering (FEM-DVR-ECS, Coulomb/Bessel,
  Crank–Nicolson, LCP/NRM). Snapshot copy.
- `libXcuda/` — recovered CUDA layer (submodule). Reference for eventual GPU kernels.

Use the `port-scout` agent to extract a method's algorithm before reimplementing it.
