# projects/femdvr_ecs

Development toy model for the FEM-DVR-ECS radial grid (lifecycle stages 1-3:
design, toy model, validate). The validated code here has been **promoted**
to the qscat standard library:

- `gll.py`, `spec.py`, `grid.py`, `kinetic.py`, `operators.py` -> `qscat.dvr`
  (`libs/qscat/qscat/dvr/`)
- the ECS coordinate map (`R0 + (x - R0) e^{i theta}`) -> `qscat.ecs.ecs_map`
  (`libs/qscat/qscat/ecs/`)

This directory remains in place as the origin/dev copy and its own test
suite (`test_gll.py`, `test_grid.py`, `test_kinetic_benchmarks.py`,
`test_ecs_benchmarks.py`) stays green independently of `libs/qscat`. New
consumers should import from `qscat.dvr` / `qscat.ecs`, not from here.

See `docs/physics/femdvr-ecs.md` for the method and validation benchmarks,
and `.superpowers/sdd/femdvr-ecs-extraction.md` for the port-scout extraction
from `reference/eMoScat` this implementation is based on.
