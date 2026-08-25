"""`n2_electronic_grid` has ONE definition, and this package re-exports it.

The factory used to be hand-duplicated: the toy-model
`projects/n2_resonance/grid_n2.py` (validated in
`projects/n2_resonance/test_grid_n2.py`) and a byte-for-byte copy in
`validation/n2/resonance.py`, guarded by a test that compared the two grids
they built. The copy is gone -- `resonance.py` imports the toy-model factory,
which is the allowed direction (`validation/` may import `projects/`; not the
reverse). What remains is an identity check: if someone re-introduces a local
copy, `resonance.n2_electronic_grid` stops BEING the toy-model function and
this fails immediately, without a grid comparison that could pass while the
two drifted in some property it did not compare.
"""

from __future__ import annotations

from projects.n2_resonance import grid_n2 as ref_grid_n2
from validation.n2 import resonance


def test_resonance_reexports_the_toy_model_grid_factory():
    assert resonance.n2_electronic_grid is ref_grid_n2.n2_electronic_grid
