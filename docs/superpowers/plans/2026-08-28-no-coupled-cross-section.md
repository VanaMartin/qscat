# Coupled Cross Section Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure what the fixed-$l$ reduction costs in an *observable* —
$\sigma_{\rm VE}(E)$ from the exact 2-D driven equation with the partial waves
coupled, against the same quantity with one partial wave, on the same deck.

**Architecture:** The coupled interaction is a sparse block matrix rather than a
diagonal, so the driven solve becomes `psi + lu.solve(V @ psi)` where `V` is
that matrix. The entrance channel vector is the existing single-channel one
embedded in one block; each exit projection reads one block of `V @ psi_plus`.
Everything else — the post-form T-matrix, the `4 pi^3 |T|^2 / 2E`
normalisation, the `SparseLU.refactor` sweep — carries over unchanged from
`qscat.core.driven`.

**Tech Stack:** Python 3.12+, numpy, scipy.sparse, `qscat.core.driven` /
`channels` / `vibrational`, `qscat.linalg.SparseLU` (MUMPS backend for the
production run), pytest, uv, Docker on sadaharu.

**Spec:** `docs/superpowers/specs/2026-08-28-no-coupled-cross-section-design.md`

## Global Constraints

- **Atomic units throughout.** Energies in Hartree, lengths in Bohr, cross
  sections in bohr².
- **Complex-safe.** `r` and `R` may be complex on the ECS tails; never coerce
  to a real dtype. Every matrix is **complex symmetric, never Hermitian** — no
  `.conj()`, no Hermitian-only routine.
- **The comparison is differential**: the coupled model (`n_channels > 1`) and
  the fixed-$l$ model (`n_channels = 1`) must run through the SAME functions on
  the SAME deck, the SAME energy mesh and the SAME entrance channel. Never a
  different grid, mesh or tolerance for one branch.
- **Entrance is $l = \Lambda = 1$ for both models** — channel index 0, since
  `channel_ells()` starts at `Lambda`.
- **Package-absolute imports only.** `validation/` may import `projects/` and
  `qscat`; `projects/` must NEVER import `validation/`. Nothing may enter
  `qscat` in this plan.
- `E741` is already ignored for `projects/no_coupled_channels/*` and
  `validation/coupled/*` in `pyproject.toml` — `l` is the partial-wave angular
  momentum. Do not add inline `# noqa`, do not rename physics variables.
- `__all__` must stay sorted (ruff RUF022).
- **Never `git commit -a`.** Stage explicit paths.
- Commit trailers on every commit:
  ```
  Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt
  ```
- **Do not run the production sweep on this machine.** Measured: SuperLU
  factorises the $N_l=2$ deck in 208 s with no analysis reuse. Tasks 1–6 use a
  small test grid; Task 7 runs on sadaharu.

## File Structure

| File | Responsibility |
|---|---|
| `projects/no_coupled_channels/model.py` (modify) | `CoupledModel.interaction_matrix` — the coupled $V_{ll'}$ alone, as a sparse block matrix |
| `projects/no_coupled_channels/scattering.py` (create) | `coupled_channel_vector`, `channel_block`, `CoupledSigma`, `coupled_ve_cross_section` |
| `projects/no_coupled_channels/test_scattering.py` (create) | the identity gates: s = 0 against the shipped solver, kappa = 0 by parity |
| `validation/coupled/energies.py` (create) | the threshold-aware mesh |
| `validation/coupled/test_energies.py` (create) | mesh shape, threshold bracketing, tolerance dedup |
| `validation/coupled/cross_section.py` (create) | the campaign: timing probe, the three sweeps, the report |
| `validation/coupled/test_cross_section.py` (create) | campaign structure on a tiny grid |
| `validation/coupled/figures.py` (modify) | the sigma(E) comparison panel |
| `docs/physics/coupled-partial-waves.md` (modify) | the cross-section result joins the existing note |
| `CLAUDE.md` (modify) | the repo-map entry gains the cross-section route |

---

### Task 1: The coupled interaction matrix

**Files:**
- Modify: `projects/no_coupled_channels/model.py`
- Test: `projects/no_coupled_channels/test_model.py`

**Interfaces:**
- Consumes: `TwoCentreWell.v_block`, `assemble_coupled`, `CoupledModel._coupling_table` (all existing).
- Produces: `CoupledModel.interaction_matrix(tgrid: TensorGrid) -> sp.csr_matrix` — the coupled interaction $V_{ll'}$ on the tensor grid, with NO kinetic energy, NO `v0(R)` and NO centrifugal term. Shape `(n_channels * tgrid.size, n_channels * tgrid.size)`.

**Why this exists:** the single-channel driven solve multiplies by
`model.interaction_diag(tgrid)`, a flat array, because the interaction is
diagonal. Coupled, it is not: the off-diagonal blocks are exactly the coupling.
So the solve needs a matrix, and it must contain the perturbation ALONE —
`v0(R)` and the centrifugal term belong to the free channel Hamiltonian that
`channel_vector`'s $F_{E,l}(r)\chi_v(R)$ already solves.

- [ ] **Step 1: Write the failing test**

Append to `projects/no_coupled_channels/test_model.py`:

```python
def test_interaction_matrix_is_the_perturbation_alone() -> None:
    """At s = 0 the coupled interaction must be block-diagonal, and its l = 1
    block must be exactly the shipped model's interaction -- v0 and the
    centrifugal term belong to the free Hamiltonian, not to the perturbation."""
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    tg = _tensor_grid()
    n = tg.size
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=3).interaction_matrix(tg))
    assert V.shape == (3 * n, 3 * n)

    # RELATIVE to the interaction's own magnitude, not absolute. The
    # off-diagonal blocks vanish by quadrature exactness, so what survives is
    # round-off on the potential -- measured, 2.9e-15 against a magnitude of
    # 6.0, i.e. 4.8e-16 relative. An absolute bound would silently tighten or
    # loosen with lambda(R), and on a weaker interaction it would start
    # rejecting a correct implementation.
    scale = float(np.max(np.abs(NO.interaction_diag(tg))))
    for i in range(3):
        for j in range(3):
            if i != j:
                blk = _block(V, i, j, n)
                assert blk.nnz == 0 or float(np.max(np.abs(blk.data))) < 1e-13 * scale

    got = np.asarray(_block(V, 0, 0, n).diagonal())
    np.testing.assert_allclose(got, NO.interaction_diag(tg), rtol=0, atol=1e-14)


def test_interaction_matrix_is_complex_symmetric() -> None:
    """Every operator here is complex symmetric, never Hermitian -- the ECS
    contour makes it so, and a Hermitian-only routine downstream would be
    silently wrong rather than loudly."""
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    tg = _tensor_grid()
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=3).interaction_matrix(tg))
    diff = (V - V.T).tocoo()
    assert diff.nnz == 0 or float(np.max(np.abs(diff.data))) < 1e-14


def test_interaction_matrix_couples_once_the_anisotropy_is_on() -> None:
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.3)
    tg = _tensor_grid()
    n = tg.size
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=2).interaction_matrix(tg))
    off = _block(V, 0, 1, n)
    assert off.nnz > 0
    assert float(np.max(np.abs(off.data))) > 1e-6
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest projects/no_coupled_channels/test_model.py -k interaction_matrix -v`
Expected: FAIL — `AttributeError: 'CoupledModel' object has no attribute 'interaction_matrix'`

- [ ] **Step 3: Write the implementation**

Add to `CoupledModel` in `projects/no_coupled_channels/model.py`:

```python
    def interaction_matrix(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """The coupled interaction `V_{ll'}` on `tgrid`, sparse, channel-outermost.

        The PERTURBATION alone: no kinetic energy, no `v0(R)` and no
        centrifugal term. Those belong to the free channel Hamiltonian, whose
        solutions `channel_vector` already supplies -- putting them here would
        drive the Lippmann-Schwinger equation with the wrong operator and
        silently produce a plausible wrong T-matrix.

        The single-channel sibling is `DiagonalChannelModel.interaction_diag`,
        which returns a flat array because one channel's interaction really is
        diagonal. Coupled it is not: the off-diagonal blocks ARE the coupling.
        """
        r, R = tgrid.points()
        diagonal = [
            sp.diags(
                np.asarray(self.well.v_block(l, l, r, R), dtype=np.complex128).ravel(),
                format="csr",
            )
            for l in self.channel_ells()
        ]
        return assemble_coupled(diagonal, self._coupling_table(r, R))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest projects/no_coupled_channels/test_model.py -v`
Expected: PASS (10 tests — the 7 existing plus these 3)

- [ ] **Step 5: Lint and commit**

```bash
uv run ruff check projects/no_coupled_channels
uv run ruff format projects/no_coupled_channels
git add projects/no_coupled_channels/model.py projects/no_coupled_channels/test_model.py
git commit -m "feat(coupled): the coupled interaction as a matrix, not a diagonal

The single-channel driven solve multiplies by a flat array because one
channel's interaction is diagonal. Coupled it is not -- the off-diagonal
blocks ARE the coupling -- so the solve needs a matrix. It carries the
perturbation alone: v0 and the centrifugal term belong to the free channel
Hamiltonian whose solutions channel_vector supplies, and putting them here
would drive Lippmann-Schwinger with the wrong operator.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

### Task 2: Channel vectors that know their block

**Files:**
- Create: `projects/no_coupled_channels/scattering.py`
- Test: `projects/no_coupled_channels/test_scattering.py`

**Interfaces:**
- Consumes: `qscat.core.channels.channel_vector`, `CoupledModel.channel_ells()`.
- Produces:
  - `coupled_channel_vector(tgrid: TensorGrid, k: float, chi_v: NDArray[complex128], ells: tuple[int, ...], channel: int, *, charge: int = 0) -> NDArray[complex128]` — length `len(ells) * tgrid.size`, the single-channel vector in block `channel`, zeros elsewhere.
  - `channel_block(vec: NDArray[complex128], channel: int, n: int) -> NDArray[complex128]` — the `channel`-th length-`n` slice of a channel-outermost vector.

**Physics the implementer needs:** $k=\sqrt{2(E-\varepsilon_v)}$ depends on the
vibrational level alone. The partial wave changes the Bessel *order* in
$F_{E,l}(r)$, not the momentum — so every channel at a given $v$ shares one
$k$, and `ells[channel]` is what selects the order.

- [ ] **Step 1: Write the failing test**

```python
# projects/no_coupled_channels/test_scattering.py
"""Channel vectors in a channel-outermost layout."""

from __future__ import annotations

import numpy as np
from qscat.core.channels import channel_vector
from qscat.core.grids import segmented_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import NO

from projects.no_coupled_channels.scattering import channel_block, coupled_channel_vector


def _tensor_grid() -> TensorGrid:
    """A deliberately small 2-D grid: these tests are about layout, not physics."""
    el = segmented_grid(((4, 8.0),), ((2, 12.0),), angle_deg=35.0, quadrature=6)
    nu = segmented_grid(((3, 4.0),), ((2, 6.0),), angle_deg=30.0, quadrature=6, x_min=1.0)
    return TensorGrid([el, nu])


def test_the_vector_lands_in_its_own_block_and_nowhere_else() -> None:
    tg = _tensor_grid()
    n = tg.size
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    ells = (1, 2, 3)
    for c in range(3):
        vec = coupled_channel_vector(tg, 0.3, chi[0], ells, c)
        assert vec.shape == (3 * n,)
        np.testing.assert_allclose(
            channel_block(vec, c, n), channel_vector(tg, 0.3, chi[0], ells[c])
        )
        for other in range(3):
            if other != c:
                assert np.all(channel_block(vec, other, n) == 0.0)


def test_one_channel_reproduces_the_shipped_vector_exactly() -> None:
    """With a single channel the coupled vector IS the shipped one -- the
    layout must add nothing at n_channels = 1."""
    tg = _tensor_grid()
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    got = coupled_channel_vector(tg, 0.3, chi[0], (1,), 0)
    np.testing.assert_array_equal(got, channel_vector(tg, 0.3, chi[0], 1))


def test_the_partial_wave_selects_the_bessel_order() -> None:
    """Same k, different l -- the momentum is shared, the order is not."""
    tg = _tensor_grid()
    n = tg.size
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    ells = (1, 2)
    a = channel_block(coupled_channel_vector(tg, 0.3, chi[0], ells, 0), 0, n)
    b = channel_block(coupled_channel_vector(tg, 0.3, chi[0], ells, 1), 1, n)
    assert float(np.max(np.abs(a - b))) > 1e-6
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest projects/no_coupled_channels/test_scattering.py -v`
Expected: FAIL — `ModuleNotFoundError: ...scattering`

- [ ] **Step 3: Write the implementation**

```python
# projects/no_coupled_channels/scattering.py
"""The coupled-channel driven solve: sigma_VE(E) with the partial waves mixed.

Mirrors `qscat.core.driven.ve_cross_section`. Three things differ, and only
three:

- the interaction is a MATRIX (`CoupledModel.interaction_matrix`), so the
  Lippmann-Schwinger step is `psi_i + lu.solve(V @ psi_i)` rather than an
  elementwise product;
- the entrance is a single-channel vector embedded in one block;
- the exit is summed over blocks, because the coupling lets the electron leave
  in a partial wave it did not enter on.

The post-form T-matrix, the non-conjugated c-product and the
`4 pi^3 |T|^2 / 2E` normalisation are unchanged.

This duplicates roughly forty lines of sweep boilerplate from `driven.py`. That
is deliberate at toy stage -- generalising a shipped solver before the coupled
shape has been used twice is what the lifecycle exists to prevent -- and
`test_scattering.py` gates the duplicate against the original at s = 0.

See docs/physics/coupled-partial-waves.md.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.channels import channel_vector
from qscat.dvr import TensorGrid

__all__ = ["channel_block", "coupled_channel_vector"]


def channel_block(
    vec: npt.NDArray[np.complex128], channel: int, n: int
) -> npt.NDArray[np.complex128]:
    """The `channel`-th length-`n` slice of a channel-outermost vector."""
    return np.asarray(vec[channel * n : (channel + 1) * n], dtype=np.complex128)


def coupled_channel_vector(
    tgrid: TensorGrid,
    k: float,
    chi_v: npt.NDArray[np.complex128],
    ells: tuple[int, ...],
    channel: int,
    *,
    charge: int = 0,
) -> npt.NDArray[np.complex128]:
    """`F_{E,l}(r) chi_v(R)` for one partial wave, embedded in its block.

    `k` is shared across channels: it depends on the vibrational level alone,
    since the partial wave changes the Bessel ORDER in `F_{E,l}` rather than
    the momentum. `ells[channel]` is what selects that order.
    """
    n = tgrid.size
    out = np.zeros(len(ells) * n, dtype=np.complex128)
    out[channel * n : (channel + 1) * n] = channel_vector(
        tgrid, k, chi_v, ells[channel], charge=charge
    )
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest projects/no_coupled_channels/test_scattering.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Lint and commit**

```bash
uv run ruff check projects/no_coupled_channels
uv run ruff format projects/no_coupled_channels
git add projects/no_coupled_channels/scattering.py projects/no_coupled_channels/test_scattering.py
git commit -m "feat(coupled): channel vectors in a channel-outermost layout

The entrance is a single-channel vector embedded in one block; k is shared
across channels because the partial wave changes the Bessel order rather than
the momentum.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

### Task 3: The coupled driven solve, gated against the shipped one

**Files:**
- Modify: `projects/no_coupled_channels/scattering.py`
- Test: `projects/no_coupled_channels/test_scattering.py`

**Interfaces:**
- Consumes: `coupled_channel_vector`, `channel_block` (Task 2); `CoupledModel.hamiltonian`, `.interaction_matrix` (Task 1), `.channel_ells()`, `.charge`.
- Produces:
  - `CoupledSigma` — frozen dataclass with `E: NDArray[float64]` shape `(n_E,)`, `total: NDArray[float64]` and `restricted: NDArray[float64]` both shape `(n_E, len(vprimes))`.
  - `coupled_ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, entrance=0, ordering="COLAMD") -> CoupledSigma`.

**The two gates this task must pass.** At `s = 0` the coupled route must
reproduce `qscat.core.driven.ve_cross_section` — that is what certifies the
duplicated sweep. At `kappa = 0` a two-channel run must reproduce a
one-channel run, because only even Legendre components survive a symmetric
well so `l = 1` cannot reach `l = 2` at ANY anisotropy. Both are identities,
not tolerances.

- [ ] **Step 1: Write the failing test**

Append to `projects/no_coupled_channels/test_scattering.py`:

```python
import pytest
from qscat.core.driven import ve_cross_section

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel, DiagonalChannelModel
from projects.no_coupled_channels.scattering import coupled_ve_cross_section

E_TEST = np.array([0.02, 0.05])
VPRIMES = [0, 1]


def _basis(tg: TensorGrid):
    return vibrational_states(tg.grids[1], NO.mu, 3, NO.v0)


def test_s0_reproduces_the_shipped_solver() -> None:
    """The duplicated sweep must give the shipped answer where the models are
    the same Hamiltonian. This is what certifies the duplication."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)

    got = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    want = ve_cross_section(
        tg, DiagonalChannelModel(well=well, l=1), eps, chi, 0, VPRIMES, E_TEST
    )
    np.testing.assert_allclose(got.total, want, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(got.restricted, want, rtol=1e-12, atol=0.0)


def test_s0_many_channels_equal_one_channel() -> None:
    """The embedding identity, carried all the way to an observable.

    At s = 0 the coupled Hamiltonian is block-diagonal, so an electron
    entering channel 0 can never reach another channel: a four-channel run
    must give exactly the one-channel answer, and its total must equal its
    restricted part because no other exit is reachable. The preceding phase
    gated this at the Hamiltonian; this gates it through the solve, the
    projection and the normalisation as well.
    """
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.5)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    four = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=4), eps, chi, 0, VPRIMES, E_TEST
    )
    np.testing.assert_allclose(four.total, one.total, rtol=1e-10, atol=0.0)
    np.testing.assert_allclose(four.total, four.restricted, rtol=1e-10, atol=0.0)


def test_kappa_zero_two_channels_equal_one() -> None:
    """Parity: only even Legendre components survive a symmetric well, so
    l = 1 cannot reach l = 2 at ANY anisotropy. An identity, at s well away
    from zero where the coupling is otherwise fully on."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.0)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    two = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, E_TEST
    )
    np.testing.assert_allclose(two.total, one.total, rtol=1e-9, atol=0.0)
    np.testing.assert_allclose(two.restricted, one.restricted, rtol=1e-9, atol=0.0)


def test_the_coupling_moves_the_cross_section() -> None:
    """With kappa on, l = 1 reaches l = 2 and the answer must change --
    otherwise the previous test proves nothing."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    two = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, E_TEST
    )
    rel = np.abs(two.total - one.total) / np.maximum(one.total, 1e-30)
    assert float(np.max(rel)) > 1e-3


def test_total_is_at_least_the_restricted_part() -> None:
    """The total sums over exit channels and the restricted one is a single
    term of that sum, so the total can never be smaller."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    out = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=3), eps, chi, 0, VPRIMES, E_TEST
    )
    assert np.all(out.total >= out.restricted - 1e-15)


def test_a_closed_channel_contributes_nothing() -> None:
    """Below its threshold an exit channel is closed and must be zero, not a
    small number from an imaginary momentum."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.3, kappa=0.5)
    tiny = np.array([1e-4])  # far below the v' = 1 threshold
    out = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, tiny
    )
    assert out.total[0, 1] == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest projects/no_coupled_channels/test_scattering.py -v`
Expected: FAIL — `ImportError: cannot import name 'coupled_ve_cross_section'`

- [ ] **Step 3: Write the implementation**

Append to `projects/no_coupled_channels/scattering.py`, and extend `__all__`
to `["CoupledSigma", "channel_block", "coupled_channel_vector", "coupled_ve_cross_section"]`:

```python
@dataclass(frozen=True)
class CoupledSigma:
    """`sigma_{v_init->v'}(E)` in bohr^2, two ways.

    `total` sums over EXIT partial waves -- what an angle-integrated
    measurement sees, and the like-for-like partner of a fixed-l model's single
    exit. `restricted` keeps only the exit channel equal to the entrance, which
    isolates how the coupling changes the entrance amplitude through virtual
    excursions into other waves. Their difference is the flux the coupling
    redistributes, and it costs nothing extra: both come from one solve.
    """

    E: npt.NDArray[np.float64]
    total: npt.NDArray[np.float64]
    restricted: npt.NDArray[np.float64]


def coupled_ve_cross_section(
    tgrid: TensorGrid,
    model: CoupledModel,
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    vprimes: list[int],
    E: npt.ArrayLike,
    *,
    entrance: int = 0,
    ordering: Ordering = "COLAMD",
) -> CoupledSigma:
    """Exact 2-D coupled-channel `sigma_{v_init->v'}(E)`.

    The analysis is done ONCE and reused: `SparseLU` is built at the first
    energy that needs a solve and `refactor`ed per subsequent energy, since
    `E_tot * I - H` keeps one sparsity pattern across the sweep. On the
    SuperLU backend that reuse does not happen (it re-runs `splu`), which is
    why the production sweep needs MUMPS.

    Energies at or below threshold return zeros without any factorisation.
    """
    E_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    ells = model.channel_ells()
    n = tgrid.size
    charge = model.charge

    H = model.hamiltonian(tgrid)
    V = model.interaction_matrix(tgrid)
    ident = sp.identity(H.shape[0], dtype=np.complex128, format="csr")

    total = np.zeros((E_arr.size, len(vprimes)), dtype=np.float64)
    restricted = np.zeros((E_arr.size, len(vprimes)), dtype=np.float64)
    lu: SparseLU | None = None

    for i, e in enumerate(E_arr):
        if e <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        a = sp.csc_matrix(e_tot * ident - H)
        if lu is None:
            lu = SparseLU(a, ordering=ordering)
        else:
            lu.refactor(a)

        k = float(np.sqrt(2.0 * e))
        psi_i = coupled_channel_vector(tgrid, k, chi[v_init], ells, entrance, charge=charge)
        psi_plus = psi_i + lu.solve(V @ psi_i)
        v_psi = V @ psi_plus

        for j, vp in enumerate(vprimes):
            excess = e_tot - eps[vp]
            if excess <= 0.0:
                continue  # closed channel
            kp = float(np.sqrt(2.0 * excess))
            for c in range(len(ells)):
                # Only block `c` of the projection is non-zero, so take that
                # block rather than building a full-length vector of zeros.
                phi = channel_vector(tgrid, kp, chi[vp], ells[c], charge=charge)
                t = c_product(phi, channel_block(v_psi, c, n))
                s = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e))
                total[i, j] += s
                if c == entrance:
                    restricted[i, j] = s

    return CoupledSigma(E=E_arr, total=total, restricted=restricted)
```

Add the imports this needs at the top of the module: `from dataclasses import
dataclass`, `import scipy.sparse as sp`, `from qscat.linalg import Ordering,
SparseLU, c_product`, and `from projects.no_coupled_channels.model import
CoupledModel`.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest projects/no_coupled_channels/test_scattering.py -v`
Expected: PASS (9 tests)

**If `test_s0_reproduces_the_shipped_solver` fails, do NOT loosen it.** The two
routes solve the same equation with the same operator; a mismatch means the
interaction matrix, the block layout or the normalisation is wrong. Report the
numbers.

- [ ] **Step 5: Lint and commit**

```bash
uv run ruff check projects/no_coupled_channels
uv run ruff format projects/no_coupled_channels
git add projects/no_coupled_channels/scattering.py projects/no_coupled_channels/test_scattering.py
git commit -m "feat(coupled): the coupled-channel driven solve

sigma_VE(E) with the partial waves mixed: the Lippmann-Schwinger step becomes
a sparse matrix product because the interaction is no longer diagonal, the
entrance sits in one block, and the exit sums over blocks because the coupling
lets the electron leave in a wave it did not enter on. Reports both the total
and the entrance-restricted cross section from one solve.

Gated where it must be: at s = 0 it reproduces qscat.core.driven's shipped
solver, which is what certifies the duplicated sweep, and at kappa = 0 two
channels reproduce one by parity at full anisotropy.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

### Task 4: The threshold-aware energy mesh

**Files:**
- Create: `validation/coupled/energies.py`
- Test: `validation/coupled/test_energies.py`

**Interfaces:**
- Consumes: `qscat.core.vibrational.vibrational_states`, `validation.diatomic.config.CONFIGS`.
- Produces:
  - `E_LO = 0.002`, `E_HI = 0.150`, `BACKGROUND = 2.5e-4`, `CLUSTER_HALF = 5.0e-4`, `CLUSTER_STEP = 5.0e-5`, `DEDUP_TOL = 1e-6` (all Hartree).
  - `vibrational_thresholds(n_levels: int = 24) -> NDArray[float64]` — `eps_v - eps_0` for every level, on NO's nuclear deck.
  - `sweep_energies() -> NDArray[float64]` — the mesh, ascending, deduplicated.

- [ ] **Step 1: Write the failing test**

```python
# validation/coupled/test_energies.py
"""The mesh has to resolve two scales an order of magnitude apart: cusps AT
each threshold, and interference spread across the range on the scale of a
resonance width."""

from __future__ import annotations

import numpy as np

from validation.coupled.energies import (
    BACKGROUND,
    CLUSTER_STEP,
    DEDUP_TOL,
    E_HI,
    E_LO,
    sweep_energies,
    vibrational_thresholds,
)


def test_the_mesh_spans_the_declared_window() -> None:
    E = sweep_energies()
    assert E[0] >= E_LO
    assert E[-1] <= E_HI
    assert np.all(np.diff(E) > 0.0), "the mesh must be strictly ascending"


def test_no_two_energies_are_wastefully_close() -> None:
    """Every solve costs ~15 s on the production deck. Two energies a
    ten-thousandth of a mHa apart are a wasted one -- dedup must use a
    tolerance, not exact equality."""
    E = sweep_energies()
    assert float(np.min(np.diff(E))) >= DEDUP_TOL


def test_every_threshold_in_range_is_bracketed_closely() -> None:
    """A cusp is non-analytic AT the threshold, so points must sit on both
    sides of it, close."""
    E = sweep_energies()
    inside = [t for t in vibrational_thresholds() if E_LO < t < E_HI]
    assert len(inside) >= 15, f"expected many thresholds in range, got {len(inside)}"
    for t in inside:
        below = E[E <= t]
        above = E[E >= t]
        assert below.size and above.size, f"threshold {t} not bracketed"
        assert t - below[-1] <= CLUSTER_STEP * 1.01
        assert above[0] - t <= CLUSTER_STEP * 1.01


def test_the_background_is_no_coarser_than_declared() -> None:
    E = sweep_energies()
    assert float(np.max(np.diff(E))) <= BACKGROUND * 1.01


def test_the_mesh_is_the_size_the_spec_costed() -> None:
    """The nine-hour estimate is per-energy times this number; a mesh that
    silently doubled would silently double the run."""
    assert 900 <= sweep_energies().size <= 1150
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest validation/coupled/test_energies.py -v`
Expected: FAIL — `ModuleNotFoundError: ...energies`

- [ ] **Step 3: Write the implementation**

```python
# validation/coupled/energies.py
"""The energy mesh for the coupled VE sweep.

Two scales have to be resolved and they differ by an order of magnitude.
Threshold cusps are non-analytic AT each channel opening, so they need points
bracketing the threshold tightly rather than a fine mesh everywhere. The
overlapping-resonance interference is spread across the whole range instead,
on the scale of a resonance width -- the vibrational levels are spaced 7-8.5
mHa while the widths run 8-135 mHa, so the resonances genuinely overlap.

Hence a background grid plus clusters. Measured, this gives 1008 energies
against 2961 for a uniform mesh at the same 0.05 mHa resolution, which would
cost three times as much to resolve twenty places.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.vibrational import vibrational_states
from qscat.model import NO

from validation.diatomic.config import CONFIGS

__all__ = [
    "BACKGROUND",
    "CLUSTER_HALF",
    "CLUSTER_STEP",
    "DEDUP_TOL",
    "E_HI",
    "E_LO",
    "sweep_energies",
    "vibrational_thresholds",
]

# Hartree. The window covers the near-threshold region, both models' full
# resonance range (coupled 0.020-0.111, fixed-l 0.040-0.150) and twenty
# vibrational thresholds. NO's DA threshold is +0.172 Ha, above all of it.
E_LO = 0.002
E_HI = 0.150
# 0.25 mHa background: 30+ points across the narrowest width (8.2 mHa).
BACKGROUND = 2.5e-4
# 21 points at 0.05 mHa spanning +-0.5 mHa around each threshold.
CLUSTER_HALF = 5.0e-4
CLUSTER_STEP = 5.0e-5
# Two energies closer than this are one energy. Dedup by TOLERANCE, not exact
# equality: rounding alone leaves pairs a ten-thousandth of a mHa apart, and
# each one is a wasted ~15 s solve on the production deck.
DEDUP_TOL = 1e-6


def vibrational_thresholds(n_levels: int = 24) -> npt.NDArray[np.float64]:
    """`eps_v - eps_0` for each neutral level: where channel `v'` opens for VE
    out of `v = 0`."""
    nuclear = CONFIGS["NO"].da_grid().grids[1]
    eps, _chi = vibrational_states(nuclear, NO.mu, n_levels, NO.v0)
    return np.asarray(eps - eps[0], dtype=np.float64)


def sweep_energies() -> npt.NDArray[np.float64]:
    """The mesh: background grid plus a dense cluster at every threshold."""
    parts = [np.arange(E_LO, E_HI + 0.5 * BACKGROUND, BACKGROUND)]
    for t in vibrational_thresholds():
        if E_LO < t < E_HI:
            parts.append(
                np.arange(
                    t - CLUSTER_HALF, t + CLUSTER_HALF + 0.5 * CLUSTER_STEP, CLUSTER_STEP
                )
            )
    E = np.sort(np.concatenate(parts))
    E = E[(E >= E_LO) & (E <= E_HI)]
    keep = np.concatenate(([True], np.diff(E) >= DEDUP_TOL))
    return np.asarray(E[keep], dtype=np.float64)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest validation/coupled/test_energies.py -v`
Expected: PASS (5 tests)

- [ ] **Step 5: Report the mesh size**

Run: `uv run python -c "from validation.coupled.energies import sweep_energies as s; E=s(); print(E.size, E[0], E[-1])"`
Record the number in your report — it multiplies the run's cost directly.

- [ ] **Step 6: Lint and commit**

```bash
uv run ruff check validation/coupled
uv run ruff format validation/coupled
git add validation/coupled/energies.py validation/coupled/test_energies.py
git commit -m "feat(coupled): the threshold-aware energy mesh

Cusps are non-analytic AT a channel opening while the overlapping-resonance
interference is spread across the range on the scale of a width, so the mesh
is a 0.25 mHa background plus 21 points at 0.05 mHa around each of the twenty
vibrational thresholds. 1008 energies against 2961 for a uniform mesh at the
same resolution.

Dedup uses a tolerance rather than exact equality: rounding alone leaves pairs
a ten-thousandth of a mHa apart, and each is a wasted 15 s solve.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

### Task 5: The campaign driver and its timing probe

**Files:**
- Create: `validation/coupled/cross_section.py`
- Test: `validation/coupled/test_cross_section.py`

**Interfaces:**
- Consumes: `coupled_ve_cross_section`, `CoupledSigma` (Task 3); `sweep_energies` (Task 4); `TwoCentreWell`, `CoupledModel`; `CONFIGS["NO"].da_grid()`.
- Produces:
  - `S_RUN = 0.3`, `KAPPA_RUN = 0.5`, `V_INIT = 0`, `VPRIMES = [0, 1, 2, 3, 4]`, `N_CHANNEL_VALUES = (1, 3, 4)`, `N_VIB = 8`, `RESULTS = Path("validation/coupled/results")`.
  - `probe_one_energy(n_channels: int) -> dict[str, float]` — times a single energy on the production deck and returns `{"unknowns", "nnz", "build_s", "factor_s", "solve_s"}`.
  - `run_one(n_channels: int, energies: NDArray[float64]) -> CoupledSigma`.
  - `main(results: Path = RESULTS, energies: NDArray[float64] | None = None, grid: TensorGrid | None = None) -> dict[str, object]` — runs every `n_channels` and writes `cross_section.json`. `energies` and `grid` exist so a fast test can drive it on two energies and a tiny grid; both default to the production mesh and deck.

**The probe exists because the cost estimate is an extrapolation.** The spec's
nine hours comes from SuperLU scaling times a published MUMPS ratio. Task 6
runs `probe_one_energy` first and reports the real figure before committing the
machine to a sweep.

- [ ] **Step 1: Write the failing test**

```python
# validation/coupled/test_cross_section.py
"""Campaign structure, on a grid small enough to run in the fast tier."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid
from qscat.dvr import TensorGrid

from validation.coupled.cross_section import (
    KAPPA_RUN,
    N_CHANNEL_VALUES,
    S_RUN,
    VPRIMES,
    main,
)


def _tiny_grid() -> TensorGrid:
    el = segmented_grid(((4, 8.0),), ((2, 12.0),), angle_deg=35.0, quadrature=6)
    nu = segmented_grid(((3, 4.0),), ((2, 6.0),), angle_deg=30.0, quadrature=6, x_min=1.0)
    return TensorGrid([el, nu])


def test_the_declared_parameters_are_the_ones_measured() -> None:
    """S_RUN is where the preceding phase's oracle is converged to 0.2 %, and
    N_l = 2 is deliberately absent: it was measured 21 % from converged."""
    assert (S_RUN, KAPPA_RUN) == (0.3, 0.5)
    assert 2 not in N_CHANNEL_VALUES
    assert set(N_CHANNEL_VALUES) == {1, 3, 4}
    assert VPRIMES == [0, 1, 2, 3, 4]


def test_main_writes_every_model_on_one_mesh(tmp_path) -> None:
    """Two energies and a tiny grid: this checks the report's shape and that
    every model saw the SAME mesh, not the physics."""
    report = main(results=tmp_path, energies=np.array([0.02, 0.05]), grid=_tiny_grid())
    assert set(report["n_channels"]) == set(N_CHANNEL_VALUES)
    meshes = {tuple(report["sigma"][str(n)]["E"]) for n in N_CHANNEL_VALUES}
    assert len(meshes) == 1, "the branches must share one energy mesh"
    for n in N_CHANNEL_VALUES:
        s = report["sigma"][str(n)]
        assert np.shape(s["total"]) == (2, len(VPRIMES))
        assert np.shape(s["restricted"]) == (2, len(VPRIMES))
    assert (tmp_path / "cross_section.json").exists()


@pytest.mark.slow
def test_probe_reports_the_real_cost() -> None:
    """The nine-hour estimate is an extrapolation; this is the measurement."""
    from validation.coupled.cross_section import probe_one_energy

    out = probe_one_energy(1)
    assert out["unknowns"] == 78804
    assert out["factor_s"] > 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest validation/coupled/test_cross_section.py -v -m "not slow"`
Expected: FAIL — `ModuleNotFoundError: ...cross_section`

- [ ] **Step 3: Write the implementation**

```python
# validation/coupled/cross_section.py
"""The coupled VE cross section: does the fixed-l reduction change the observable?

Runs the exact 2-D driven solve at (s, kappa) = (0.3, 0.5) -- the anisotropy
where the preceding phase's channel truncation is converged to 0.2 % and the
width difference is 58 % over all 41 comparable R -- for N_l = 1 (fixed-l), 3
and 4, on one deck and one energy mesh.

N_l = 2 is deliberately absent from the sweep: it was measured 21 % from
converged, so it would be neither the oracle nor a useful convergence step. It
appears only in the parity identity gate, which lives in
`projects/no_coupled_channels/test_scattering.py`.

The production run needs MUMPS and does not belong on a laptop: measured,
SuperLU factorises the N_l = 2 deck in 208 s and gives no analysis reuse. Run
`probe_one_energy` first -- the cost estimate behind this campaign is an
extrapolation, and one measurement is cheaper than a wrong night.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel
from projects.no_coupled_channels.scattering import CoupledSigma, coupled_ve_cross_section
from validation.coupled.energies import sweep_energies
from validation.diatomic.config import CONFIGS

__all__ = [
    "KAPPA_RUN",
    "N_CHANNEL_VALUES",
    "N_VIB",
    "RESULTS",
    "S_RUN",
    "VPRIMES",
    "V_INIT",
    "main",
    "probe_one_energy",
    "run_one",
]

S_RUN = 0.3
KAPPA_RUN = 0.5
V_INIT = 0
VPRIMES = [0, 1, 2, 3, 4]
N_CHANNEL_VALUES = (1, 3, 4)
N_VIB = 8
RESULTS = Path("validation/coupled/results")


def _model(n_channels: int) -> CoupledModel:
    return CoupledModel(
        well=TwoCentreWell(base=NO, s=S_RUN, kappa=KAPPA_RUN), n_channels=n_channels
    )


def probe_one_energy(n_channels: int) -> dict[str, float]:
    """Measure one energy on the production deck: build, factor, solve.

    The campaign's cost estimate is an extrapolation from SuperLU scaling and
    a published MUMPS ratio. This is the measurement that replaces it, and it
    costs one energy rather than a night.
    """
    tgrid = CONFIGS["NO"].da_grid()
    model = _model(n_channels)

    t0 = time.perf_counter()
    H = model.hamiltonian(tgrid)
    build_s = time.perf_counter() - t0

    ident = sp.identity(H.shape[0], dtype=np.complex128, format="csr")
    a = sp.csc_matrix(0.05 * ident - H)
    t0 = time.perf_counter()
    lu = SparseLU(a)
    factor_s = time.perf_counter() - t0

    b = np.zeros(H.shape[0], dtype=np.complex128)
    b[0] = 1.0
    t0 = time.perf_counter()
    lu.solve(b)
    solve_s = time.perf_counter() - t0

    return {
        "n_channels": float(n_channels),
        "unknowns": float(H.shape[0]),
        "nnz": float(H.nnz),
        "build_s": build_s,
        "factor_s": factor_s,
        "solve_s": solve_s,
    }


def run_one(
    n_channels: int,
    energies: npt.NDArray[np.float64],
    grid: TensorGrid | None = None,
) -> CoupledSigma:
    """One model's sweep, on a deck the caller supplies.

    `main` builds the deck ONCE and threads it into every model, so that "one
    deck" is structural rather than three deterministic reconstructions that
    happen to agree. The `None` default exists only so a caller running a
    single model on its own need not build one.
    """
    tgrid = CONFIGS["NO"].da_grid() if grid is None else grid
    eps, chi = vibrational_states(tgrid.grids[1], NO.mu, N_VIB, NO.v0)
    return coupled_ve_cross_section(
        tgrid, _model(n_channels), eps, chi, V_INIT, VPRIMES, energies
    )


def main(
    results: Path = RESULTS,
    energies: npt.NDArray[np.float64] | None = None,
    grid: TensorGrid | None = None,
) -> dict[str, object]:
    """Every model on ONE mesh and ONE deck; writes `cross_section.json`."""
    E = sweep_energies() if energies is None else np.asarray(energies, dtype=np.float64)
    # ONE deck, built here and threaded into every model. Rebuilding it inside
    # the loop would give three deterministic reconstructions that agree by
    # luck of purity rather than by construction -- and on the production path
    # nothing would ever check that they had.
    tgrid = CONFIGS["NO"].da_grid() if grid is None else grid
    report: dict[str, object] = {
        "s": S_RUN,
        "kappa": KAPPA_RUN,
        "v_init": V_INIT,
        "vprimes": VPRIMES,
        "n_channels": list(N_CHANNEL_VALUES),
        "sigma": {},
    }
    for n_ch in N_CHANNEL_VALUES:
        t0 = time.perf_counter()
        out = run_one(n_ch, E, grid=tgrid)
        elapsed = time.perf_counter() - t0
        report["sigma"][str(n_ch)] = {  # type: ignore[index]
            "E": out.E.tolist(),
            "total": out.total.tolist(),
            "restricted": out.restricted.tolist(),
            "wall_clock_s": elapsed,
        }
        print(f"[coupled] N_l={n_ch}: {E.size} energies in {elapsed / 60:.1f} min")

    results.mkdir(parents=True, exist_ok=True)
    (results / "cross_section.json").write_text(json.dumps(report, indent=1))
    print(f"[coupled] wrote {results / 'cross_section.json'}")
    return report


if __name__ == "__main__":
    import sys

    if "--probe" in sys.argv:
        for n in N_CHANNEL_VALUES:
            print(probe_one_energy(n))
    else:
        main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest validation/coupled/test_cross_section.py -v -m "not slow"`
Expected: PASS (2 tests)

- [ ] **Step 5: Register the validation suite**

`.github/workflows/validation.yml` already has a `coupled` entry covering
`projects/no_coupled_channels validation/coupled`, so the new `@slow` test is
already inside a `validate:*` suite. Confirm rather than assume:

Run: `uv run pytest tests/test_validation_suites.py -v`
Expected: PASS

- [ ] **Step 6: Lint and commit**

```bash
uv run ruff check validation/coupled
uv run ruff format validation/coupled
git add validation/coupled/cross_section.py validation/coupled/test_cross_section.py
git commit -m "feat(coupled): the cross-section campaign and its timing probe

Runs N_l = 1, 3 and 4 at (s, kappa) = (0.3, 0.5) on one deck and one mesh.
N_l = 2 is deliberately absent -- measured 21 % from converged, so neither
oracle nor useful step.

probe_one_energy exists because the nine-hour estimate is an extrapolation
from SuperLU scaling and a published MUMPS ratio. One measured energy is
cheaper than a wrong night.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

### Task 6: Probe the real cost on sadaharu

**Files:** none — this task produces a measurement and a decision.

**Interfaces:**
- Consumes: `probe_one_energy` (Task 5).
- Produces: the measured per-energy cost, and a go/no-go on the sweep as specified.

- [ ] **Step 1: Build the MUMPS image on sadaharu**

The production deck needs the `test-deps` target, which carries MUMPS but does
not run the suite. From a checkout of this branch on sadaharu:

```bash
docker/build.sh test-deps
```

- [ ] **Step 2: Probe every channel count**

```bash
docker run --rm -v "$PWD":/work -w /work qmodeling:test-deps \
  uv run --no-sync python -m validation.coupled.cross_section --probe
```

Expected shape: three dicts, for `n_channels` 1, 3 and 4, each with
`unknowns`, `nnz`, `build_s`, `factor_s`, `solve_s`. At `n_channels=4`,
`unknowns` must be 315216.

- [ ] **Step 3: Compare against the estimate and decide**

The spec assumed ~15 s per energy at `N_l = 4`, giving 4.2 h for that model and
about nine hours for all three. The per-energy cost is roughly
`factor_s + len(VPRIMES) * n_channels * solve_s` — the factorisation dominates.

Multiply the measured figure by the mesh size (Task 4 Step 5) and compare.

- **Within about 2× of the estimate:** proceed to Task 7 as specified.
- **More than 2× over:** stop and report. Do not silently thin the mesh or drop
  a channel count — the mesh resolves cusps that are the point of the run, and
  `N_l = 3` is the convergence check without which `N_l = 4` cannot be quoted.
  The decision about what to trade is not the implementer's.

Also check `OMP_NUM_THREADS`. sadaharu has 32 cores and MUMPS threads its
factorisation, so leaving it unset is probably right here — but the repository
has a recorded case where a 32-thread OpenBLAS was ~400× SLOWER on tiny
eigenproblems, so measure with it set to 1 and unset, and report both.

- [ ] **Step 4: Record**

Put the measured numbers in your report: per-energy cost per channel count,
the projected total, whether threading helped, and the go/no-go.

---

### Task 7: Run the campaign

**Files:**
- Create: `validation/coupled/results/cross_section.json` (generated)

- [ ] **Step 1: Run it**

```bash
docker run --rm -v "$PWD":/work -w /work qmodeling:test-deps \
  uv run --no-sync python -m validation.coupled.cross_section
```

Expect one progress line per channel count. If it exceeds twice the projection
from Task 6, stop it and report rather than letting it run.

- [ ] **Step 2: Sanity-check the output before trusting it**

```bash
uv run python -c "
import json, numpy as np
d = json.load(open('validation/coupled/results/cross_section.json'))
for n in d['n_channels']:
    s = d['sigma'][str(n)]
    tot = np.asarray(s['total'])
    print(n, 'shape', tot.shape, 'finite', np.isfinite(tot).all(),
          'elastic max', tot[:, 0].max(), 'in', s['wall_clock_s']/60, 'min')
"
```

Every entry must be finite and non-negative, and the elastic channel must be
non-zero. A cross section of exactly zero everywhere means the entrance never
coupled to the exit — a bug, not a result.

- [ ] **Step 3: Commit the artifact**

```bash
git add validation/coupled/results/cross_section.json
git commit -m "data(coupled): the coupled VE cross-section campaign

N_l = 1, 3 and 4 at (s, kappa) = (0.3, 0.5), <N> energies on NO's production
deck, MUMPS on sadaharu, <T> hours.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

Replace `<N>` and `<T>` with the measured mesh size and wall-clock.

- [ ] **Step 4: Report the result**

In your report, give for each channel count: the elastic and $0\to1$ cross
sections at a few energies, and — the point of the whole campaign — **a measure
of how much structure each curve has.** The prediction under test is that the
fixed-$l$ curve is substantially SMOOTHER. A serviceable measure is the total
variation $\sum_i |\sigma_{i+1} - \sigma_i|$ divided by the mean, computed per
channel; report it for $N_l = 1$ against $N_l = 4$, and say whether the
prediction held.

Report the convergence gate in the same breath: the median relative difference
between $N_l = 3$ and $N_l = 4$ over the sweep, per channel. The $N_l = 4$
numbers cannot be quoted as an oracle without it, and the preceding phase
measured $N_l = 2$ as inadequate precisely by this test.

---

### Task 8: The figure, the note and the repo map

**Files:**
- Modify: `validation/coupled/figures.py`
- Modify: `docs/physics/coupled-partial-waves.md`
- Modify: `CLAUDE.md`
- Create: `docs/physics/figures/no-coupled-cross-section.png` (generated)

- [ ] **Step 1: Add the cross-section figure**

Append to `validation/coupled/figures.py`:

```python
CROSS_SECTION_FIGURE = "docs/physics/figures/no-coupled-cross-section.png"


def cross_section_figure() -> str:
    """Fixed-l against coupled, one panel per channel, thresholds marked.

    The thin vertical lines are the vibrational thresholds: with them a reader
    can tell a cusp (pinned to a line) from a resonance (not), which is the
    whole reason the mesh clusters where it does.
    """
    import json

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from validation.coupled.energies import E_HI, E_LO, vibrational_thresholds

    d = json.loads((RESULTS / "cross_section.json").read_text())
    E = np.asarray(d["sigma"]["1"]["E"], dtype=float)
    thresholds = [t for t in vibrational_thresholds() if E_LO < t < E_HI]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for vp, ax in enumerate(axes.ravel()):
        for n_ch, style, colour in (("1", "--", "tab:orange"), ("4", "-", "tab:blue")):
            sig = np.asarray(d["sigma"][n_ch]["total"], dtype=float)[:, vp]
            ax.plot(E, sig, style, color=colour, lw=1.0, label=f"$N_l$ = {n_ch}")
        for t in thresholds:
            ax.axvline(t, color="0.85", lw=0.5, zorder=0)
        ax.set(ylabel="$\\sigma$ ($a_0^2$)", title=f"$0 \\to {vp}$")
        ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("$E$ (Ha)")
    fig.suptitle(
        "NO vibrational excitation: coupled against fixed-$l$, "
        f"$s$ = {d['s']}, $\\kappa$ = {d['kappa']}",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(CROSS_SECTION_FIGURE, dpi=130)
    plt.close(fig)
    print(f"[coupled] wrote {CROSS_SECTION_FIGURE}")
    return CROSS_SECTION_FIGURE
```

and call it from the module's `main()` alongside the existing figure.

Run: `uv run python -m validation.coupled.figures` and **look at the result** —
a matplotlib test skips in CI, so a green CI run is not evidence it rendered.

- [ ] **Step 2: Extend the note**

Add a section to `docs/physics/coupled-partial-waves.md` — the cross section is
the same story's next chapter and belongs in the same note, not a new one.

It must carry: the measured structure comparison and whether the smoothing
prediction held; the numbers with their channel counts; the $N_l = 3$ against
$N_l = 4$ convergence; the run cost; and the limits — one anisotropy, one
entrance channel, nothing compared with experiment, the anisotropy geometric
rather than fitted, and the forty duplicated lines of sweep with their
promotion path. Embed the figure inline.

Constraints the repo's tests enforce: **no Greek inside backticks** (use
`$...$`), **no MyST directives** in `docs/physics/`, and **no citation of the
spec, the plan, a PR or a task number** — the note must stand alone for a
reader with only the clone.

- [ ] **Step 3: Update the repo map**

`CLAUDE.md`'s `validation/coupled/` entry gains the cross-section route:
`energies.py`, `cross_section.py`, the run recipe, and the headline result.

- [ ] **Step 4: Full verification**

```bash
uv run pytest projects/no_coupled_channels validation/coupled -v -m "not slow"
uv run pytest tests/ -q
uv run ruff check .
uv run ruff format validation projects
uv run sphinx-build -b html docs docs/_build/html -W --keep-going
```

- [ ] **Step 5: Commit**

```bash
git add validation/coupled/figures.py docs/physics/coupled-partial-waves.md docs/physics/figures/no-coupled-cross-section.png CLAUDE.md
git commit -m "docs(coupled): the cross-section result

<one line: did the smoothing prediction hold, and by how much>

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01V3EYxZUKnE39YBX3mRWVjt"
```

---

## Self-review notes (for the executor, not steps to run)

Three things this plan deliberately does, so nobody undoes them:

1. **The interaction matrix carries the perturbation alone** — no `v0`, no
   centrifugal term. Those belong to the free channel Hamiltonian that
   `channel_vector` solves. Adding them would drive Lippmann-Schwinger with
   the wrong operator and produce a plausible wrong T-matrix with no error.
2. **`N_l = 2` is absent from the sweep on purpose.** It is 21 % from
   converged. It appears only in the parity identity gate, where being a
   *different* channel count is the whole point.
3. **The duplicated sweep is gated, not apologised for.** Task 3's
   `test_s0_reproduces_the_shipped_solver` is what makes forty duplicated
   lines acceptable at toy stage; if that test is ever weakened, the
   duplication stops being defensible.

One risk worth naming: `test_kappa_zero_two_channels_equal_one` asserts an
identity at `rtol=1e-9` through a full driven solve rather than a matrix
comparison, so it rides on the linear solve being deterministic. If it fails at
the tenth digit while `test_s0_reproduces_the_shipped_solver` passes, that is a
solver-conditioning observation, not a coupling bug — report the numbers rather
than loosening either.
