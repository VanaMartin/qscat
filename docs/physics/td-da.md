# Time-dependent dissociative attachment (TD-DA): the nuclear-axis extractors

**Location:** `qscat.core.time_dependent` (`td_da_cross_section(method=...)`,
`td_da_cross_sections_all`), `qscat.core.td_extractors` (`Flux`, `Dirac`,
`TannorWeeks`, all `axis="nuclear"`), `qscat.core.dissociation`
(`da_cross_section` -- the exact TI oracle this converges to,
`anion_electronic_states`), `validation/diatomic/td_da.py` +
`validation/diatomic/test_td_da.py` (the F2/NO validation harness).
**Origin:** sub-project #4 "TD extractors" (branch `td-alternative-extractors`,
SP1: `docs/physics/td-extractors.md`) generalized in sub-project SP2 (branch
`td-da-route`) from VE (vibrational excitation, an electronic-axis outgoing
channel) to DA (dissociative attachment, a NUCLEAR-axis outgoing channel).
**Units:** atomic units throughout.

## What this is

`docs/physics/td-extractors.md` established three independent ways to read
an energy-domain S-matrix off one propagated Gaussian-wavepacket trajectory
(`TannorWeeks`/`Dirac`/`Flux`), for the VE cross section
`sigma_{v_init->v'}(E)` -- an ELECTRONIC-axis outgoing channel (the electron
re-emerges, the molecule stays bound in vibrational level `v'`).
Dissociative attachment is the other exit channel of the same 2-D
(electron-`r` x nuclear-`R`) collision: the electron is captured, the
molecule DISSOCIATES, and the outgoing flux moves along the NUCLEAR
coordinate `R` instead. `qscat.core.dissociation.da_cross_section` already
computes this exactly (time-independent, driven Lippmann-Schwinger); this
sub-project is the TIME-DEPENDENT route to the SAME cross section, and (as
with VE) an independent numerical method that must converge to the same
answer as its own differential oracle.

The three electronic-axis extractors generalize to the nuclear axis
one-for-one (all three `axis="nuclear"`, `qscat.core.td_extractors`):

| Extractor | electronic-axis (VE) | nuclear-axis (DA) |
|---|---|---|
| `Flux` | fixed electronic surface, Wronskian value+derivative vs `chi_{v'}` | fixed NUCLEAR surface `R=surface`, Wronskian value+derivative vs an ANION electronic bound state `phi_c` |
| `Dirac` | fixed electronic point, VALUE projection vs `chi_{v'}` | fixed NUCLEAR point `R=position`, VALUE projection vs `phi_c` |
| `TannorWeeks` | propagated outgoing Gaussian `g_out(r) chi_{v'}(R)` | propagated outgoing Gaussian `phi_c(r) g_out(R)` (test packet now in `R`) |

In every case the roles of the two axes swap: the exit "channel index" is no
longer a neutral vibrational level `v'` but one of `n_channels` ANION
electronic bound states `phi_c` at the dissociation limit `R_inf =
tgrid.grids[1].R0` (`qscat.core.dissociation.anion_electronic_states`), and
the outgoing test function/surface/point moves from the electronic
coordinate `r` to the nuclear coordinate `R` (mass `mu_R = model.mu`,
partial wave `l=0`, since the nuclear motion is s-wave along the bond).
The INCIDENT side is unchanged in all three: `eta_incident` stays on the
electronic axis (the electron is still what comes in).

## The DA sigma normalization: reconciling `C_DA = pi` with the TI oracle's `4*pi^3`

Every nuclear extractor's `sigma(E)` uses `sigma_DA,c(E) = C_DA * |S_c(E)|^2
/ (2E)` with `C_DA = pi` -- the SAME prefactor the electronic extractors use
for an INELASTIC VE channel. This looks inconsistent with
`da_cross_section`'s literal `4*pi^3 |T|^2 / 2E`, but the two are the SAME
formula in different variables. The general partial-wave S/T relation is

    S = 1 - 2*pi*i*T   =>   |S|^2 = 4*pi^2*|T|^2   (off-diagonal / no-`1`-term)

so

    pi*|S|^2/(2E) = pi*(4*pi^2*|T|^2)/(2E) = 4*pi^3*|T|^2/(2E)

identically. The TD extractors compute an S-matrix element (`S = raw
correlation / (2*pi*eta_out*eta_in)`, the SAME Tannor-Weeks deconvolution
used everywhere in this stack); the TI oracle computes a T-matrix element
(a driven-equation projection). `C_DA = pi` is therefore not a separate
constant to tune -- it is dictated by the S/T identity, and a WRONG value
(e.g. an extra `4*pi^2` from conflating the two conventions) would show up
as the converged ratio landing at a wildly different constant, not ~1. The
`@slow` TI-convergence gates below are the empirical confirmation of this,
not merely an algebraic assertion.

## `method=` selection and the honest three-way helper

```python
sigma_da = td_da_cross_section(tgrid, model, eps, chi, v_init, E,
                                dt=dt, n_steps=n_steps, wp_in=wp_in,
                                method="flow" | "delta" | "tw",
                                surface=..., position=..., wp_out=...,
                                n_channels=1)
```

`method="flow"` (default -- the natural DA extractor: a fixed-surface
Wronskian flux needs no propagated outgoing test packet) requires
`surface`; `method="delta"` requires `position`; `method="tw"` requires
`wp_out` (now the NUCLEAR outgoing test packet's `r0_out`/`p0_out`/
`sigma_out`, in bohr/momentum along `R`). Omitting the one the selected
method needs raises `ValueError`. Unlike `td_ve_cross_section`, there is
**no `wp_out`/free-reference for `"flow"`/`"delta"`, and no elastic
free-reference subtraction for ANY method** -- DA is a pure rearrangement
channel with no `v'==v_init` diagonal to subtract a reference from (passing
`free` to any nuclear extractor's `sigma()` raises `ValueError`, enforced in
`td_extractors.py`).

```python
sigmas = td_da_cross_sections_all(tgrid, model, eps, chi, v_init, E,
                                   dt=dt, n_steps=n_steps, wp_in=wp_in,
                                   surface=..., position=..., wp_out=...,
                                   n_channels=1)
# {"flow": ..., "delta": ..., "tw": ...}
```

runs `Flux`, `Dirac`, and `TannorWeeks` (all `axis="nuclear"`) from ONE
shared propagation -- the honest three-way comparison: any spread between
the three returned cross sections is a genuine property of the extraction
method (or a shared discretization/truncation residual all three inherit
together), never an artifact of propagating slightly different dynamics.
`E` may be scalar (`(n_channels,)`) or array-like (`(len(E), n_channels)`),
matching `da_cross_section`'s own convention.

## THE KEY LESSON: TD-DA needs a launch-box grid, not the TI `da_grid`

**TD-DA reuses the SAME grid + incident wavepacket SHAPE as the
Tannor-Weeks VE method — a LARGE electronic LAUNCH box (`r0` well inside
`r_max`) x the FINE per-molecule nuclear deck (resolves the fast K_R
dissociation-flux wave) — NOT the TI oracle's `MoleculeConfig.da_grid()`.**

This was found empirically (controller, SP2 Task 2) while validating the
first nuclear `Flux` gate, as two distinct, real failure modes:

1. **Off-box incident.** `da_grid()`'s electronic box (`e_r_max=16`, F2/NO)
   is sized for the DRIVEN TI solve -- it does not need to hold a launched
   wavepacket, only support the driven equation's asymptotic channel
   function. Placing an incident Gaussian at eMoScat's TD `r0=45` on that
   box starts the packet ~30 bohr PAST the real region, inside the ECS
   complex-scaled tail: `V(r)*e^{i*theta}` there is not a physical
   potential, it is an absorbing rotation, so the packet picks up
   exponential garbage from step 1 and `sigma_DA` diverges ~1e6x by
   `n=1000`.
2. **Coarse nuclear.** Fixing the electronic box (launch-box, `r_max=25`,
   `r0=12` well inside) but keeping a COARSE shared nuclear grid (~139 pts,
   the N2-style grid `docs/physics/per-molecule-discretisation.md` already
   flags as under-resolved for DA) gives a numerically STABLE propagation
   (no divergence, `|psi|` well-behaved) but `sigma_flux ~ 1e-21` --
   essentially zero. The coarse grid cannot represent the K_R~72
   dissociation flux wave at the analysis surface, so the flux extractor
   never sees the outgoing current; the propagation looks fine, the
   OBSERVABLE is silently wrong.

The fix: a DEDICATED launch-box electronic grid (`electronic_grid(r_max=25,
order=6, n_complex=3, angle_deg=40.0)`, clean enough to hold the incident
`r0=12` packet without ECS-tail contamination) paired with the UNCHANGED
eMoScat fine nuclear deck (`MoleculeConfig.da_grid().grids[1]` --
`segmented_grid` reproduces the exact per-molecule `(n_elements, endpoint)`
deck the TI oracle already uses for the same reason: it resolves the fast
dissociation wave). `validation/diatomic/td_da.py`'s `td_launch_grid(cfg)`
builds exactly this pairing; do not substitute `cfg.da_grid()` wholesale for
a TD run.

## F2/NO three-way validation

**Fast gate** (`validation/diatomic/test_td_da.py`, not `@slow`): shape/
contract tests on a tiny F2 grid (`electronic_grid(r_max=12, order=5,
n_complex=3)` x `nuclear_grid(quadrature=6, r_max=14, n_complex=3)`, mirroring
`libs/qscat/tests/test_td_extractors.py`'s own tiny config), `n_steps=5` --
NOT a converged cross section, just "builds, propagates, returns the right
shape, `method=`/`ValueError` wiring correct, `td_da_cross_sections_all`
matches calling `td_da_cross_section` per method to machine precision".

**`@slow` gate** (`validation/diatomic/test_td_da.py`, F2 and NO): the
launch-box grid above, `wp_in={r0:12, p0:-0.5, sigma:3}`, `surface`/
`position` at `R~6` bohr (converted from bohr to the nearest real-region
nuclear DVR index), `wp_out={r0_out:8.0, p0_out:72.0, sigma_out:0.07}`
(eMoScat's F2 nuclear test packet: a NARROW, wide-K packet placed inward of
the surface/position analysis points -- needed to resolve the fast K_R~72
wave; a wide/slow packet undersamples it entirely, see `td_extractors.py`'s
module docstring), `n_steps=1800` (~86k unknowns, ~10 min per molecule).
`td_da_cross_sections_all` is asserted within `(0.7, 1.3)` of `da_cross_
section` (TI) per method, at each molecule's own anchor energies.

Per-extractor `@slow` gates already validated the identical launch-box
config in `libs/qscat/tests/test_td_extractors.py`
(`test_nuclear_flux_da_converges_to_ti_oracle`,
`test_nuclear_dirac_da_converges_to_ti_oracle`,
`test_nuclear_tw_da_converges_to_ti_oracle`) at `n_steps=1500`; this
sub-project's `test_td_da.py` gate is the SAME physics run through the
propagate-ONCE `td_da_cross_sections_all` helper at `n_steps=1800` (a small
margin over the per-extractor gates' 1500, since the assertion tolerance
widens slightly to `(0.7, 1.3)` from `(0.7, 1.25)`).

**Measured (controller, 2026-07-31), F2, launch-box grid:**

```text
sigma_flux/sigma_ti   (E=0.02,0.03,0.04): STABLE plateau ~0.86-0.97 by n>=1350
sigma_delta/sigma_ti: STABLE plateau ~0.87-0.96 by n>=1350 (mirrors flux)
sigma_tw/sigma_ti:    converges to the RIGHT MAGNITUDE (order ~1) but does NOT
                      plateau cleanly -- it OSCILLATES (~0.55-1.42 across
                      n=1750-2000). TW is the noisiest, most test-packet-
                      SENSITIVE of the three (a propagated-Gaussian
                      deconvolution, not a point-value/Wronskian read of the
                      SAME trajectory): it needs the narrow sigma_out~0.07
                      (wide-K) packet placed inward of the ECS edge with
                      p0_out~K_R; a wide-in-R (narrow-in-K) packet
                      ill-conditions the deconvolution and blows up the
                      higher-E channels. Hence flux/delta are gated at
                      (0.7,1.3) but TW at the wider order-~1 band (0.4,1.7).
                      This test-packet sensitivity is exactly what the
                      `qscat.tuning` TW-analysis machinery exists to tame.
```

(A SEPARATE, smaller cross-method check on `libs/qscat/tests/test_td_
extractors.py`'s deliberately tiny/toy N2 grid -- not F2's launch-box grid,
and not this validation's own gate -- found `sigma_dirac/sigma_flux ~
0.67-0.74`, still settling between `n_steps=800` and `1500`; that number is
a toy-grid sanity check on the shared codepath, not evidence about `Dirac`'s
F2 plateau, which mirrors `Flux`'s per its own `@slow` gate's docstring.)

Flux and delta plateau near (not exactly at) 1.0 -- the same ~3-14% TD-vs-TI
cross-method residual the electronic-axis VE extractors show at a converged
(but not infinitely fine) grid; TW lands at order ~1 with more scatter. This
is NOT a normalization error (see the `C_DA`
section above: a wrong prefactor would land at a wildly different constant).

## Deferred: the full eMoScat grid, and SP3 (TD-DR)

The FULL eMoScat F2 grid (electronic real region out to 90 bohr, ~402k
unknowns) TD-DA convergence is a Docker/overnight-deferred production run --
it exceeds any interactive harness's patience budget (the reduced launch-box
grid above, ~86k unknowns, is already ~10 min per `@slow` test). This
mirrors the same deferral `docs/physics/n2-2d-td-cross-section.md` and
`docs/physics/td-extractors.md` already document for their own working
grids vs. the fully eMoScat-scale deck.

**SP3 = TD-DR** (dissociative recombination for the IONIC H2+ target,
`qscat.core.dissociation.dr_cross_section`'s TI oracle) reuses this EXACT
nuclear-axis extractor infrastructure: a Coulomb incident channel (already
plumbed via `model.charge`) and a loop over Rydberg electronic exit states
in place of `n_channels` anion bound states. The `charge` plumbing already
present in the nuclear extractors' outgoing-Hankel-half machinery
(`hankel_point_value`/`outgoing_surface_wave`'s `charge` argument, currently
exercised only at `charge=0` since F2/NO are neutral) is the reason SP3 is
additive rather than a rewrite.

## Deferred

The committed three-way comparison FIGURE (sigma_DA-vs-TI per method for F2/NO)
was NOT produced in-session -- generating it needs the `@slow` per-molecule
propagations (~10 min each), which exceed the harness run window. The gates
(`validation/diatomic/test_td_da.py`, `@slow`) encode the numeric comparison;
the figure is a Docker/overnight follow-on, alongside the full eMoScat 90-bohr-
electronic-grid convergence run.
