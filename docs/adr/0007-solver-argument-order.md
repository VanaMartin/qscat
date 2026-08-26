# 7. Solver argument order

Date: 2026-08-25

## Status

Accepted

## Context

The public solvers grew across eight sub-projects and the 2026-08-25 release
review asked whether their argument orders follow a rule or an accident.
Reading them all: they follow a rule with one deliberate exception, but the
rule was never written down, so every new solver re-derived it.

## Decision

1. **Observable solvers** — functions computing a physical observable on a
   prepared discretisation (`*_cross_section`, level/state solvers that take a
   vibrational basis) — take the discretisation first:

       f(grid(s), model, eps, chi, v_init, <per-call physics>, *, options)

   `<per-call physics>` is what varies between calls on one problem
   (`vprimes`, `E`, `phi_d`); everything tunable-but-stable is keyword-only.
   Examples: `ve_cross_section`, `da_cross_section`, `dr_cross_section`, the
   four `td_*` solvers, `nrm_ve_cross_section`, `nrm_da_cross_section`.

2. **Model-derived builders** — functions whose output is derived from the
   model itself (potential curves, resonance-level pipelines, proposed grids,
   incident placement) — take the model first, then the grids they need:
   `resonance_levels`, `local_complex_potential`, `exact_resonance_states`,
   `bo.resonance_curve`, `propose_grid`, `tw_analysis`.

3. **`lcp_da_cross_section` is a documented exception to (1)**: it takes
   `(nuclear_grid, mu, Vd, Gamma, eps, chi, v_init, E, ...)` — a bare reduced
   mass and a curve, not a model — because the LCP equation contains no model:
   its physics input IS the curve, which legitimately arrives from
   `resonance_levels(return_curve=True)`, from a fit, or from a file. A
   model-accepting signature would either hide an expensive electronic pole
   walk inside a cross-section call or carry a redundant argument. The
   `ScatteringProblem.lcp_da_cross_section` method supplies
   `mu`/`eps`/`chi`/`v_init` from its bundle, so the exception costs facade
   users nothing. The exception is scoped to the LCP family — curve-input
   solvers with no model to lead their argument list — not to this one
   function; `lcp_resonance_levels` is the same case (it diagonalizes a
   supplied `Vd`/`Gamma` curve, not a model) and follows the same order for
   the same reason.

4. New solvers follow (1) or (2); a new exception needs its own recorded
   reason, in its docstring and here.

## Consequences

- A reader can predict any solver's leading arguments from its kind, and the
  `ScatteringProblem` facade can bundle the shared `(grid, model, eps, chi,
  v_init)` group mechanically.
- The rule is also recorded in the `qscat-conventions` skill for lookup.
