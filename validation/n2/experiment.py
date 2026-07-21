"""N₂ LCP benchmark harness. Prints a PASS/PENDING/FAIL table; exits non-zero on FAIL.

Run: python validation/n2/experiment.py  (or in Docker: docker run --rm qmodeling:runtime
python validation/n2/experiment.py)
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import loader
import model
import resonance
from qscat.units import HARTREE_TO_EV

import reference


def _load_local(module_name: str, internal_name: str):  # type: ignore[no-untyped-def]
    """Load a sibling module by file path under a private internal name.

    `projects/n2_ti_cross_section` also has a top-level module literally
    named `cross_section.py`. When both that project's tests and this
    harness are collected in the same pytest session, whichever one a bare
    `import cross_section` sees FIRST wins the ambient `sys.modules
    ['cross_section']` cache -- a plain `import cross_section` here would
    then silently pick up the wrong module depending on test collection
    order. Loading THIS file explicitly by path under a distinct name
    sidesteps that entirely (mirrors the same technique `cross_section.py`
    itself uses to load the *other* `cross_section.py`, in reverse).
    """
    path = Path(__file__).resolve().parent / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(internal_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[internal_name] = module  # dataclasses/pickling need this registered
    spec.loader.exec_module(module)
    return module


cross_section = _load_local("cross_section", "_n2_harness_cross_section")

Check = tuple[str, str, str, str]  # (group, name, status, detail)
# Status values: PASS, FAIL, PENDING (no result yet), or NOTE (a result exists
# but is a documented, non-gating observation -- never counted as a FAIL).


def run_checks() -> list[Check]:
    checks: list[Check] = []

    # Group A — closed-form model (green now)
    for name, ok, detail in model.model_checks():
        checks.append(("A model", name, "PASS" if ok else "FAIL", detail))

    # Group C1–C4 — golden-data integrity (green now)
    for name, ok, detail in loader.integrity_checks():
        checks.append(("C data", name, "PASS" if ok else "FAIL", detail))

    # Group B — resonance position via two-angle ECS matching (green now)
    lo, hi = reference.LITERATURE["E_res_eV"]
    residual_tol = 1e-3  # Ha; angle-stability threshold for a genuine pole
    try:
        E_pole, residual = resonance.e_res_at_R0()
    except Exception as e:
        checks.append(("B resonance", "B1 E_res(R0) in literature window", "FAIL",
                        f"pole computation failed: {e}"))
    else:
        E_res_eV = E_pole.real * HARTREE_TO_EV
        b1_ok = lo <= E_res_eV <= hi and residual < residual_tol
        detail = (f"E_res={E_res_eV:.3f} eV (expect {lo}-{hi} eV; "
                  f"match residual={residual:.2e} Ha, expect <{residual_tol:.0e} Ha)")
        if not b1_ok and lo <= E_res_eV <= hi:
            detail += " -- residual too large for an angle-stable pole"
        checks.append(("B resonance", "B1 E_res(R0) in literature window",
                       "PASS" if b1_ok else "FAIL", detail))

    # Group C5 — cross-section value anchors vs Houfek data, via the TI solver
    # (projects/n2_ti_cross_section). GATED anchors (a VE channel clear of its own
    # threshold) get a real PASS/FAIL at reference.ANCHOR_FACTOR; DOCUMENTED-LIMITED
    # anchors (elastic, or too close to their own threshold) are a NOTE -- their ratio
    # is always printed, but they never fail the harness (see cross_section.py).
    try:
        anchor_results = cross_section.compute_anchor_results()
    except Exception as e:
        checks.append(("C anchors", "C5 sigma anchors (6, TI solver)", "FAIL",
                        f"TI solver failed: {e}"))
    else:
        for r in anchor_results:
            lbl = "elastic" if r.channel == 0 else f"v=0->{r.channel}"
            name = f"C5 sigma({r.energy_ha:.4g} Ha, {lbl})"
            core = (f"computed={r.sigma_computed:.4e} bohr^2, houfek={r.sigma_houfek:.4e} "
                    f"bohr^2, ratio={r.ratio:.3f}")
            if r.gated:
                ok = 1.0 / reference.ANCHOR_FACTOR <= r.ratio <= reference.ANCHOR_FACTOR
                detail = f"{core} (LCP vs Houfek 2D, factor={reference.ANCHOR_FACTOR:.1f})"
                checks.append(("C anchors", name, "PASS" if ok else "FAIL", detail))
            else:
                checks.append(("C anchors", name, "NOTE", f"{core} -- {r.mechanism}"))

    # Group D — time-dependent model: PENDING (later)
    checks.append(("D time-dependent", "D1 TD cross sections", "PENDING",
                   "needs time-dependent LCP propagation"))
    return checks


def main() -> int:
    checks = run_checks()
    width = max(len(f"{g}: {n}") for g, n, _s, _d in checks)
    print("N2 LCP benchmark harness")
    print("=" * (width + 30))
    for group, name, status, detail in checks:
        print(f"[{status:7}] {group}: {name:<{width - len(group) - 2}}  {detail}")
    n_pass = sum(c[2] == "PASS" for c in checks)
    n_pend = sum(c[2] == "PENDING" for c in checks)
    n_note = sum(c[2] == "NOTE" for c in checks)
    n_fail = sum(c[2] == "FAIL" for c in checks)
    print("=" * (width + 30))
    print(f"{n_pass} PASS, {n_pend} PENDING, {n_note} NOTE, {n_fail} FAIL")
    # NOTE rows are documented, non-gating observations (known LCP-model
    # limitations) -- they never count toward FAIL / the exit code.
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
