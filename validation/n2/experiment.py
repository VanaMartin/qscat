"""N₂ LCP benchmark harness. Prints a PASS/PENDING/FAIL table; exits non-zero on FAIL.

Run: python validation/n2/experiment.py  (or in Docker: docker run --rm qmodeling:runtime
python validation/n2/experiment.py)
"""

from __future__ import annotations

import sys

import loader
import model
import resonance
from qscat.units import HARTREE_TO_EV

import reference

Check = tuple[str, str, str, str]  # (group, name, status, detail)


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
    E_pole, residual = resonance.e_res_at_R0()
    E_res_eV = E_pole.real * HARTREE_TO_EV
    b1_ok = lo <= E_res_eV <= hi
    checks.append(("B resonance", "B1 E_res(R0) in literature window",
                   "PASS" if b1_ok else "FAIL",
                   f"E_res={E_res_eV:.3f} eV (expect {lo}-{hi} eV; "
                   f"match residual={residual:.2e} Ha)"))

    # Group C5 — cross-section value anchors vs Houfek data (needs TI solver): PENDING
    for e, ch, ref in reference.anchors():
        lbl = "elastic" if ch == 0 else f"v=0->{ch}"
        checks.append(("C anchors", f"C5 sigma({e:.4g} Ha, {lbl})", "PENDING",
                       f"ref={ref:.4e} bohr^2, rtol={reference.RTOL:.0%}; needs TI solver"))

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
    n_fail = sum(c[2] == "FAIL" for c in checks)
    print("=" * (width + 30))
    print(f"{n_pass} PASS, {n_pend} PENDING, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
