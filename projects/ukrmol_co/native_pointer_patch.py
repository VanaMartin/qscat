"""Add CSF identity guards and Class1 progress markers to a fresh diagnostic tree."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from projects.ukrmol_co.native_debug_patch import CONTRACTED, ORIGINAL_SHA256, instrument

PREFIX = "source/mpi-ci-diag/"
CSF = PREFIX + "Modules/Hamiltonian/CSF_module.F90"
HAMILTONIAN = PREFIX + "Modules/Hamiltonian/Hamiltonian_module.f90"
PINS = {
    CSF: "a8a162fed5b323e01e297b22b520304dd2467e7e197dba209fa06d5045e9736e",
    HAMILTONIAN: "b8e65092a303b60465c5481686fc9f870f81c8590b3452b60f30da77b4770ffe",
    PREFIX + "mpi_scatci.f90": "0017a480c17faf30abbab40ea6cf743237a11e38561d13e92e7ecb804e7776e0",
    CONTRACTED: ORIGINAL_SHA256,
}

AUDIT = """
    subroutine audit_csf_links (manager, values, stage, retained)
        class(CSFManager), target, intent(in) :: manager
        class(CSFObject), target, intent(in) :: values(:)
        character(*), intent(in) :: stage
        class(CSFObject), pointer, optional, intent(in) :: retained(:)
        integer :: i, j, ndet
        if (present(retained)) then
            if (.not. associated(retained, values)) then
                write(stdout,*) 'CO_PTR_FAIL retained-array', trim(stage), 'rank=', grid % grank
                flush(stdout)
                error stop 86
            end if
        end if
        ndet = size(manager % coeffdtrs)
        do i = 1, size(values)
            if (.not. associated(values(i) % manager, manager)) then
                write(stdout,*) 'CO_PTR_FAIL manager', trim(stage), 'rank=', grid % grank, 'csf=', i
                flush(stdout)
                error stop 86
            end if
            if (.not. allocated(values(i) % orbital)) error stop 86
            if (size(values(i) % orbital) /= values(i) % num_orbitals) error stop 86
            if (values(i) % id /= i) error stop 86
            do j = 1, values(i) % num_orbitals
                if (.not. associated(values(i) % orbital(j) % csf, values(i))) then
                    write(stdout,*) 'CO_PTR_FAIL parent', trim(stage), 'rank=', grid % grank, &
                        'csf=', i, 'det=', j, 'actual=', loc(values(i) % orbital(j) % csf), &
                        'expected=', loc(values(i))
                    flush(stdout)
                    error stop 86
                end if
                if (values(i) % orbital(j) % id_dtr < 1 .or. &
                    values(i) % orbital(j) % id_dtr > ndet) error stop 86
            end do
        end do
        write(stdout,*) 'CO_PTR_PASS', trim(stage), 'rank=', grid % grank, &
            'csfs=', size(values), 'determinants=', ndet
        flush(stdout)
    end subroutine audit_csf_links

"""


def replace(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise ValueError(f"No unique source anchor: {old[:80]!r}")
    return text.replace(old, new)


def patches(originals: dict[str, str]) -> dict[str, str]:
    """Keep computational statements intact; assert identity before pointer dereferences."""
    csf = PREFIX + "Modules/Hamiltonian/CSF_module.F90"
    result = dict(originals)
    result[csf] = replace(
        result[csf],
        "    public CSFObject, CSFManager, CSFOrbital",
        "    public CSFObject, CSFManager, CSFOrbital, audit_csf_links",
    )
    result[csf] = replace(result[csf], "\ncontains\n", "\ncontains\n" + AUDIT)
    main = PREFIX + "mpi_scatci.f90"
    result[main] = replace(
        result[main],
        "only: CSFObject, CSFOrbital, CSFManager",
        "only: CSFObject, CSFOrbital, CSFManager, audit_csf_links",
    )
    result[main] = replace(
        result[main],
        '        call master_timer % stop_timer("Construct CSFs")',
        "        call audit_csf_links(Configuration_manager, CSFS, 'after-create')\n"
        '        call master_timer % stop_timer("Construct CSFs")',
    )
    anchor = "            call enrgmx % initialize(ci_rmat)"
    result[main] = replace(
        result[main],
        anchor,
        anchor + "\n"
        "            call audit_csf_links(Configuration_manager, CSFS, 'before-build', "
        "enrgmx % csfs)\n"
        "            if (loc(enrgmx % rmat_ci(1)) /= loc(ci_rmat(1))) error stop 86\n"
        "            if (loc(enrgmx % orbitals) /= loc(Orbitals)) error stop 86\n"
        "            if (loc(enrgmx % options) /= loc(SCATCI_input % opts(i))) error stop 86\n"
        "            write(stdout,*) 'CO_PTR_PASS retained-context rank=', process_grid % grank\n"
        "            flush(stdout)",
    )
    contracted = instrument(result[CONTRACTED])
    anchor = (
        "            call this % evaluate_class_1(i1, "
        "this % options % num_target_state_sym(i1), matrix_elements)"
    )
    contracted = replace(
        contracted,
        anchor,
        "            write(stdout,*) 'CO_PTR_CLASS1_BEGIN rank=', grid % grank, 'target=', i1\n"
        "            flush(stdout)\n" + anchor + "\n"
        "            write(stdout,*) 'CO_PTR_CLASS1_END rank=', grid % grank, 'target=', i1\n"
        "            flush(stdout)",
    )
    result[CONTRACTED] = contracted
    ham = PREFIX + "Modules/Hamiltonian/Hamiltonian_module.f90"
    result[ham] = replace(
        result[ham],
        "module Hamiltonian_module\n",
        "module Hamiltonian_module\n\n    use const_gbl, only: stdout\n",
    )
    for side, index, anchor in (
        ("a", "dtr_idx_a", "        do dtr_idx_a = 1, csf_a % num_orbitals"),
        ("b", "dtr_idx_b", "            do dtr_idx_b = 1, csf_b % num_orbitals"),
    ):
        prefix = f"csf_{side} % orbital({index})"
        block = (
            f"\n                if (loc({prefix} % csf) /= loc(csf_{side})) then\n"
            f"                    write(stdout,*) 'CO_PTR_FAIL slater-{side} rank=', "
            "grid % grank, &\n"
            f"                        'csf=', csf_{side} % id, 'det=', {index}, &\n"
            f"                        'actual=', loc({prefix} % csf), "
            f"'expected=', loc(csf_{side}), &\n"
            f"                        'id_dtr=', {prefix} % id_dtr\n"
            "                    flush(stdout)\n"
            "                    error stop 86\n"
            "                end if"
        )
        result[ham] = replace(result[ham], anchor, anchor + block)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    originals = {}
    for name, digest in PINS.items():
        data = (args.source / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"Unrecognized or already patched native source: {name}")
        originals[name] = data.decode()
    changed = patches(originals)
    record = {
        "scope": "CSF identity/index assertions and progress markers; no engine correction",
        "files": {},
    }
    for name, text in changed.items():
        data = text.encode()
        (args.source / name).write_bytes(data)
        record["files"][name] = {
            "original_sha256": PINS[name],
            "instrumented_sha256": hashlib.sha256(data).hexdigest(),
        }
    with args.output.open("x") as stream:
        stream.write(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
