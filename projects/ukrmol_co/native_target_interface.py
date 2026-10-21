"""Bounded native target I/O over saved CIDATA and integral files; no solver calls."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

from projects.ukrmol_co.a2_channel_preflight import stage


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compiler_environment(environment):
    """Override stale compiler relocation paths exactly as source-build.sh does."""
    return environment | {
        "OMPI_FC": "/opt/compiler/bin/gfortran",
        "OMPI_CC": "/opt/compiler/bin/gcc",
    }


def load_contract(root, path):
    """An environment-only successor inherits every original scientific setting."""
    contract = json.loads(path.read_text())
    if "parent_contract" not in contract:
        return contract
    parent = root / contract["parent_contract"]
    if digest(parent) != contract["input_sha256"][contract["parent_contract"]]:
        raise ValueError("Changed original native-interface contract")
    base = json.loads(parent.read_text())
    return base | {
        "input_sha256": base["input_sha256"] | contract["input_sha256"],
        "predecessor_execution": contract["predecessor_execution"],
        "predecessor_failure": contract["predecessor_failure"],
        "predecessor_compile_log": contract["predecessor_compile_log"],
    }


def integral_block(path):
    """Read the seven-record diagnostic format, not the original GBTOlib format."""
    with FortranFile(path) as stream:
        header = stream.read_record("<i8")
        if len(header) != 3 or header[0] != 13 or header[2] not in (2, 3):
            raise ValueError("Unexpected target integral diagnostic header")
        indices = stream.read_record("<i8")
        if len(indices) != 21:
            raise ValueError("Missing native inventory/target addresses")
        geometry = stream.read_record("<f8").reshape((4, header[2]), order="F").T
        overlap = stream.read_record("<f8").reshape((13, 13), order="F")
        h = stream.read_record("<f8").reshape((13, 13), order="F")
        eri = stream.read_record("<f8").reshape((13,) * 4, order="F")
        nuclear = stream.read_record("<f8")
        if len(nuclear) != 1 or stream._fp.tell() != path.stat().st_size:
            raise ValueError("Incomplete or trailing native integral diagnostic")
    return dict(
        total=int(header[1]),
        counts=indices[:8],
        mapping=indices[8:],
        geometry=geometry,
        overlap=overlap,
        h=h,
        eri=eri,
        nuclear=float(nuclear[0]),
    )


def compare_integrals(source, full):
    """Conservative common-basis operator bound, not a scattering error estimate."""
    for block, counts, mapping in (
        (source, [7, 3, 3, 0, 0, 0, 0, 0], list(range(1, 14))),
        (full, [58, 37, 37, 22, 0, 0, 0, 0], [*range(1, 8), 59, 60, 61, 96, 97, 98]),
    ):
        if block["total"] != sum(counts):
            raise ValueError("Unexpected total native orbital count")
        np.testing.assert_array_equal(block["counts"], counts)
        np.testing.assert_array_equal(block["mapping"], mapping)
        for name in ("geometry", "overlap", "h", "eri", "nuclear"):
            if not np.isfinite(block[name]).all():
                raise ValueError(f"Nonfinite saved integral {name}")
        np.testing.assert_allclose(block["overlap"], np.eye(13), atol=1e-10, rtol=0)
        np.testing.assert_allclose(block["h"], block["h"].T, atol=1e-10, rtol=0)
        for permutation in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
            np.testing.assert_allclose(
                block["eri"], block["eri"].transpose(permutation), atol=1e-10, rtol=0
            )
    np.testing.assert_allclose(source["geometry"][:, 3], [6, 8], atol=0, rtol=0)
    physical = full["geometry"][full["geometry"][:, 3] != 0]
    if len(full["geometry"]) != 3 or len(physical) != 2:
        raise ValueError("Missing zero-charge continuum centre")
    np.testing.assert_allclose(physical, source["geometry"], atol=1e-12, rtol=0)
    np.testing.assert_allclose(source["overlap"], full["overlap"], atol=1e-10, rtol=0)
    dh, dv = [np.abs(full[key] - source[key]) for key in ("h", "eri")]
    nuclear_delta = abs(full["nuclear"] - source["nuclear"])
    if max(dh.max(), dv.max()) > 1e-10:
        raise ValueError("Target integral element difference exceeds 1e-10 Hartree")
    bound = nuclear_delta + 2 * dh.sum() + 2 * dv.sum()
    if bound > 1e-7:
        raise ValueError("Target Hamiltonian operator-difference bound exceeds 1e-7 Hartree")
    return dict(
        target_integral_equivalence_verdict=True,
        one_electron_entries=169,
        overlap_entries=169,
        eri_ordered_entries=28561,
        max_h_difference_hartree=float(dh.max()),
        max_eri_difference_hartree=float(dv.max()),
        h_difference_l1_hartree=float(dh.sum()),
        eri_difference_l1_hartree=float(dv.sum()),
        nuclear_difference_hartree=float(nuclear_delta),
        operator_bound_hartree=float(bound),
        max_source_overlap_defect=float(np.max(abs(source["overlap"] - np.eye(13)))),
        max_full_overlap_defect=float(np.max(abs(full["overlap"] - np.eye(13)))),
        overlap_change=float(np.max(abs(full["overlap"] - source["overlap"]))),
        scope="Difference of saved finite target Hamiltonians under common orthonormal "
        "basis identification, not absolute accuracy or a scattering error bound.",
    )


def compare_target(native, cidata, csfdir, expected):
    """Compare the actual CIRMAT output with all saved roots and qualified signs."""
    reports = []
    units = (710, 730, 711, 731, 712, 732, 713, 733)
    with FortranFile(native) as result, FortranFile(cidata) as source:
        for index, want in enumerate(expected):
            source.read_record("u1")
            header = source.read_record("u1").tobytes()
            centers, dimension, roots = np.frombuffer(header, "<i8", count=3, offset=136)
            if centers != 2 or (dimension, roots) != (want["dimension"], want["roots"]):
                raise ValueError("Changed native source CIDATA inventory")
            for _ in range(centers):
                source.read_record("u1")
            energy_record = source.read_record("u1").tobytes()
            source_phases = np.frombuffer(energy_record, "<i8", count=dimension)
            energies = np.frombuffer(energy_record, "<f8", count=roots, offset=8 * dimension)
            with FortranFile(csfdir / f"fort.{units[index]}") as csfs:
                csfs.read_record("u1")
                full_phases = csfs.read_record("<i8")
            signs = source_phases * full_phases
            if (
                hashlib.sha256(signs.astype("<i8").tobytes()).hexdigest()
                != want["sign_int64_sha256"]
            ):
                raise ValueError("Native input phase map differs from qualified correspondence")
            vectors = []
            for root in range(roots):
                record = source.read_record("u1").tobytes()
                if np.frombuffer(record, "<i8", count=1)[0] != root + 1:
                    raise ValueError("Changed native source root ID")
                vectors.append(np.frombuffer(record, "<f8", offset=8))
            np.testing.assert_array_equal(result.read_record("<i8"), [index + 1, dimension, roots])
            native_energies = result.read_record("<f8")
            native_vectors = result.read_record("<f8")
            mapped = np.asarray(vectors) * signs
            if (
                native_energies.tobytes() != energies.tobytes()
                or native_vectors.tobytes() != mapped.tobytes()
            ):
                raise ValueError("Native CIRMAT differs from saved energy/phase-rule coefficients")
            reports.append(
                dict(
                    sector=want["sector"],
                    roots=int(roots),
                    dimension=int(dimension),
                    exact_native_coefficients=True,
                    exact_native_energies=True,
                    source_nuclear_reference_hartree=float(
                        np.frombuffer(header, "<f8", count=1, offset=192)[0]
                    ),
                    native_vector_sha256=hashlib.sha256(native_vectors.tobytes()).hexdigest(),
                )
            )
        if result._fp.tell() != native.stat().st_size or source._fp.tell() != cidata.stat().st_size:
            raise ValueError("Trailing or incomplete native/CIDATA data")
    return reports


def execute(root, output, contract):
    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed native interface input: {name}")
        for name, want in contract["engine_library_sha256"].items():
            if digest(Path(name)) != want:
                raise ValueError(f"Changed installed engine library: {name}")

    verify()
    if "predecessor_execution" in contract:
        predecessor = json.loads((root / contract["predecessor_execution"]).read_text())
        failure = json.loads((root / contract["predecessor_failure"]).read_text())
        assert predecessor["original_owner_exit"] == 1
        assert (
            predecessor["steps"][0]["original_command_exit"]
            == predecessor["steps"][0]["original_container_exit"]
            == 1
        )
        assert [(s["name"], s["original_exit"], s["timed_out"]) for s in failure["stages"]] == [
            ("compile", 1, False)
        ]
        assert (
            "unable to find the specified compiler"
            in (root / contract["predecessor_compile_log"]).read_text()
        )
    qualified = json.loads((root / contract["map_result"]).read_text())
    assert qualified["original"]["signed_target_csf_correspondence_verdict"] is True
    assert qualified["reconstruction"]["canonical_expansion_reconstruction_verdict"] is True
    assert qualified["original"]["roots_verified"] == 74
    output.mkdir(exist_ok=False)
    report = dict(
        stages=[],
        native_consumer_execution_verdict=None,
        target_integral_equivalence_verdict=None,
        scattering_channel_map_verdict=None,
        retained_subspace_dominance_verdict=False,
        new_qc_calls=0,
        target_diagonalizations=0,
        property_exports=0,
        scattering_inner_solves=0,
        input_sha256=contract["input_sha256"],
    )
    empty = output / "empty.inp"
    empty.write_text("")
    try:
        source = Path(__file__).with_name("native_target_io.f90")
        executable = output / "native_target_io"
        command = [
            "mpifort",
            "-fdefault-integer-8",
            "-ffree-line-length-none",
            "-O0",
            "-g",
            "-fcheck=all",
        ]
        command += [
            "-I/opt/ukrmolp/include/UKRmol+/" + p
            for p in ("gbtolib", "ci-diag", "global", "utilities")
        ]
        libraries = ["libscatci.a", "libukplus_global.a", "libukplus_utilities.a", "libGBTO.a"]
        command += [str(source), "-o", str(executable), "-Wl,--start-group"]
        command += ["/opt/ukrmolp/lib.double/" + n for n in libraries]
        command += [
            "-Wl,--end-group",
            "-L/opt/ukrmolp/lib",
            "-lscalapack",
            "-larpack",
            "-lopenblas",
            "-fopenmp",
        ]
        report["compile_command"] = command
        report["compiler_environment"] = compiler_environment({})
        stage(
            output, report, "compile", command, empty, 120, compiler_environment(os.environ.copy())
        )
        report["executable_sha256"] = digest(executable)
        report["module_sha256"] = {
            str(p): digest(p) for p in Path("/opt/ukrmolp/include/UKRmol+").rglob("*.mod")
        }
        common = [str(executable)]
        for key, name in (("source_integrals", "source.moints"), ("full_integrals", "full.moints")):
            shutil.copy2(root / contract[key], output / name)
            if digest(output / name) != digest(root / contract[key]):
                raise ValueError("Integral read-copy differs from frozen input")
        stages = [
            (
                "target",
                [
                    *common,
                    "target",
                    str(root / contract["cidata"]),
                    str(root / contract["full_csf_directory"]),
                    str(output / "target.bin"),
                ],
            ),
            (
                "source_integrals",
                [
                    *common,
                    "integrals",
                    str(output / "source.moints"),
                    "unused",
                    str(output / "source-integrals.bin"),
                ],
            ),
            (
                "full_integrals",
                [
                    *common,
                    "integrals",
                    str(output / "full.moints"),
                    "unused",
                    str(output / "full-integrals.bin"),
                ],
            ),
        ]
        report["native_commands"] = dict(stages)
        for label, command in stages:
            stage(output, report, label, command, empty, 30, os.environ.copy())
            marker = (
                "CO_NATIVE_TARGET_READ_COMPLETE"
                if label == "target"
                else "CO_NATIVE_INTEGRAL_READ_COMPLETE"
            )
            if (output / f"{label}.out").read_text().count(marker) != 1:
                raise ValueError(f"Missing unique native completion marker: {label}")
        report["target_sectors"] = compare_target(
            output / "target.bin",
            root / contract["cidata"],
            root / contract["full_csf_directory"],
            qualified["original"]["sectors"],
        )
        report["native_consumer_execution_verdict"] = True
        source_block = integral_block(output / "source-integrals.bin")
        report["integral_comparison"] = compare_integrals(
            source_block, integral_block(output / "full-integrals.bin")
        )
        np.testing.assert_allclose(
            [s["source_nuclear_reference_hartree"] for s in report["target_sectors"]],
            source_block["nuclear"],
            atol=1e-12,
            rtol=0,
        )
        report["target_integral_equivalence_verdict"] = True
        for key, name in (("source_integrals", "source.moints"), ("full_integrals", "full.moints")):
            if digest(output / name) != digest(root / contract[key]):
                raise ValueError("Native reader modified its integral read-copy")
        verify()
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except BaseException as error:
        report["error"] = repr(error)
        (output / "failure.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        raise
    finally:
        verify()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    execute(args.root, args.output, load_contract(args.root, args.contract))


if __name__ == "__main__":
    main()
