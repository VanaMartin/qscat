"""Authenticate saved aug-TZ target reuse and identify missing64-state properties."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.stretched_basis_feasibility import (
    SECTORS,
    native_tokens,
    property_inventory,
    read_cidata,
)
from projects.ukrmol_co.verify_stretched_properties import files_in, native_energies, require, sha

RUN = (
    "runs/co-r2500-augtz-sa12-slepc-roots8-spaces80-240-qualified-target"
    "-mpi12-framequalified-20261010"
)
SEED = "runs/co-r2500-augtz-sa12-qc-dz80-residual-qualified-start"
OWNER = "prepared/supported-import-frame-source-20261010"
SUPPORTED_SHA = "3513b97900945be5cdd33366d592d91003252fcbd4b7c0f24f91454ab6acab5e"


def checked_files(archive, collection):
    require(
        sha(archive) == collection["archive"]["sha256"]
        and len(archive) == collection["archive"]["bytes"],
        "Saved archive bytes/hash",
    )
    files = files_in(archive)
    require(set(files) == set(collection["files"]), "Exact captured input inventory")
    for name, pin in collection["files"].items():
        require(
            len(files[name]) == pin["bytes"] and sha(files[name]) == pin["sha256"],
            "Changed saved input " + name,
        )
    return files


def supported(archive):
    require(sha(archive) == SUPPORTED_SHA, "Published supported-space proof")
    prefix = "state-averaged-evidence/"
    files = {n.removeprefix(prefix): v for n, v in files_in(archive).items()}
    pins = json.loads(files["file-index.json"])["files"]
    require(set(files) == set(pins) | {"file-index.json"}, "Supported proof inventory")
    for name, pin in pins.items():
        require(
            sha(files[name]) == pin["sha256"] and len(files[name]) == pin["bytes"],
            "Changed public proof member " + name,
        )
    result = json.loads(files["data/result.json"])
    (scan,) = result["full_scans"]
    require(
        scan["run"] == Path(SEED).name
        and scan["all_roots_pass"]
        and result["all_probes_complete"]
        and len(scan["probes"]) == 48,
        "Supported scan",
    )
    combinations = {(p["spin"], p["irrep"], p["roots"], p["space"]) for p in scan["probes"]}
    require(
        combinations
        == {
            (s, i, n, k)
            for s in [0, 2]
            for i in ["A1", "B1", "B2", "A2"]
            for n in [5, 8]
            for k in [80, 160, 240]
        },
        "Exact48 supported probes",
    )
    groups = {}
    for probe in scan["probes"]:
        require(
            probe["checkpoint_sha256"] == sha(files["qc/" + Path(SEED).name + "/co.casscf.chk"])
            and probe["cycles"] == 600
            and probe["all_roots_pass"]
            and all(probe["ci_converged"])
            and probe["ci_vector_orthogonality_error"] <= 1e-9,
            "Original supported flags/checkpoint",
        )
        require(len(probe["states"]) == probe["roots"], "Probe root inventory")
        for state in probe["states"]:
            require(np.isfinite(list(state.values())).all(), "Nonfinite saved probe")
            require(
                abs(state["spin_square"] - (0 if probe["spin"] == 0 else 2)) <= 1e-6
                and all(
                    0 <= state[k] <= 1e-9
                    for k in [
                        "physical_residual_norm_hartree",
                        "penalized_residual_norm_hartree",
                        "spin_penalty_action_norm_hartree",
                    ]
                ),
                "Saved independent action/spin gate",
            )
        groups[(probe["spin"], probe["irrep"], probe["roots"], probe["space"])] = np.array(
            [s["energy_hartree"] for s in probe["states"]]
        )
    for spin, irrep in [(s, i) for s in [0, 2] for i in ["A1", "B1", "B2", "A2"]]:
        for count in [5, 8]:
            require(
                np.ptp([groups[(spin, irrep, count, k)] for k in [80, 160, 240]], axis=0).max()
                <= 1e-7,
                "Supported trial-space energy stability",
            )
        for space in [80, 160, 240]:
            np.testing.assert_allclose(
                groups[(spin, irrep, 5, space)],
                groups[(spin, irrep, 8, space)][:5],
                atol=1e-7,
                rtol=0,
            )
    owner = json.loads(files["owner-source/execution.json"])
    require(
        owner["original_controller_exit"] == 1
        and [s["original_exit"] for s in owner["steps"]] == [1, 0, 0],
        "Original supported owner/outer failure preserved",
    )
    return files, groups


def reconstruct(archive, collection, public_proof):
    files = checked_files(archive, collection)
    proof, groups = supported(public_proof)

    def load(name):
        return json.loads(files[name])

    frozen = load(OWNER + "/frozen-inputs.json")
    owner = load(OWNER + "/execution.json")
    require(
        sha(files[OWNER + "/frozen-inputs.json"]) == owner["frozen_inputs_sha256"],
        "Original import frozen index",
    )
    original_pins = 0
    for category in ["source", "inputs"]:
        for name, pin in frozen[category].items():
            path = OWNER + "/" + name if category == "source" else name
            if path in files:
                require(sha(files[path]) == pin, "Original frozen input/source " + name)
                original_pins += 1
    for name in ["config.json", "target.json", "target-diagnostics.json", "co.casscf.chk"]:
        require(
            files[SEED + "/" + name] == proof["qc/" + Path(SEED).name + "/" + name],
            "Actual covered seed " + name,
        )
    cfg, target, diagnostics = [
        load(RUN + "/" + n) for n in ["config.json", "target.json", "target-diagnostics.json"]
    ]
    require(
        cfg["basis"] == "aug-cc-pVTZ"
        and cfg["bond_length"] == 2.5
        and cfg["active_orbitals"] == [6, 3, 3, 0]
        and cfg["virtual_orbitals"] == [0] * 4
        and cfg["sa_singlet_roots"] == cfg["sa_triplet_roots"] == [5] * 4
        and cfg["target_roots"] == 8
        and cfg["target_states_used"] == 40,
        "Common40 objective/CAS12",
    )
    for t, d in [
        (target, diagnostics),
        (load(SEED + "/target.json"), load(SEED + "/target-diagnostics.json")),
    ]:
        require(
            t["converged"]
            and d["orbital_converged"]
            and d["ci_converged"] == [True] * 8
            and d["fresh_ci_converged"] == [True] * 8
            and len(t["states"]) == 40
            and d["iterations"][-1]["orbital_gradient_norm"] <= 1e-7
            and abs(d["iterations"][-1]["energy_change_hartree"]) <= 1e-11,
            "Ordinary objective convergence",
        )
    geometry = RUN + "/output/CO/geom1/"
    require(
        sha(files[OWNER + "/aug/projects/ukrmol_co/target.py"]) == target["source_sha256"]
        and sha(files[geometry + "co.molden"]) == target["molden_sha256"],
        "Executed target source/Molden",
    )
    native = read_cidata(files[geometry + "fort.26"])
    independent, norm_error = native_energies(files[geometry + "fort.26"])
    np.testing.assert_allclose(
        [native[s]["energies"] for s in SECTORS],
        [independent[s] for s in SECTORS],
        atol=1e-12,
        rtol=0,
    )
    errors = {}
    for index, sector in enumerate(SECTORS, 1):
        spin, irrep = sector.split(".")
        number, tokens = native_tokens(
            files[geometry + "outputs/target.scatci." + sector + ".out"].decode()
        )
        require(number == index, "Own native set order")
        np.testing.assert_allclose(tokens, native[sector]["energies"], atol=1e-9, rtol=0)
        error = float(
            abs(
                native[sector]["energies"] - groups[(0 if spin == "singlet" else 2, irrep, 8, 160)]
            ).max()
        )
        require(error <= 1e-7, "Independent all-eight target roots")
        errors[sector] = error
    verification = load(OWNER + "/" + Path(RUN).name + ".verification.json")
    historical = json.loads(
        Path("projects/ukrmol_co/progress-completion-20261010T2045.json").read_text()
    )
    require(
        verification == historical["augtz"]["verification"],
        "Unchanged independently recorded import verification",
    )
    require(
        verification["all_64_roots_pass"]
        and verification["imported_checkpoint_sha256"] == sha(files[geometry + "co.casscf.chk"])
        and verification["reference_checkpoint_sha256"] == sha(files[SEED + "/co.casscf.chk"])
        and verification["verifier_source_sha256"]
        == sha(files[OWNER + "/aug/verify-cas12-native-import.py"])
        and min(v for vs in verification["subspace_overlap_singular_values"].values() for v in vs)
        >= 1 - 1e-8,
        "Saved independent physical frame qualification",
    )
    old = load(SEED + "/target.json")
    seed_error = float(
        abs(
            np.array([s["energy_hartree"] for s in target["states"]])
            - [s["energy_hartree"] for s in old["states"]]
        ).max()
    )
    require(seed_error <= 1e-7, "Original seed-preservation gate")
    properties = property_inventory(files[geometry + "fort.24"].decode(), native)
    property_lines = [line.split() for line in files[geometry + "fort.24"].decode().splitlines()]
    dipole = -next(
        float(f[8].replace("D", "E"))
        for f in property_lines
        if len(f) > 8 and f[:8] == ["1", "1", "0", "1", "0", "3", "1", "0"]
    )
    require(abs(dipole - target["ground_dipole_au"][2]) <= 1e-5, "Native ground dipole")
    baseline = json.loads(
        Path("projects/ukrmol_co/stretched-basis-properties-contract-20261026.json").read_text()
    )
    require(len(collection["retained_csf_hashes"]) == 8, "Eight retained CSF definitions")
    for name, pin in collection["retained_csf_hashes"].items():
        relative = "output/CO/geom1/" + Path(name).name
        require(
            pin["sha256"] == baseline["models"]["tz"]["input_sha256"][relative],
            "Byte-identical qualified CAS12 CSF " + name,
        )
    stages = [line.split("\t") for line in files[RUN + "/stages.tsv"].decode().splitlines()]
    require(
        len(stages) == 19 and all(s[-1] == "0" for s in stages), "Original19 native stage exits"
    )
    require(
        owner["original_controller_exit"] == 1
        and [s["original_exit"] for s in owner["steps"]] == [0, 0, 0, 1],
        "Original import/analysis/verifier successes and later queue failure",
    )
    resources = load(RUN + "/resources.json")
    samples = list(csv.DictReader(io.StringIO(files[RUN + "/resources.csv"].decode())))
    require(
        resources["exit_code"] == 0
        and max(int(s["container_memory_bytes"]) for s in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Saved actual import resources",
    )
    return dict(
        reuse_feasibility_verdict=True,
        source_inputs_verified=len(files),
        original_frozen_pins_verified=original_pins,
        supported_public_payloads_verified=269,
        supported_saved_probes=48,
        original_coverage_owner_exit=1,
        original_import_owner_exit=1,
        original_import_analysis_verifier_exits=[0, 0, 0],
        native_root_count=64,
        maximum_binary_to_supported_root_error_hartree=max(errors.values()),
        native_vector_norm_error=norm_error,
        maximum_seed_preservation_error_hartree=seed_error,
        imported_checkpoint_sha256=verification["imported_checkpoint_sha256"],
        covered_checkpoint_sha256=verification["reference_checkpoint_sha256"],
        saved_physical_frame_singular_values=verification["subspace_overlap_singular_values"],
        csf_files_byte_identical_to_qualified_tz=8,
        properties=properties,
        native_ground_dipole_z_au=dipole,
        original_denprop_seconds=float(next(s[-2] for s in stages if s[1] == "denprop")),
        original_import_resources=resources,
        original_import_resource_samples=len(samples),
        missing_property_exports=1,
        missing_new_preparation_and_B1_2048_export=True,
        new_ci_calls=0,
        new_qc_calls=0,
        new_native_calls=0,
        new_maps=0,
        new_hamiltonian_actions=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        scope=(
            "Saved native64 energies/norms/properties reconstructed against public supported-space "
            "probes;physical action/subspace metrics are authenticated saved evidence,not fresh "
            "contractions;no TZ/aug-TZ target correspondence,resonance identity,scattering effect "
            "or accuracy verdict"
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["archive", "collection", "supported_proof"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            reconstruct(
                args.archive.read_bytes(),
                json.loads(args.collection.read_text()),
                args.supported_proof.read_bytes(),
            ),
            indent=2,
            allow_nan=False,
        )
    )
