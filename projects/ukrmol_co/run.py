"""Prepare and profile a CO deck inside the UKRmol+ reference image.

Only the standard library is needed, including on the image's Python 3.11.
Run from the repository root with ``python3 -m projects.ukrmol_co.run``.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import time
import urllib.request
import zipfile
from pathlib import Path

SCRIPTS_URL = "https://zenodo.org/api/records/7851856/files/UKRmol-scripts-release-1.0.zip/content"
SCRIPTS_MD5 = "f2f37391a26280d212eec8393f779830"


def acquire_scripts(directory: Path) -> Path:
    """Download the pinned archive, verify it, and extract into a fresh directory."""
    directory.mkdir(parents=True, exist_ok=True)
    # Cold-cache geometry workers must not share a partial download or extraction.
    with (directory / ".acquire.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        archive = directory / "UKRmol-scripts-release-1.0.zip"
        if not archive.exists():
            temporary = archive.with_suffix(".part")
            urllib.request.urlretrieve(SCRIPTS_URL, temporary)
            temporary.rename(archive)
        if hashlib.md5(archive.read_bytes()).hexdigest() != SCRIPTS_MD5:
            raise ValueError(f"Upstream archive checksum mismatch: {archive}")
        extracted = directory / "UKRmol-scripts-release-1.0"
        if not extracted.exists():
            with zipfile.ZipFile(archive) as bundle:
                bundle.extractall(directory)
    return extracted.resolve()


def instrument_library(library: Path) -> None:
    """Record exact stage wall times and propagate failures suppressed upstream."""
    source = library.read_text()
    before = '    if (system("$command") != 0) {'
    after = """    my $stage_start = Time::HiRes::time();
    my $stage_status = system("$command");
    open(my $stage_log, ">>", $ENV{'UKRMOL_STAGE_LOG'}) or die $!;
    print $stage_log join("\\t", $task, $program, $statespin, $str_ir,
                         Time::HiRes::time() - $stage_start, $stage_status), "\\n";
    close($stage_log);
    if ($stage_status != 0) {"""
    if source.count(before) != 1:
        raise ValueError("UKRmol-scripts run_code layout has changed")
    source = source.replace(before, after)
    before = (
        '    if ($redirected_master_output) { &run_system("$mv_cmd log_file.0 $output", $r_par); }'
    )
    after = (
        """    die "UKRmol stage failed: $command\\n" if $stage_status != 0;
"""
        + before
    )
    if source.count(before) != 1:
        raise ValueError("UKRmol-scripts failure-handling layout has changed")
    library.write_text("use Time::HiRes ();\n" + source.replace(before, after))


def enable_state_average(library: Path) -> None:
    """Replace the single-state upstream QC adapter with the explicit PySCF backend."""
    source = library.read_text()
    before = '$command = "$dir$bs$program$ext_exe --input $input --output $output 2> $error";'
    after = (
        '$command = "python3 -m projects.ukrmol_co.target '
        "--config $ENV{'UKRMOL_RUN_DIR'}/config.json --output $output 2> $error\";"
    )
    if source.count(before) != 1:
        raise ValueError("UKRmol-scripts Psi4 command layout has changed")
    source = source.replace(before, after)
    before = "sub read_psi4_output {\n  my ($r_par) = @_;"
    if source.count(before) != 1:
        raise ValueError("UKRmol-scripts Psi4 reader layout has changed")
    source = source.replace(before, before + "\n  return &read_state_average_output($r_par);")
    library.write_text(source + "\n" + Path(__file__).with_name("state_average.pm").read_text())


def select_target_diagonalizer(
    template: Path, solver: str, tolerance: float, max_cycles: int
) -> None:
    """Select the pinned target path explicitly without changing scattering."""
    igh, force_serial = {"davidson-serial": (0, 1), "slepc": (-1, 0)}[solver]
    source = template.read_text()
    before = ">>>IGHT<<<igh = >>>IGH<<<,"
    if source.count(before) != 1:
        raise ValueError("No unique target diagonalizer setting in the pinned template")
    template.write_text(
        source.replace(
            before,
            f"  igh = {igh}, forse = {force_serial}, "
            f"crite = {tolerance:.16g}, maxiter = {max_cycles},",
        )
    )


def allocated_bytes(directory: Path) -> int:
    """Count allocated file blocks once, excluding symlinks and hardlink duplicates."""
    seen = set()
    total = 0
    for root, _, files in os.walk(directory):
        for name in files:
            try:
                entry = (Path(root) / name).lstat()
            except FileNotFoundError:
                continue
            key = (entry.st_dev, entry.st_ino)
            if key not in seen:
                seen.add(key)
                total += entry.st_blocks * 512
    return total


def cgroup_bytes(name: str) -> int:
    """Read aggregate container memory; return zero when cgroup v2 is unavailable."""
    path = Path("/sys/fs/cgroup") / name
    return int(path.read_text()) if path.exists() else 0


def cgroup_anon_bytes() -> int:
    """Read anonymous memory separately from the container's filesystem page cache."""
    path = Path("/sys/fs/cgroup/memory.stat")
    if not path.exists():
        return 0
    values = dict(line.split() for line in path.read_text().splitlines())
    return int(values["anon"])


def cgroup_cpu_seconds() -> float:
    """Read aggregate CPU time used by all processes in the container."""
    path = Path("/sys/fs/cgroup/cpu.stat")
    if not path.exists():
        return 0.0
    values = dict(line.split() for line in path.read_text().splitlines())
    return int(values["usage_usec"]) / 1e6


def run_profiled(workdir: Path, env: dict[str, str], command: list[str] | None = None) -> int:
    """Sample container memory and persistent disk usage throughout the calculation."""
    started = time.monotonic()
    baseline_memory = cgroup_bytes("memory.current")
    baseline_cpu = cgroup_cpu_seconds()
    peak_memory = peak_anon = peak_disk = peak_scratch = 0
    disk = scratch = 0
    next_disk = 0.0
    with (
        (workdir / "run.log").open("w") as log,
        (workdir / "resources.csv").open("w") as resource_log,
    ):
        writer = csv.writer(resource_log)
        writer.writerow(
            [
                "elapsed_s",
                "container_memory_bytes",
                "anonymous_memory_bytes",
                "run_disk_bytes",
                "scratch_bytes",
            ]
        )
        process = subprocess.Popen(
            ["perl", "main.pl", "co.pl"] if command is None else command,
            cwd=workdir,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        while True:
            elapsed = time.monotonic() - started
            memory = cgroup_bytes("memory.current")
            anon = cgroup_anon_bytes()
            if elapsed >= next_disk:
                disk = allocated_bytes(workdir)
                scratch = allocated_bytes(workdir / "scratch")
                next_disk = elapsed + 2.0
            peak_memory = max(peak_memory, memory)
            peak_anon = max(peak_anon, anon)
            peak_disk = max(peak_disk, disk)
            peak_scratch = max(peak_scratch, scratch)
            writer.writerow([round(elapsed, 3), memory, anon, disk, scratch])
            resource_log.flush()
            if process.poll() is not None:
                break
            time.sleep(0.2)
    peak_disk = max(peak_disk, allocated_bytes(workdir))
    summary = {
        "exit_code": process.returncode,
        "wall_seconds": time.monotonic() - started,
        "sampled_peak_container_memory_bytes": peak_memory,
        "baseline_container_memory_bytes": baseline_memory,
        "sampled_peak_container_anon_bytes": peak_anon,
        "kernel_peak_container_memory_bytes": cgroup_bytes("memory.peak"),
        "container_cpu_seconds": cgroup_cpu_seconds() - baseline_cpu,
        "sampled_peak_run_allocated_bytes": peak_disk,
        "sampled_peak_scratch_allocated_bytes": peak_scratch,
        "memory_sample_interval_seconds": 0.2,
        "disk_sample_interval_seconds": 2.0,
    }
    (workdir / "resources.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return process.returncode


def main() -> None:
    """Create a new, self-contained run directory and execute its CO input deck."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--scripts-cache", type=Path, default=Path("/work/upstream"))
    parser.add_argument("--bond-length", type=float, default=2.1323, help="C-O distance in bohr")
    parser.add_argument("--model", choices=["SE", "SEP", "CAS-A"], default="SEP")
    parser.add_argument("--orbitals", choices=["HF", "natural", "state-averaged"], default="HF")
    parser.add_argument(
        "--target-only", action="store_true", help="Run QC and UKRmol target checks"
    )
    parser.add_argument(
        "--qc-only", action="store_true", help="Run state-averaged QC checks before UKRmol"
    )
    parser.add_argument("--sa-singlet-roots", type=int, nargs=4, default=[5, 5, 5, 5])
    parser.add_argument("--sa-triplet-roots", type=int, nargs=4, default=[5, 5, 5, 5])
    parser.add_argument("--target-memory-mb", type=int, default=8000)
    parser.add_argument("--target-max-cycles", type=int, default=100)
    parser.add_argument("--target-energy-tolerance", type=float, default=1e-9)
    parser.add_argument("--target-gradient-tolerance", type=float, default=1e-5)
    parser.add_argument("--target-ci-tolerance", type=float, default=1e-10)
    parser.add_argument("--target-ci-residual-tolerance", type=float)
    parser.add_argument("--target-ci-fresh-start", action="store_true")
    parser.add_argument("--target-ci-lindep", type=float, default=1e-14)
    parser.add_argument(
        "--target-optimizer",
        choices=["one-step", "newton"],
        default="one-step",
        help="Orbital optimizer; newton is an experimental diagnostic",
    )
    parser.add_argument("--target-ah-lindep", type=float, default=1e-14)
    parser.add_argument("--target-ah-tolerance", type=float, default=1e-12)
    parser.add_argument("--target-ah-start-tolerance", type=float, default=2.5)
    parser.add_argument("--target-verbosity", type=int, choices=range(10), default=4)
    parser.add_argument(
        "--target-initial-checkpoint",
        type=Path,
        help="Project core/active MOs from a retained PySCF CASSCF checkpoint",
    )
    parser.add_argument(
        "--active-orbitals", type=int, nargs=4, help="CAS active counts in A1,B1,B2,A2"
    )
    parser.add_argument("--target-roots", type=int, default=5, help="CAS roots per spin/irrep")
    parser.add_argument(
        "--target-singlet-roots", type=int, nargs=4, help="Override A1,B1,B2,A2 roots"
    )
    parser.add_argument(
        "--target-triplet-roots", type=int, nargs=4, help="Override A1,B1,B2,A2 roots"
    )
    parser.add_argument("--target-states-used", type=int, help="Lowest CAS target states retained")
    parser.add_argument(
        "--congen-workspace",
        type=int,
        default=200000,
        help="CONGEN NDIMX entries; CDIMX and NODIMX are set to one tenth of this size",
    )
    parser.add_argument("--scatci-memory-gib", type=float, default=2.5)
    parser.add_argument(
        "--target-diagonalizer",
        choices=["auto", "davidson-serial", "slepc"],
        default="auto",
        help="UKRmol target solver; selected-root paths are experimental",
    )
    parser.add_argument("--target-diagonalizer-tolerance", type=float, default=1e-12)
    parser.add_argument("--target-diagonalizer-max-cycles", type=int, default=500)
    parser.add_argument(
        "--basis",
        choices=["cc-pVDZ", "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"],
        default="aug-cc-pVDZ",
    )
    parser.add_argument(
        "--frozen-orbitals",
        type=int,
        choices=[2, 4],
        default=2,
        help="Freeze the two 1s orbitals, or also the two occupied 2s-derived sigma orbitals",
    )
    parser.add_argument("--virtual-orbitals", type=int, nargs=4, default=[10, 4, 4, 2])
    parser.add_argument("--radius", type=int, choices=[10, 13, 15, 18], default=18)
    parser.add_argument("--maxl", type=int, choices=range(1, 7), default=4)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--precision", choices=["double", "quad"], default="double")
    parser.add_argument("--deletion-threshold", type=float, default=1e-7)
    parser.add_argument(
        "--propagation-radius", type=float, default=100.0, help="Outer radius in bohr"
    )
    parser.add_argument("--energy-start-ev", type=float, default=0.1)
    parser.add_argument("--energy-step-ev", type=float, default=0.01)
    parser.add_argument("--energies", type=int, default=491)
    args = parser.parse_args()
    if args.qc_only and (args.target_only or args.orbitals != "state-averaged"):
        parser.error("QC-only requires state-averaged orbitals and excludes target-only")
    if args.target_diagonalizer != "auto" and (args.qc_only or args.model != "CAS-A"):
        parser.error("Selected-root diagonalization requires a UKRmol CAS target")
    positive = (
        args.bond_length,
        args.ranks,
        args.energy_start_ev,
        args.energy_step_ev,
        args.energies,
        args.deletion_threshold,
        args.propagation_radius,
        args.congen_workspace,
        args.scatci_memory_gib,
        args.target_diagonalizer_tolerance,
        args.target_diagonalizer_max_cycles,
        args.target_memory_mb,
        args.target_max_cycles,
        args.target_energy_tolerance,
        args.target_gradient_tolerance,
        args.target_ci_tolerance,
        args.target_ci_lindep,
        args.target_ah_lindep,
        args.target_ah_tolerance,
        args.target_ah_start_tolerance,
    )
    if not all(math.isfinite(value) and value > 0 for value in positive):
        parser.error("Bond length, MPI ranks, deletion threshold and energy grid must be positive")
    if args.target_ci_residual_tolerance is not None and not (
        math.isfinite(args.target_ci_residual_tolerance) and args.target_ci_residual_tolerance > 0
    ):
        parser.error("CI residual tolerance must be finite and positive")
    if args.congen_workspace < 1000:
        parser.error("CONGEN workspace must be at least 1000 entries")
    if min(args.virtual_orbitals) < 0 or args.virtual_orbitals[1] != args.virtual_orbitals[2]:
        parser.error("Virtual counts must be nonnegative with equal B1/B2 counts")
    if args.model == "CAS-A":
        if args.frozen_orbitals != 2:
            parser.error("CAS diagnostics use two frozen core orbitals and ten active electrons")
        args.active_orbitals = args.active_orbitals or [4, 2, 2, 0]
        if (
            min(args.active_orbitals) < 0
            or args.active_orbitals[0] < 3
            or args.active_orbitals[1] < 1
            or args.active_orbitals[1] != args.active_orbitals[2]
        ):
            parser.error("CAS must contain the occupied valence orbitals and equal B1/B2 spaces")
        if args.target_roots < 1:
            parser.error("CAS target roots must be positive")
        args.target_singlet_roots = args.target_singlet_roots or [args.target_roots] * 4
        args.target_triplet_roots = args.target_triplet_roots or [args.target_roots] * 4
        for counts in (args.target_singlet_roots, args.target_triplet_roots):
            if min(counts) < 0 or counts[1] != counts[2]:
                parser.error("Target roots must be nonnegative with equal B1/B2 counts")
        if args.target_singlet_roots[0] < 1:
            parser.error("Target roots must include the A1 singlet ground state")
        computed_states = sum(args.target_singlet_roots + args.target_triplet_roots)
        if args.target_states_used is None:
            args.target_states_used = computed_states
        if not 1 <= args.target_states_used <= computed_states:
            parser.error("CAS target-state count must be positive and not exceed computed roots")
    else:
        if (
            args.active_orbitals is not None
            or args.orbitals != "HF"
            or args.target_singlet_roots is not None
            or args.target_triplet_roots is not None
        ):
            parser.error("SE/SEP use the occupied HF space; active counts apply only to CAS")
        args.active_orbitals = [5 - args.frozen_orbitals, 1, 1, 0]
        args.target_singlet_roots = [1, 0, 0, 0]
        args.target_triplet_roots = [0, 0, 0, 0]
        args.target_states_used = 1
    if args.radius == 18 and args.maxl > 5:
        parser.error("The upstream radius-18 Gaussian continuum is available only through l=5")
    if args.propagation_radius <= args.radius:
        parser.error("Propagation radius must exceed the R-matrix sphere radius")
    if args.orbitals == "state-averaged":
        if args.sa_singlet_roots[0] < 1:
            parser.error("State averaging must include the A1 singlet ground state")
        for counts, targets in (
            (args.sa_singlet_roots, args.target_singlet_roots),
            (args.sa_triplet_roots, args.target_triplet_roots),
        ):
            if min(counts) < 0 or counts[1] != counts[2]:
                parser.error("State-average roots must be nonnegative with equal B1/B2 counts")
            if any(count > target for count, target in zip(counts, targets, strict=True)):
                parser.error("UKRmol target roots must cover every state in the orbital ensemble")
    if args.target_initial_checkpoint is not None:
        if args.orbitals != "state-averaged" or not args.target_initial_checkpoint.is_file():
            parser.error("Initial checkpoint requires state-averaged orbitals and an existing file")
    scripts = acquire_scripts(args.scripts_cache)
    workdir = args.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=False)
    (workdir / "scratch").mkdir()
    shutil.copy(scripts / "scripts/main.pl", workdir / "main.pl")
    shutil.copy(Path(__file__).with_name("co.pl"), workdir / "co.pl")
    shutil.copytree(scripts / "lib", workdir / "lib")
    shutil.copytree(scripts / "input.templates", workdir / "templates")
    congen_template = workdir / "templates/congen.inp"
    congen_template.write_text(
        congen_template.read_text().replace(
            "&state",
            f"&state\n  ndimx = {args.congen_workspace},"
            f"\n  cdimx = {args.congen_workspace // 10},"
            f"\n  nodimx = {args.congen_workspace // 10},",
            1,
        )
    )
    for name in ("target.scatci.inp", "scattering.scatci.inp"):
        template = workdir / "templates" / name
        source = template.read_text()
        if source.count("memp = 2.5,") != 1:
            raise ValueError(f"No unique MPI-SCATCI memory setting in {name}")
        template.write_text(source.replace("memp = 2.5,", f"memp = {args.scatci_memory_gib},"))
    if args.target_diagonalizer != "auto":
        select_target_diagonalizer(
            workdir / "templates/target.scatci.inp",
            args.target_diagonalizer,
            args.target_diagonalizer_tolerance,
            args.target_diagonalizer_max_cycles,
        )
    psi4_template = workdir / "templates/psi4.inp"
    psi4_template.write_text(
        "set scf_type pk\nset e_convergence 1e-10\nset d_convergence 1e-10\n"
        + psi4_template.read_text()
    )
    instrument_library(workdir / "lib/ukrmollib.pm")
    if args.orbitals == "state-averaged":
        enable_state_average(workdir / "lib/ukrmollib.pm")
        psi4_template.write_text(
            "# PySCF state-averaged backend: the input is the retained run config.json.\n"
            "# The upstream quantum-chemistry slot/output name is psi4.\n"
        )
    config = vars(args) | {
        "workdir": str(workdir),
        "scf_type": "conventional-exact" if args.orbitals == "state-averaged" else "pk",
        "scripts_url": SCRIPTS_URL,
        "scripts_md5": SCRIPTS_MD5,
    }
    if args.target_initial_checkpoint is not None:
        checkpoint = workdir / "initial.casscf.chk"
        shutil.copy(args.target_initial_checkpoint, checkpoint)
        config["target_initial_source_checkpoint"] = str(args.target_initial_checkpoint)
        config["target_initial_checkpoint"] = str(checkpoint)
    (workdir / "config.json").write_text(json.dumps(config, default=str, indent=2) + "\n")
    env = os.environ | {
        "UKRMOL_SCRIPTS": str(scripts),
        "UKRMOL_RUN_DIR": str(workdir),
        "UKRMOL_STAGE_LOG": str(workdir / "stages.tsv"),
        "UKRMOL_BOND_LENGTH": str(args.bond_length),
        "UKRMOL_MODEL": args.model,
        "UKRMOL_BASIS": args.basis,
        "UKRMOL_FROZEN_ORBITALS": str(args.frozen_orbitals),
        "UKRMOL_ORBITALS": args.orbitals,
        "UKRMOL_TARGET_ONLY": str(int(args.target_only)),
        "UKRMOL_ACTIVE_ORBITALS": ",".join(map(str, args.active_orbitals)),
        "UKRMOL_TARGET_ROOTS": str(args.target_roots),
        "UKRMOL_TARGET_SINGLET_ROOTS": ",".join(map(str, args.target_singlet_roots)),
        "UKRMOL_TARGET_TRIPLET_ROOTS": ",".join(map(str, args.target_triplet_roots)),
        "UKRMOL_TARGET_STATES_USED": str(args.target_states_used),
        "UKRMOL_VIRTUAL_ORBITALS": ",".join(map(str, args.virtual_orbitals)),
        "UKRMOL_RADIUS": str(args.radius),
        "UKRMOL_MAXL": str(args.maxl),
        "UKRMOL_RANKS": str(args.ranks),
        "UKRMOL_PRECISION": args.precision,
        "UKRMOL_DELETION_THRESHOLD": str(args.deletion_threshold),
        "UKRMOL_PROPAGATION_RADIUS": str(args.propagation_radius),
        "UKRMOL_ENERGIES": str(args.energies),
        "UKRMOL_ENERGY_START": str(args.energy_start_ev),
        "UKRMOL_ENERGY_STEP": str(args.energy_step_ev),
        "TMPDIR": str(workdir / "scratch"),
        "PSI_SCRATCH": str(workdir / "scratch"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        # The upstream driver changes into a geometry directory before invoking
        # the backend; child Python must resolve the read-only source snapshot.
        "PYTHONPATH": os.pathsep.join(
            filter(None, [str(Path(__file__).resolve().parents[2]), os.environ.get("PYTHONPATH")])
        ),
    }
    # The entrypoint selects double libraries. Quad executables require the
    # matching GBTOlib ABI, not just a change to the executable directory.
    library_paths = os.environ.get("LD_LIBRARY_PATH", "").split(":")
    library_paths = [
        path
        for path in library_paths
        if path not in ("/opt/ukrmolp/lib.double", "/opt/ukrmolp/lib.quad", "")
    ]
    env["LD_LIBRARY_PATH"] = ":".join([f"/opt/ukrmolp/lib.{args.precision}", *library_paths])
    env["PATH"] = f"/opt/ukrmolp/bin.{args.precision}:" + env["PATH"]
    print(f"CO {args.model} at R={args.bond_length} bohr: {workdir}", flush=True)
    if args.qc_only:
        started = time.monotonic()
        status = run_profiled(
            workdir,
            env,
            [
                "python3",
                "-m",
                "projects.ukrmol_co.target",
                "--config",
                str(workdir / "config.json"),
                "--output",
                str(workdir / "target.out"),
            ],
        )
        (workdir / "stages.tsv").write_text(
            f"target\tpyscf\t\t\t{time.monotonic() - started}\t{status}\n"
        )
        raise SystemExit(status)
    raise SystemExit(run_profiled(workdir, env))


if __name__ == "__main__":
    main()
