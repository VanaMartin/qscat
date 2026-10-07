# Fresh-target CI coverage and residuals

All three passing [fresh-basis QC targets](../fresh-basis/README.md) pass an
independent final-orbital coverage scan under the **600-cycle** numerical
contract. Each target has 48 probes: both spins, all four C2v sectors, five/eight
roots and forty/eighty/160-vector trial spaces. Orbitals, tolerances, spin
penalty and ensemble weights stay fixed. The forty orbital-ensemble components
remain distinct from the 64 roots available in eight-root sector probes.

| Fresh target | Eigenpairs checked across refinements | Maximum physical residual (Hartree) | Maximum penalized residual (Hartree) | Maximum trial-space energy difference (Hartree) |
|---|---:|---:|---:|---:|
| CAS(10,11), cc-pVTZ | 312 | 6.66848e-10 | 9.98075e-10 | 2.27374e-13 |
| CAS(10,12), cc-pVTZ | 312 | 6.93428e-10 | 9.99432e-10 | 2.55795e-13 |
| CAS(10,12), aug-cc-pVTZ | 312 | 6.96009e-10 | 9.95671e-10 | 3.12639e-13 |

Every convergence flag, spin, physical/spin-penalized residual and CI-vector
orthogonality gate passes for all **144 probes / 936 eigenpair evaluations**.
Maximum first-five ensemble-energy difference is 2.84e-13 Hartree; changing
five→eight requested roots changes their common first five by at most
3.98e-13 Hartree. Maximum CI-vector orthogonality error is 6.42e-14.

The residual measurement first passes an analytic two-site Hubbard-dimer
energy/residual benchmark and rejects a perturbed eigenvector. The normalized
physical residual uses the unpenalized Hamiltonian and its physical Rayleigh
energy; the penalized residual uses the solved eigenvalue. Their limits both
remain **1e-9 Hartree**. The scan uses PySCF 2.11.0, CI energy tolerance 1e-12
Hartree, `lindep=1e-22`, spin-penalty shift one and four threads. The ordinary
target builder and `ci_probe.py` retain their 200-cycle defaults.

Original exit is **zero**, wall **1941.87 seconds / 32.36 minutes**, kernel
memory peak **781508608 bytes / 0.728 GiB**, and aggregate CPU time
5645.22 seconds. The CPU-4–7 lease acknowledges its waiting owner before
launch and resumes it afterward. The queued fresh aug-TZ repair then acquires
the slot independently.

Raw spectra/spins and all 144 whole-probe convergence blocks reconstruct from
`run.log`. Individual `casci.log` files are empty because PySCF retains its
own stdout handle; the aggregate raw log is the numerical evidence. The
verifier recomputes spectral, spin, orthogonality and cross-space/root-count
gates and checks the recorded Hamiltonian-action residuals. The frozen source
and analytic controls document how those residuals were measured.

The archive includes unchanged parent QC checkpoint inputs, target/configuration
records, all probe records, raw logs, resources, frozen source, lease/resumption
and the executable verifier. Parent checkpoints match the checksum-pinned
fresh-basis companion. All **352 payloads** verify, with byte-identical
repackaging and a byte-for-byte public fetch. Competing orbital starts and
independent UKRmol all-root/dipole imports remain required for qualification.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fresh-coverage
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/fresh-coverage/fresh-basis-coverage.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-fresh-coverage.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

To recompute the scans on Sadaharu with the recorded image installed, use a
fresh persistent run root. Set `CPUS` to four available physical cores:

```bash
D="$EVIDENCE_DIR/state-averaged-evidence/diagnostics/fresh-basis-ci-residual-coverage"
RUN_ROOT=$(mktemp -d "$HOME/ukrmol/co-fresh-coverage-recompute.XXXXXX")
mkdir -p "$RUN_ROOT/runs" "$RUN_ROOT/diagnostics/fresh-basis-ci-residual-coverage"
cp -R "$D/inputs/." "$RUN_ROOT/runs/"
docker run --rm --cpuset-cpus="$CPUS" --memory=8g --memory-swap=8g \
  -v "$RUN_ROOT:/work" -v "$D/source:/opt/qmodeling:ro" -w /opt/qmodeling \
  -e PYTHONPATH=/opt/qmodeling -e OMP_NUM_THREADS=1 \
  -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  --entrypoint /opt/ukrmolp/entrypoint.sh \
  sha256:6b3b0ededa85494076229b111efe969985407023e404bff888f17cafb4630666 \
  python3 /opt/qmodeling/fresh-basis-residual-coverage.py --child
```

The [continuation](../../CONTINUATION.md) records the released
[fresh CAS(10,11) TZ import](../../calibration-sa11-fresh-tz-import.json)
and the remaining basis/start and continuum checks.
