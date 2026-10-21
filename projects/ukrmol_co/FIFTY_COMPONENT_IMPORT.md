# Equilibrium fifty-component target: native pass, seed-frame gate rejects

**21 October 2026:** the target calculation and standard analyzer complete with
original exits **0/0**. The independent verifier exits **1/1**, owner **1**,
because reoptimization changed six of the fifty averaged-root energies beyond
the declared **1e-7-Hartree seed-preservation gate**. The overall import verdict
remains unset. There is no automatic retry, tolerance change or scattering solve.

The [recipe](calibration-sa11-fifty-import-20261021.json),
[contract](fifty-import-contract-20261021.json) and verifier freeze `7827867`
before execution. The model is equilibrium R=2.1323bohr /cc-pVDZ /CAS(10,11),
two frozen A1 spatial orbitals, fifty equal-weight orbital components
`[10,5,5,5]` per spin. Eight native sectors compute `[13,8,8,8]` roots per spin,
**74 total**; DENPROP uses50 states. This is target input preparation, not a
50-versus40/48-channel scattering comparison.

## Two distinct comparisons

1. **Native import against the newly optimized QC frame passes:** all50 averaged
   roots agree within **5.29468e-10 Hartree**; ground dipole error is
   **−1.12928e-11 a.u.** All requested native sector solves finish, as do their
   original stage commands. Final QC gradient **7.43813e-8** passes1e-7, with
   ordinary200-cycle CI and all original QC/fresh-CI flags passing.
2. **New QC frame against the previously covered seed rejects:** maximum root
   change **3.38450e-7 Hartree** is3.3845 times the1e-7 gate. It is approximately
   **0.00921meV**, small compared with the physical50-meV position budget, but it
   cannot be used to waive an independently declared frame-preservation gate.

| Averaged root | New minus seed energy / Hartree |
|---|---:|
| singlet.A1.5 |−3.3844958125e-7|
| singlet.A1.7 |−1.2644662206e-7|
| singlet.A1.8 |−1.0982691379e-7|
| triplet.A1.3 |−2.4702406165e-7|
| triplet.A1.7 |−1.3971720136e-7|
| triplet.A1.10 |−3.1770626663e-7|

The ensemble objective changes by only **−7.10543e-14 Hartree**. Ground-energy
change3.19420e-8 Hartree and dipole change1.28731e-7 a.u. pass their respective
budgets. An unchanged ensemble mean does not guarantee unchanged individual
roots. This is a genuine failed comparison, not the tuple/list schema bug from
the earlier radial verifier.

The original independent verifier stops before its **74-root binary/vector/EOF
and core/active-subspace checks**. Those checks are not claimed as passing.
The original48-probe fixed-orbital coverage belongs to the old checkpoint;
qualification of the new frame remains separate.

## Measured cost and preserved evidence

The target finishes **13:45:21UTC**, before its three-hour hard stop, after
**10436.23s /2.899h**. Sampled CPU43709.40s /12.14CPUh; kernel memory peak
**5.518GiB**. The owner records52,122 resource samples. Physical cores4–11,
eight native MPI ranks, single-threaded native libraries and eight-threaded QC
use the unchanged24GiB container cap.

| Stage | Measured wall |
|---|---:|
| QC/reoptimization and fresh CI |8101.37s /135.02min|
| Eight native target eigensolvers |199.79s /3.33min|
| DENPROP |2131.91s /35.53min|

The one-hour forecast underestimated QC and property work. The existing
checkpoint did not eliminate reoptimization cost; extra MPI ranks do not
parallelize all stages. These timings are not an isolated rank-scaling benchmark.

All **12 frozen sources and7 original inputs** verify before/after. The
[public collection](qualification-evidence/fifty-import/README.md) contains
**153 original files /292,424,932 uncompressed bytes**, including the entire run,
binary outputs, checkpoints, templates, logs, resource samples, owner/source
archive, seed and coverage reference. Its SHA256 is
`dc055c3964aa52474715c70ea405185c7f75b8bd5cd1907bdc8e0577010dd389`.
All four public payloads fresh-download-verify. Saved-table reconstruction
independently reproduces all50 averaged native comparisons and all six
seed-frame rejections; it does not qualify the unvisited74-root binary gate.

Original directories remain under `/home/kooza/ukrmol/co-sa-20261006/`:
`prepared/fifty-import-source-20261021/` and
`runs/co-eq-ccdz-sa11-fifty-slepc-target-mpi8-20261021/`.
Only unrelated ComfyUI remains running at13:47UTC. No scattering inner solve
or repeated target calculation was submitted.

## Next decision

Do **not** repeat three hours of QC just to force agreement or relabel the seed
gate as a pass. The least-wasteful candidate is a separately frozen fixed-orbital
coverage calculation on the **completed new checkpoint**, then native all-root
reconstruction against that checkpoint. Previous two-frame96-probe coverage cost
1280.58s /.658GiB; a one-frame successor can be forecast at roughly10–20min,
but needs its own explicit reference, scope and stop rule. It must preserve the
original rejection rather than inherit old-frame coverage or relax tolerances.
An exact-old-frame import is another option, requiring its own input route.
Neither successor is submitted here.

After a qualified frame, determine physical channel/subspace matching and cost
the matched orbital-ensemble scattering contrast. The [production decisions](PRODUCTION_DECISIONS.md)
still require independent basis/active-space/channel uncertainty, reliable
width extraction, continuity, neutral correlation/domain and a restartable pilot.
