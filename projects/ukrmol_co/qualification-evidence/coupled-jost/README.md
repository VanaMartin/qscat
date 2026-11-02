# Coupled canonical-Jost endpoint control independently passes

[Self-contained construction and gates](../../../../docs/physics/coupled-jost-boundary-control.md).
Pushed `16f76fb` before2November08:59:09UTC launch;54 integrations finish08:59:10UTC,
all six science/control/command/container/owner/process exits0. Six analytic/
synthetic runtime controls pass first, no ODE before freeze.4 physical cores/
2GiB/one library thread/no network; seven minimal source files, no CO input mount.
Persistent capsule `/home/kooza/ukrmol/co-sa-20261006/prepared/coupled-jost-source-20261102`.

Independent reader `e8efa4f` uses the full hyperbolic fundamental seed and F,F'
instead of the integrator's rank-one superpotential formula; scalar one-row
least squares replaces its lstsq call. No ODE import or execution. Eight fresh
public payloads authenticate12 raw files/seven executed Git sources/three reader
Git sources; all54 saved maps/metrics/gates/exits/image/mount/resource caps qualify.

|rtol|Maximum absolute Y/A error, exact starts|
|---|---:|
|1e-8|1.40548e-8|
|1e-10|1.40328e-10|
|1e-12|1.39615e-12|

Finest free-tail errors versus the independent coupled oracle:

|Angle|Length4|Length8|Length12 /bohr|
|---|---:|---:|---:|
|20degrees|2.23166e-5|1.19378e-8|6.21718e-12|
|40degrees|8.75296e-5|1.91535e-7|4.18703e-10|
|50degrees|2.35071e-4|1.37672e-6|8.04356e-9|

All finest exact starts and length12 free starts pass preregistered map1e-7/
flux1e-7/Wronskian1e-6 gates; finite-tail reductions exceed100 at every angle.
Shorter free starts intentionally expose truncation error; they are not claimed
qualified maps. `s_identity_error` in raw rows is distance from free S=1, **not**
an error against the coupled exact S or a gate; coupled scattering is nontrivial.
Science0.576349s/owner1.236477s; cgroup high-water56,229,888bytes (0.052368GiB)
recorded at science completion. No timeout/OOM/retry or changed parameters.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/coupled-jost
uv run --no-sync python -m projects.ukrmol_co.verify_coupled_jost projects/ukrmol_co/qualification-evidence/coupled-jost/evidence.tar.gz projects/ukrmol_co/qualification-evidence/coupled-jost/collection.json
```

**Decision relevance:** canonical wave/derivative normalization and coupled
transport work here; actual endpoint-tail error can be separated from the ODE
floor. Convergence above the diagonal45degree crossover shows that a crossover
alone neither certifies nor rejects actual finite-tail initialization. This
exponential-tail control does not certify the CO inverse-power expansion or fix
the original failed guard. Original six exits2/native-export qualifications and
114 untouched holdouts remain. No active job or automatic CO successor.

A [separate three-map CO study proposal](../../stretched-contour-study-proposal-20261102.json)
is **unreleased pending explicit authorization**; it changes the contour and
therefore cannot be smuggled in as a retry of the frozen40degree guard.
