# Saved aug-TZ thresholds block the unchanged inner/export prescription

The [qualified own preparation](../stretched-augtz-preparation/README.md) remains
successful. Before allocating the hours-long B1/2048 solve, a saved-only screen
reads its authenticated CIDATA and the unchanged RSOLVE deck. Five freshly fetched
public payloads and seven Git source pins reproduce the result; seven no-native
tests cover closed/crossing/endpoint and changed-grid/invalid-input cases.

| Quantity | Hartree | Physical eV |
|---|---:|---:|
| triplet.B1 root1 threshold |0.18159458918989912|4.941440506627862|
| triplet.B2 root1 threshold |0.18159458919126337|4.941440506664985|
| Unchanged native grid maximum |0.18375|5.000092222700295|

The requested grid is `nescat=99`, `einc=0.1,0.05`, `ieunit=2`. The pinned native
conversion is0.0735Rydberg per requested eV; modern reported eV uses `qscat.units`.
Both thresholds lie below grid indices **97/98** (zero-based), requested
**4.95/5.00eV**. `native_full_reconstruction._model_from_interface` explicitly
requires99×6×6 saved K matrices and every remaining channel threshold strictly
above the full grid. The prerequisite therefore **rejects**. Screen execution
exit0 means the check ran, not that the prerequisite passed. No precise new open-
channel count is claimed without a native interface export.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-augtz-inner-preflight
uv run --no-sync python -m projects.ukrmol_co.stretched_augtz_inner_preflight projects/ukrmol_co/qualification-evidence/stretched-augtz-inner-preflight/fort.26 --deck projects/ukrmol_co/qualification-evidence/stretched-augtz-inner-preflight/rsolve.inp
```

The compact bundle supplies identical CIDATA/deck bytes, source and input hashes,
and result. CIDATA SHA256 agrees with the full public preparation collection.
An initial generic Python HTTP fetch received403; the established artifact-fetch
User-Agent retrieved all five authenticated payloads. No science rerun occurred.

**Stop the unchanged recipe before expensive execution.** No grid narrowing,
relaxed gate, target retry, native/CI/Hamiltonian/map call, aug-TZ inner contract
or allocation was made. This is an export/reference-prescription blocker, **not**
a rejection of aug-TZ electronic validity or its qualified properties/preparation.
Next requires an explicit threshold-aware export/reference design and qualification
scope before a separately frozen costly solve. No automatic B2/geometry/pole/
fit/sweep successor; all114 S holdouts and original failures stay untouched.
