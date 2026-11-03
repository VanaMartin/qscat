# Matched5Z atoms complete;finite-basis insensitivity rejects

All four aug-cc-pV5Z atoms finish **3November12:09:56UTC** from source`e5f8e74`.
The actual-image preflight and all eight runtime controls pass. Same calculation
function bytes as QZ,same states/occupations/core/ROHF-reference UCCSD(T),unchanged
numerical/state/stability/memory gates. All original exits0,no timeout/OOM/retry.
The earlier [zero-atom startup stop](../atomic-5z-launch-stop/README.md) remains
separate;the corrected source-byte control does not requalify that failed attempt.

| Species |5Z total / Hartree|
|---|---:|
| C(3P)|−37.78825666178488|
| O(3P)|−75.00059248661933|
| O-(2P)|−75.0530193069947|
| C-(4S)|−37.83386256840955|

The5Z fragment sums are C+O **−112.78884914840421Ha**,C+O-
**−112.84127596877957Ha**,C-+O **−112.83445505502888Ha**.

| Diagnostic |QZ|5Z|5Z−QZ|Predeclared finite-basis gate|
|---|---:|---:|---:|---|
| O EA /eV|1.403072709|1.426606459|+23.533750meV|10meV:reject|
| C EA /eV|1.235049235|1.240999940|+5.950705meV|10meV:pass|
| C+O offset from own molecularR=2.1323 reference /eV|11.115011294|11.167927472|+52.916179meV|20meV:reject|
| C+O- offset,same reference convention /eV|9.711938585|9.741321013|+29.382428meV|20meV:reject|
| C-+O offset,same reference convention /eV|9.879962059|9.926927532|+46.965473meV|20meV:reject|

**These offsets are not De,D0 or initial-v DA thresholds.** The molecular zero
is fixedR=2.1323,not an independently optimized minimum or vibrational level.
The independently qualified well-shape QZ/5Z diagnostic remains distinct:small
relative/interpolation errors do not establish a stable atomic alignment.

Eleven fresh public payloads reconstruct48 full raw files,21 executed Git sources,
12 source prerequisites,573 resource samples and the original QZ bundle. Raw
SCF/CC/triples logs give the totals;QZ and5Z sums/EAs and matching-basis molecular
reference offsets are independently recomputed. State/spin/occupation gates are
source-authenticated diagnostics,not fresh wavefunction/integral calculations;
full correlated L2 purity and absolute atomic accuracy remain unproved.

Science114.268448s,owner **114.673835s**,CPU302.461764s;container kernel cumulative
high-water **12,760,150,016bytes =11.883816GiB**,below unchanged16GiB allocation.
Four physical cores0–3;conditional1–5minute forecast contains actual cost,no
speedup claim. Immediate terminal observation precedes first600s periodic check;
no campaign quantum job/monitor remains active.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints
uv run --no-sync python -m projects.ukrmol_co.verify_atomic_5z projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/evidence.tar.gz projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/collection.json projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/evidence.tar.gz projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/collection.json
```

**Decision:** stop endpoint basis-insensitivity/production-readiness inference.
Both finite datasets remain useful;two bases do not give a CBS error bound.
No6Z,automatic extrapolation or gate relaxation. Next one bounded neutral/
correlation/reference-treatment decision that addresses the missing molecular
dissociation connection and reference error,rather than another basis ladder.
Augmented scattering stays paused;all114 holdouts and prior failures stay intact.
