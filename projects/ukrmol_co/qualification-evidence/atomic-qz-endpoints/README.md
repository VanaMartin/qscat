# Four calculated atomic endpoints, independently reconstructed

The first aug-QZ endpoint pilot finishes **3November10:53:10UTC**. All four
occupation-constrained high-spin ROHF-reference UCCSD(T) cases complete, one
frozen1s spatial orbital per atom. Original four profile exits,16 stage exits
and all command/container/owner/service/launcher/control exits are0; no timeout,
OOM, orbital retry or peer fallback. Source/inputs remain unchanged.

| Target term | Total energy / Hartree |
|---|---:|
| C(3P) |−37.78675917689432|
| O(3P) |−74.99511693614352|
| O-(2P) |−75.0466789071488|
| C-(4S) |−37.83214639914157|

| Free-fragment channel | Sum / Hartree | Relative to C+O / physical eV |
|---|---:|---:|
| C+O |−112.78187611303784|0|
| C+O- |−112.83343808404312|−1.4030727086292332|
| C-+O |−112.8272633352851|−1.2350492352024132|

These fill the previously missing **finite-method calculated** atomic sums.
They do not establish absolute EA accuracy, a CBS limit or molecular dissociation.
Relative to rounded Dora2016 p4 Table3 external bindings1.461/1.262eV, the computed
O/C EAs differ by approximately−57.93/−26.95meV. This is diagnostic context,
not a precise pass/fail comparison: theoretical nonrelativistic spin-free and
experimental lowest-J conventions/corrections have not been fully reconciled.
The ordering is correct, but no absolute-accuracy pass follows.

Raw SCF/UCCSD/triples log arithmetic and original convergence/stability messages
independently reproduce all totals and sums. Nine fresh public payloads replay
47 raw files,14 original Git sources,nine source prerequisites and91 resource
samples. Full hashes include saved ROHF/UCCSD checkpoints. State/spin/occupation/
angular gates are source-authenticated producer diagnostics, not newly regenerated
wavefunctions/integrals. High-spin valence character is checked; full correlated
L2 purity is not proved. Unsupported external ROHF stability remains null.

Science17.478060s,owner **17.913795s**; summed profiled CPU45.812191s. Container
kernel high-water mark **2,086,653,952bytes =1.943348GiB**; this is cumulative
across the four atoms, not an independent peak per stage. Four physical cores0–3/
16GiB allocation; no parallel-speedup claim. Durable owner scheduled600s status
checks and immediate completion checks. It finished before the first periodic
check; the single terminal observation is17.906786s after owner start, within the
user's900s maximum requested interval. No running campaign quantum job/monitor.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints
uv run --no-sync python -m projects.ukrmol_co.verify_atomic_references projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/evidence.tar.gz projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/collection.json
```

[Exact release](../../atomic-asymptote-release-20261103.json),
[monitoring policy](../../atomic-pilot-release-20261103.md),
[compact result](../../progress-atomic-qz-endpoints-result-20261103.json) and
[physics conventions](../../../../docs/physics/co-atomic-energy-references.md)
remain tracked. The original input ledger intentionally retains unset result
fields: it is an immutable pre-execution input, not the terminal result.

Next one bounded basis/reference-convention decision: determine whether one
matched5Z endpoint control can distinguish basis effects from the remaining
correlation/term-reference uncertainty and connect the existing5Z neutral zero.
No automatic second basis, larger ladder, molecular continuation, augmented
scattering or production release. All114 holdouts and prior failures are intact.
