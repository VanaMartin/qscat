# Fixed-grid Numerov algebra control

**Status: analytic qualification passes; ten-solve CO contract specified.**
The [mathematics and prospective gates](../../docs/physics/co-numerov-increment-control.md)
derive an identity-defect recurrence that avoids cancellation in initialisation,
propagation and matching. `exterior_numerov.py` remains the original oracle.

The [analytic verification](verification-numerov-increment-analytic-20261019.json)
passes16 frozen checks and49 focused tests. At the finest scalar grid, original
versus high-precision discrete K drift is1.34450e-4; the increment form gives
9.82188e-10. The separate continuous discretisation error is6.11308e-8.
The rotated two-open/one-closed case and short noncommuting coupled controls
also pass their defined criteria. This validates the tested analytic regime,
not native CO propagation or a complex pole.

Analytic source freezes: `381a201`, dependency-only supplement `e10e387`.
The first attempt exits1/1 before calculation because the tools image lacks
mpmath; its original logs and archive are preserved. The separately pinned
wheel enables the same numerical source in a fresh directory. The passing run
costs10.4293 s /9.826654 sampled CPU-seconds /46,272,512-byte recorded kernel
high-water,50 samples at0.2 s;30 sources and242 engine artifacts verify
before/after. No native data or engine invocation occurs.

The [fixed-grid CO contract](numerov-increment-co-contract.json) permits only
the existing two stretched1-eV identities at the existing three spacings and
three endpoints: ten candidate solves,368160 intervals,500-s forecast/800-s
upper/1200-s cap,four physical cores/4 GiB. All final K/S/phase and invariant
budgets remain1e-6. Stop submissions on false/unset invariants; in-flight
bundles finish. No new energy, inner/native solve or automatic successor.

The [original two-point rejection](NATIVE_K_ATTRIBUTION.md) remains false,
exits1/1, as does the six-case reconstruction.132/594 unique-energy coverage,
462 uncomputed energies and all pole/production gates remain unchanged.
