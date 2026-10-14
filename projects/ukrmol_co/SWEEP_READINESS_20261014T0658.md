# CO sweep readiness — 14 October, 06:58 UTC

**The scientific 2048→4096-root pair-stability comparison passes in both sectors.**
Owner 1246007 finishes at **13 October 22:55:24 UTC** with all native/outer
stages, independent reader and owner exits 0. The observation-only observer ends
normally after 16 immutable captures. No additional native calculation is running.

## Measured result

The [completed comparison](progress-audited-4096-20261014T0658.json) independently
checks both sectors against the qualified audited2048baseline, on the same exact
engine/direct-MPI/residual-telemetry runtime:

| Check | Measured result | Declared budget |
|---|---:|---:|
| Shared 2048-prefix energies | 0 Hartree difference | 1e-7 Hartree |
| Full-grid phase difference | 1.67e-5 rad moduloπ | 0.05 rad |
| Maximum fitted position difference | 2.2383495e-6 eV (0.00224 meV) | 0.05 eV |
| Maximum fitted full-width difference | 4.1943029e-5 eV (0.04195 meV) | 5%, with 0.001-eV floor |
| Maximum native action residual | 1.3276448e-10 Hartree | 1e-7 Hartree |
| Pi phase difference | 1.0000001e-9 rad moduloπ | 0.002 rad |

Position/full-width maxima cover background orders 1–4 over the same requested-
energy window. Every fit passes with successful optimization and no parameter
bound; native fits are nonempty and not MAXFIT. Inputs differ only in `nstat`;
actual Hamiltonian data, outer decks, target map, selection and channel/threshold/
multipole interfaces match. Both sectors export 4096 finite ordered CI/boundary
poles, 4096×1410 continuum coefficients and 4096×250 boundary amplitudes.
CI/boundary energies agree exactly. The emitted range extends to approximately
−108.9627014 Hartree. SLEPc converges 4104/4102 roots in 159/160 iterations,
reason 1, retaining selection 1; the requested/exported count is 4096.

**Native cost: 14.924481 hours / 21.736774 GiB kernel peak**, including both
matrix builds, solves and outer stages, excluding target construction and the
separate independent reader. The matched2048control cost was 11.565930 hours /
21.964687 GiB. These are measured run costs, not a general scaling law.

## Reconstruction and evidence

[Local verification](verification-audited-4096-20261014.json) re-reads 40 checksum-
verified fetched files: binary poles/vectors/amplitudes, residual telemetry, grids,
decks and fits for both root counts and sectors. It independently passes the
original scientific budgets. Full-grid and shared-energy differences match the
remote reading; cross-runtime optimizer fit-difference drift is at most
2.713e-10 eV. The record preserves an extra local parity assertion at 1e-10 eV
that rejected this tiny drift on arm64/NumPy2.5.1/SciPy1.18. That assertion was
supplementary to the declared scientific budgets; no native retry or budget change
was made, and measured drift is reported without a replacement parity tolerance.

The completion capture reverifies **43 source / 500 input-evidence hashes**,
**242 installed artifacts**, both sectors' actual native/outer inputs, final
original exits and raw output hashes. It pins **665 immutable remote records**.
All persistent artifacts remain under
`prepared/cas12-audited-engine-4096-source-20261013/`. Original release/debug/
instrumentation failures remain preserved. Completion establishes neither their
cause nor an engine correction.

## Next scientific decisions

The pair shows very small sensitivity between these two resolutions, supporting
reuse of the tested2048baseline for matched controls. **Absolute omitted-spectrum
error remains unqualified.** Any third resolution or independent omission test
needs a separate declaration. Public raw reconstruction remains pending.

Then complete explicit-A2/channel and matched three-anchor basis/active-space/
ensemble controls, continuity and extraction qualification, neutral correlation/
domain uncertainty and a frozen restartable production pilot. Model and pole
identity remain unresolved; the full sweep and held-out potential fitting remain
gated. Completion schedules no additional native calculation.
