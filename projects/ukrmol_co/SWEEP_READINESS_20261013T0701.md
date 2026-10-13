# CO sweep readiness — 13 October, 07:01 UTC

**B1 completes all native/outer stages and passes independent same-root
reproduction; B2 has completed its matrix classes and L2 update.**
Owner1228004 remains active on Sadaharu at **07:01:25UTC**, after **11.4172hours**
of the two-sector control. Current sampling is **801.16% CPU /9.743GiB memory**
under the declared eight-core/80-GiB allocation. These are ongoing usage samples,
not final cost; no whole-run resource result exists yet.

## Independently verified B1

SCATCI,rsolve,eigenp andreson all exit0 without timeout. Native SCATCI takes
**21416.5942s** and rsolve123.6785s; eigenp/reson are each about0.016s.
The original and candidate B1 native/outer decks match byte-for-byte, and
the native target map is unchanged.

Independent local reading verifies2048finite ordered CI/boundary poles,
complete continuum/boundary coefficients and the full99-point phase grid.
SLEPc requests2048roots, converges2050 with reason1 in75iterations, retains
selection1/mpisbaij/krylovschur and reports maximum absolute action residual
**1.3199200374804417e-10Hartree**. CI and boundary energies match exactly.

Against the original2048-root baseline:

- **Maximum exported-energy difference:0Hartree.**
- **Maximum full-grid phase difference:0rad moduloπ.**
- Channel/threshold/multipole interface digests match.
- Identical-window fits for background orders1–4 all pass, with **0eV
  position and full-width differences**. Native fits are nonempty and notMAXFIT.

The [07:01 capture](progress-engine-baseline-20261013T0701.json) includes the
independent reader source hashes, raw-source checksums, solver coverage,
native original exits and per-background fit diagnostics. Its **B1-only
engine-reproduction verdict is true**. Full two-sector/Pi/engine reproduction
remains unset until B2 and the declared final reader finish.

## B2 and remaining gates

All eight B2 rank logs record Class1,Class3,Class4,Class567,Class2+8 and L2
completion, with4Class5/4Class6/24Class7pair exits per rank. Rank0 prints
the SLEPc/KrylovSchur solve setup for dimension86352 and2048requested roots.
No B2 eigenpair/boundary export or residual telemetry exists at the cutoff.
Buffered tails do not resolve the subsequent live solve operation. No
pointer-failure, intercepted-signal or bounds-error marker appears in either
sector. B1's normal exit preserves IEEE underflow/denormal notices and zero
native warning/error message counts; the generic SLEPc symmetry warning is
also preserved and does not replace the pending Pi check.

The [same-root contract](cas12-audited-engine-baseline-contract.json), exact
audited build/runtime and0.1meV position/full-width /1e-6rad phase reproduction
budgets remain in force. B1 agreement establishes reproduction for that sector,
not a crash correction or root-count/omitted-spectrum convergence.

All **41source/343input-evidence hashes**, **242installed engine artifacts**
and both sectors' five native inputs reverify against the original baseline.
The capture pins107immutable remote records plus its frozen observation.
Observer1228557 remains active with12hourly captures19:36UTC12October through
06:37UTC13October. Source/run/observer outputs remain under
`prepared/cas12-audited-engine-baseline-source-20261012/`; owner and observer
services are active. Stage/whole-native/owner limits remain24/50/72hours.
`STOP` in the fresh source directory ends the active step.

**Next:** complete B2 native/outer stages and the final independent two-sector
reader, including Pi and both sectors' reproduction checks. A passing control
precedes a separately declared scientific4096-root comparison. Earlier crash
causes, engine correction, omitted-spectrum/model/channel/identity/extraction/
neutral and frozen restartable-production-pilot gates remain separate.
