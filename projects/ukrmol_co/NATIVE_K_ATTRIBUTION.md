# Stretched 1-eV numerical-error attribution

**Status: the resource-corrected experiment is complete and rejected, original
runner/container exits 1/1.** Four native invocations exit0 and all ten independent
invariant checks pass, but both sectors fail the separate step, endpoint and
native-reference K gates. The original180-s attempt remains unset/exits1/1. The
[original six-case rejection](NATIVE_FULL_RECONSTRUCTION.md) is independently
reproduced and remains false, exits 1/1. The large-tangent-coordinate diagnosis
explains amplification but leaves the underlying numerical error unresolved.

The [design](../../docs/physics/co-native-k-attribution.md) and
[contract](native-k-attribution-contract.json) restrict the successor to the
same stretched B1/B2 1-eV identities. It reuses four saved independent solves
and the native16/32 matrices, then performs:

1. Two native32 single-energy reproductions, requiring identical energy bits,
   complete records and K/S/phase differences ≤1e-10 before refinement.
2. Two native64 single-energy controls, requiring separate native32→64
   K/S/phase differences ≤1e-6 before any independent control.
3. Ten independent controls: third step 0.00625 bohr, endpoints 200/800 bohr,
   and tail orders 4/6 in each sector. Existing step0.025/0.0125, endpoint400,
   order8 rows complete three-rung one-factor families.

The fixed family gates and native comparison retain separate **1e-6 K/S/phase**
vetoes and unsymmetrized invariant checks. A false or unset gate stops without
retry or automatic full-grid continuation. All 298 channels, 2048 saved poles,
the electronic frame, inner boundary and validated Numerov method stay fixed.

The saved step pair costs 59.75–61.59 / 114.72–115.73 s per sector. The original,
incorrect forecast was **600 s wall time**, upper1200 s, hard cap1800 s, **four single-threaded
physical cores0–3**, **4 GiB**, four native invocations and ten independent
solves. No new inner solve or native complex-pole claim is involved.

**29 focused tests pass**, including exact energy serialization, B1/B2 output
unit preservation, separate K vetoes, missing-control rejection, failed-slice
stopping and preservation of zero native exits when a scientific gate rejects.

```bash
python -m projects.ukrmol_co.native_k_attribution \
  --raw /native --prior /prior \
  --engine /opt/ukrmolp/bin.double/rsolve \
  --prerequisite /source/projects/ukrmol_co/verification-native-full-public-20261017.json \
  --contract /source/projects/ukrmol_co/native-k-attribution-contract.json \
  --output /outputs/result
```

The matching source, input, saved-reference and prerequisite hashes must be
verified before execution. Full native K reconstruction, CO poles and production
remain gated regardless of this two-identity result.

## Preserved timeout and cost correction

Source/design/contract were frozen and pushed at `94fccee` before execution.
The original run lasted **180.617567 s**, 22:27:40–22:30:40 UTC on 18 October.
Sampled container CPU was **360.070375 s**, recorded kernel memory high-water
**904,531,968 bytes / 0.84241 GiB**, with 900 samples at 0.2 s. All242 installed
engine artifacts, eight inputs, eight saved references and nineteen source files
verify before/after execution. Partial logs are 16,508 bytes in each sector;
buffering prevents locating the exact executing routine at termination.

The original native forecast of about twenty CPU-seconds per invocation was
wrong. Checksum-pinned `rsolve.f` calls `ASYM1` before opening K output and before
the scattering-energy loop. `RPROP1` in `rpropg.f` constructs and diagonalizes
one **298 × 10 = 2980**-dimensional Hamiltonian per radial subrange, independently
of NESCAT. Saved stretched native32 full-grid invocations cost **866.028 /
920.267 s**, rather than twenty seconds. A one-energy deck removes per-energy
work but retains this setup. This explains why the short cap was unjustified;
the logs do not prove which routine was executing at timeout.

[Public evidence](qualification-evidence/native-k-attribution-timeout/README.md)
preserves the original archive and exits. Fresh public-byte verification passes
**12 payloads / 19 source files / 2 engine sources**, the separate saved-decision
replay and all29 focused tests. It verifies the unset decision and lack of K
output; it does not reexecute the timeout or assert deterministic timeout walls.

A separately frozen [resource successor](../../docs/physics/co-native-k-attribution-resources.md)
uses the same matrices, equations, deck changes, gates, stopping rules and
producer implementation with an explicitly corrected setup-cost allowance.
The original contract, timeout, unset result and exits remain separate evidence.

## Resource-successor launch — 18 October, 22:37 UTC

The resource successor was frozen and pushed at **`4a558b3`** before launch at
**22:37:20.572451 UTC**. The [launch observation](progress-native-k-resource-launch-20261019.json)
confirms both native32 slices executing at nearly100% of one CPU each, with
native64/independent stages unsubmitted and all verdicts unset. Source archive
SHA256 is`a376ca82be4e3dd8a9f2a6954e603f2b4c079aebf759b24639b126e42fb64f66`.
All23 source files,242 installed artifacts,eight inputs/eight saved references,
the published timeout archives and the exact scientific-contract comparison pass
preflight. The resource-only contract changes four timing allowances.

Allocation: four physical cores0–3,4 GiB,four possible native calls and ten
independent controls; forecast3400 s,upper5400 s,batch cap9000 s,native-stage
cap3600 s. The observation's718,954,496-byte recorded high-water and28.343329
CPU-seconds are early samples, not final costs. Sadaharu has61.76 GiB free
persistent space. No inner solve runs.

Persistent directory:
`/home/kooza/ukrmol/co-sa-20261006/prepared/native-k-resource-successor-source-20261019/`.
Owner74564 records0.2-s resources,original exits and post-run hashes. A process-exit
notification delivered the completed result. The original launch observation
remains an explicitly dated running-state record.

## Completed resource successor — 18 October, 23:05 UTC

The fixed study finishes at **23:05:10.231254 UTC** with **false / original exits
1/1**. Both native32 one-energy slices reproduce the saved full-grid record
exactly, including energy bits and all K/S/phase components. Native32→64 passes
the separate1e-6 gates; all four native subprocesses exit0. All ten additional
independent controls finish and pass their unsymmetrized invariant gates.

The [separate saved-decision replay](verification-native-k-resource-20261019.json)
reconstructs four native records, ten new and four reused independent rows and
all thirty case checks. **24/30 checks pass; six K vetoes reject.**

| Final comparison | B1 max entrywise K change | B2 max entrywise K change | Budget / verdict |
|---|---:|---:|---|
| Native32→64 | 1.30702745e-7 | 4.12325107e-9 | 1e-6 / pass |
| Step0.0125→0.00625, endpoint400, order8 | 6.42620334e-5 | 6.64505339e-5 | 1e-6 / reject |
| Endpoint400→800, step0.0125, order8 | 7.50181019e-6 | 7.47786424e-6 | 1e-6 / reject |
| Tail order6→8, step0.0125, endpoint400 | 4.32009983e-11 | 4.37125891e-11 | 1e-6 / pass |
| Step0.00625 /endpoint400 /order8 versus native64 | 7.25612354e-5 | 7.25604985e-5 | 1e-6 / reject |

Every corresponding S/phase comparison and every invariant gate passes.
Finest/native64 S differences are4.83250e-9/5.15136e-9 and determinant-phase
differences7.95845e-9/8.15810e-9 rad moduloπ. These passes do not waive K.

### What the fixed controls establish

- **No observed native subrange failure at these points:**32→64 changes are
  below budget and much smaller than the finest independent/native discrepancy.
  The native16→32→64 K sequence is itself nonmonotone in B1, so this is a
  conditional resolution check, not an absolute native-reference error bound.
- **The independent step sequence is not in a resolved fourth-order regime:**
  the previous changes1.05931e-5/1.22856e-5 grow to6.42620e-5/6.64505e-5 after
  halving the step. The S changes also grow, but their formal trend checks pass
  only because every change lies below the frozen1e-8 floor. That gate is not
  evidence of fourth-order decay here. Richardson extrapolation is unqualified.
- **Tail order is not the dominant measured sensitivity:**4→6 changes are
  about2.5816e-8 K and6→8 about4.4e-11. This is finite-order sensitivity, not
  an absolute omitted-tail bound. Endpoint changes remain much larger despite
  fixed step/order; endpoint variation includes propagation accumulation as
  well as matching/tail effects and does not isolate one by itself.
- **Choosing a favourable row would hide the failure:** at endpoint800 with
  step0.0125, native64 K differences are7.97392e-7 (B1) and1.36790e-6 (B2).
  This does not replace the frozen finest-step/native comparison or the failed
  endpoint family. No best-row selection is performed.
- **Local algebra does not certify propagated accuracy:** equilibrated matching
  condition estimates are about174 at400 bohr,204 at800 and804 at200; local
  residuals are below8.2e-17. Large-channel projections continue to explain almost
  all K discrepancy. Neither those residuals nor the Cayley identity bounds
  accumulated propagation/reference error.

Roundoff/cancellation, endpoint matching and native-reference error remain
possible contributors. The measurements reject a simple monotone spacing cure
and distinguish small observed tail-order changes from larger step/endpoint
sensitivities. They do not uniquely attribute the remaining error.

### Resources, preservation and public replay

Measured wall time **1669.647451 s /27.8275 min**, sampled container CPU
**4005.621716 s**, recorded kernel memory high-water **2,119,385,088 bytes /
1.973831 GiB**,8338 samples at0.2 s. Four physical cores0–3,four single-threaded
workers,4-GiB limit. Native32 calls take428.44–428.89 s; native64 calls take
848.52–849.42 s. The corrected3400-s forecast is conservative for this execution.
All242 installed artifacts /eight inputs /eight prior files /23 sources verify
before and afterward. The39-file result collection preserves original outputs.

The owner's post-verification Docker cleanup failed because its loop variable
overwrote the container name; the original traceback is retained in`owner.log`.
This occurred after final scientific exits and post-run hashes were recorded.
It neither changes the false verdict nor discards the completed container.
The original180-s attempt has the same separately preserved cleanup error.

[Public evidence](qualification-evidence/native-k-attribution-resource/README.md)
publishes the fixed inputs, references, producer, four native outputs and all
ten controls. Fresh downloads verify all39 collected files and23 source files,
binary/report agreement, disk/report row identity, the independently expressed
false decision and all29 focused tests. The final public bundle has **15
payloads**. Numerical public replay of exactly the ten independent controls
reproduces **every stored K and S component exactly**, with identical case
checks and false verdicts. It uses the same frozen producer; native matrices
are checked from their saved complete records. Replay exits **0/0** remain
separate from the original scientific exits **1/1**.

Replay finishes at23:15:08.917600 UTC in **392.657157 s**,1449.435336 sampled
CPU-seconds,231,624,704-byte recorded kernel high-water,1959 samples at0.2 s,
four physical cores/4 GiB. Its17-file archive preserves all ten numerical rows,
the worker/owner, reports and resources. The
[23:15 UTC result observation](progress-native-k-resource-result-20261019.json)
reverifies **56 collected hashes /23 sources /eight inputs /eight prior files /
242 installed artifacts**, with no running scientific container and61.71 GiB
free persistent disk.

No further numerical control is authorized automatically. A successor would
need a fixed-step precision/algebra control or a genuinely independent
propagation formulation, validated on analytic near-large-K cases before native
use. In particular, the recurrence's`12 B⁻¹ − 10 I` can be rewritten exactly as
`2 I + h² B⁻¹ W` for`B = I − h² W/12`; whether that changes accumulated floating-
point error is an untested hypothesis, not an explanation already established.
The two identities do not increase the original132/594 unique-energy coverage;
462 energies, full native K, CO poles and production remain unqualified.

## Subsequent fixed-grid algebra control — 19 October

The [increment-form control](NUMEROV_INCREMENT_CONTROL.md) now tests the
previously untested cancellation hypothesis after analytic discrete-oracle
qualification. Its two-point gates pass24/24 at unchanged physical inputs,
grids and1e-6 budgets. The original two-point rejection documented here remains
false/exits1/1. The candidate's endpoint800/native diagnostic still exceeds
budget, and the result is not an all-row, absolute-error or full-grid claim.
