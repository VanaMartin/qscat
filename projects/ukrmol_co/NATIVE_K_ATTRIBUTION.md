# Stretched 1-eV numerical-error attribution

**Status: original attempt stopped on its native32 timeouts; verdict unset,
original runner/container exits 1/1.** Both native processes were killed with
exit −9 at the frozen 180-s cap. No K file was produced, no native64 ran and no
independent control was submitted. The
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

The saved step pair costs 59.75–61.59 / 114.72–115.73 s per sector. New work
forecasts **600 s wall time**, upper1200 s, hard cap1800 s, **four single-threaded
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

## Resource successor running — 18 October, 22:37 UTC

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
notification is attached. On completion, collect outputs, independently verify
all decisions and publish/replay the bounded evidence before qualification.
