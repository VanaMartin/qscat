# Stretched 1-eV numerical-error attribution

**Status: prospective fixed experiment.** The
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
