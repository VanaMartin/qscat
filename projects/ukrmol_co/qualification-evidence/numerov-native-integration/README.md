# Native residual coupled integration evidence

[Qualification and design](../../NUMEROV_NATIVE_INTEGRATION.md): producer
**538/538 checks**, 50 propagations / 1,737,200 mixed solves; 1,735,263 native
residual evaluations agree with NumPy, with zero decision disagreements.
Integrated marginal speedup is **2.26–2.27×** on the fixed synthetic 298-channel
families; complete N128 speedup is **1.84–1.88×**.
Saved reconstruction and full public numerical replay both pass: all 50 rows,
residual observations and 538 decisions reproduce exactly, original replay
exits 0/0. All 14 final payloads fresh-download and verify.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-native-integration
```

Extract `source.tar.gz` and `verification-support.tar.gz` into `source/`,
`results.tar.gz` into `collected/`, `native-wheel.whl` into `extension/` and the
mpmath wheel into `dependencies/`. Source freezes at `9ed0123`, verification
support at `09a4e83`. The result archive contains three pinned prior reports,
all numerical rows, original owners/commands/exits and sampled resource records.

From `source/`, in the qualified pinned Linux runtime:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_native_integration_replay --report ../collected/result/report.json --priors ../collected/priors --output ../saved-decisions.json
```

Exact derived metric reconstruction uses the pinned runtime's NumPy/BLAS
arithmetic. Other hosts can differ in last-bit K-to-S metrics. For complete
numerical replay, use a fresh output directory and append `--reexecute`:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_native_integration_replay --report ../collected/result/report.json --priors ../collected/priors --output ../fresh-replay --reexecute
```

Numerical rows and residual-audit aggregates must reproduce exactly; fresh
timings/hypotheses are observations. Both nested wheel components and the whole
wheel are pinned in the [replay contract](../../numerov-native-integration-public-replay-contract.json).
This reuses the validated binary; an independent native rebuild is outside the
numerical replay scope. Actual CO conditioning/throughput and native K accuracy
remain separate gates.
