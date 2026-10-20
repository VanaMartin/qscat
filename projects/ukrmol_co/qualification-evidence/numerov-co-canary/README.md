# Saved-CO recurrence canary evidence

[Design and result](../../NUMEROV_CO_CANARY.md): **303/303 checks**, 36 prefixes /
3,996 mixed solves across six saved CO cases. All 180 audited residual matrices
match NumPy; zero decision disagreements. Early-region marginal speedup is
**2.23–2.24×**, before 18.2 bohr. No tail matching or native K qualification.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-co-canary
```

Extract `source.tar.gz` and `verification-support.tar.gz` into `source/`,
`results.tar.gz` into `collected/`, `raw.tar.gz` into `raw/`, `native-wheel.whl`
into `extension/`, and the mpmath wheel into `dependencies/`. Producer source
freezes at `8520bbb`, separate verification support at `646257e`. Verify all
source-manifest/support hashes and all 24 raw input hashes before execution.
The source archive includes the original input contract and selection.

From `source/`, on the qualified x86-64 Linux host with the pinned runtime:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_co_canary_replay --report ../collected/result/report.json --operands ../collected/result/operands --output ../saved-verification.json
```

This reconstructs all saved decisions, checks twelve full captured operand
triples against actual solve hashes, recomputes their native/NumPy residuals,
and recalculates 36 paired 80/110-digit residual references. Adding `--raw ../raw`
reexecutes exactly six audited native N8 prefixes / 90 mixed solves; every final
delta, actual operand hash, condition estimate, diagnostic and residual aggregate
must reproduce exactly under the **original strict contract**. That replay fails
exits1/1 on last-bit auxiliary reciprocal-condition estimates, preserved in
`failed-public-replay.tar.gz`. A fixed one-prefix diagnosis finds all science
fields and A/B/X entries exact; only eight rcond fields drift <=8.88e-16 absolute.

Extract `condition-support.tar.gz` into `source/` after original/support archives.
The [separate supplement](../../numerov-co-canary-condition-replay-supplement.json)
keeps scientific records exact and compares rcond observations separately at
32 binary64 epsilons relative. Run its successor module:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_co_canary_condition_replay --report ../collected/result/report.json --operands ../collected/result/operands --raw ../raw --output ../successor-verification.json
```

Original producer/strict verifier remain unchanged; `diagnosis.tar.gz` preserves
the bounded diagnosis and its raw comparison. Timing rows are not reexecuted.
The separate successor passes exits0/0: all303 saved decisions and all six fresh
scientific prefix records match exactly, maximum rcond relative drift4.55336e-16
within7.10543e-15. All19 final public payloads fresh-download and verify.

The [public replay contract](../../numerov-co-canary-public-replay-contract.json)
sets 90/180/300-s expected/upper/hard limits, one worker in reserved physical
cores0–3 and 4 GiB. Whole wheel and both nested initializer/binary hashes are
checked; no rebuild is part of the numerical replay. Early-prefix conditioning
and cost observations do not release a long-range interior comparison or alter
original interior false/exits1/1, historical132/594 coverage or model gates.
