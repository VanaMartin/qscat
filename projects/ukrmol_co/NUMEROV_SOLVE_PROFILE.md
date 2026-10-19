# Recurring mixed-solve profile and reusable operands

The [design](../../docs/physics/co-numerov-solve-profile.md) and
[finite contract](numerov-solve-profile-contract.json) were frozen/pushed at
**`bbbc3d8c64fde2322910bfecc8d2a6eeb3bee86e`** before execution. The control
passes **76/76 checks / original exits0/0**. Eight reexecutions of the previous
N=32 synthetic298-channel rows add256 radial intervals and no new CO/native/
inner solve or unique energy. Coupled numerical modules remain byte-identical
to`4a9baa1`; the scaling input builder remains identical to`1d8cfac`.

## Exact numerics and measured statement cost

Four untraced rows run free/varying and binary64/longdouble, followed by four
traced rows in reversed order. **Every decimal/binary64 K component and complete
numerical diagnostic reproduces the archived scaling repetition0 exactly**.
All shapes, arithmetic formats, solve roles, bounded correction counts and
unchanged32-eps backward-error gates pass. Measured significands are53/64 bits.

The tracer measures parent-frame executable-line intervals in the unchanged
`refined_solve`. Nested call time is included once in the parent statement.
Handler accounting is excluded from line intervals and remains in the enclosing
solve timer. Recurring extended solve attribution:

| Model | Recurring solve / s | Residual evaluation / s | Residual fraction | Residual / (factor+initial+correction) | Traced/untraced ratio |
|---|---:|---:|---:|---:|---:|
| Free | 4.172260 | 3.946039 | **94.578%** | **27.5415×** | .999729 |
| Varying | 8.334918 | 7.889105 | **94.651%** | **27.5711×** | 1.000074 |

All six timing hypotheses pass independently of numerical qualification:
each residual fraction>=80%, each residual/LU ratio>=10 and both trace ratios
within[.8,1.2]. Varying potential and defect roles separately take4.162012/
4.172906 s; residual statements take3.942917/3.946188 s. Binary64 residual
fractions are about22%; its initial solves dominate instead.

The measured target is **`residual = b - a @ x`**, including matrix product,
subtraction and allocation. This is source-statement timing, not isolation of one
elementary instruction or a cause of CO/native numerical error. The next
optimization should target this extended residual path, retaining the validated
Python implementation as its differential oracle.

## Six captured solve operands

At recurrence interval16 the producer retains actual A,B and accepted X for
free defect and varying potential/defect, in both arithmetic formats. These
six captures have298×298 arrays, measured shape/dtype/significand metadata and
SHA256 pins. Copies occur after solve timing and files are written after each
propagation. Captures add no solve.

The public results archive contains`result/operands-<model>-<backend>-<role>.npz`.
Load with`numpy.load` on a compatible host. Host longdouble padding is not a
numeric identity; compare loaded values and dtype, not NPZ byte identity across
machines. These synthetic matrices exercise a small active block within a dense
propagation and do not qualify arbitrary CO conditioning or cancellation.

The independent [saved verifier](numerov_solve_profile_replay.py), frozen at
`f396ed3`, recomputes all76 decisions,392 solve-call shape/dtype/role audits,
source-statement timing partitions and all timing summaries/hypotheses. It
checks the six capture hashes/metadata. Actual wider-precision operand values
and accepted backward errors are checked by the separately frozen full
[public replay](numerov-solve-profile-public-replay-contract.json).

## Resources and next implementation boundary

Producer finishes **2026-10-19T14:07:13.855601UTC**: **49.634564 s /49.061244
sampled CPU s /155,426,816-byte kernel peak**,246 samples at0.2 s, one worker,
reserved physical cores0–3/4 GiB, single-threaded BLAS. Below55-s forecast/
90-s upper/180-s cap. All71 frozen sources,242 installed artifacts and pinned
dense-scaling report verify before/after. Original evidence remains at
`prepared/numerov-solve-profile-20261019/` under the campaign root on Sadaharu.

Next design and benchmark a precision-preserving compiled residual primitive on
the captured operands plus independently specified dense cancellation/stride
controls. A binary64 BLAS product alone cannot substitute for a64-bit-significand
residual. Resolve the native arithmetic/NumPy format boundary explicitly, preserve
the Python fallback and backward gates, and measure a speedup before integrating
and requalifying the coupled propagator. The six captures avoid another
propagation just to obtain representative operands.

Original interior **false/exits1/1**,14 candidate-qualified CO identities,
historical132/594 coverage and462 uncomputed remain. Full native reconstruction,
complex poles and production remain unqualified.

[Public evidence](qualification-evidence/numerov-solve-profile/README.md) ·
[Saved verification](verification-numerov-solve-profile-saved-20261019.json) ·
[Public verification](verification-numerov-solve-profile-public-20261019.json).

## Complete public numerical and operand replay

The separate verifier/replay scope was frozen/pushed at **`f396ed3`**. Complete
reexecution from fresh public bytes finishes **2026-10-19T14:10:51.257197UTC /
exits0/0**: all eight decimal/binary64 K matrices, full numerical diagnostics
and76 decisions reproduce **exactly**. All six captured A/B/X triples reproduce
loaded numeric values and dtypes exactly and pass independently recomputed32-eps
backward errors. Largest extended capture backward error is4.21951e-20.

Fresh timing observations again pass all six hypotheses: extended residual
fractions94.588%/94.664%, residual/LU ratios27.56/27.62. These timing identities
are observations, not numerical replay gates. Replay costs **49.832828 s /
49.263742 sampled CPU s /155,291,648-byte kernel peak**,247 samples at0.2 s,
reserved cores0–3/4 GiB, below55-s forecast/90-s upper/180-s cap. All71 sources,
one separate verifier,242 installed artifacts and prior report verify before/after.
The14-file replay collection (including six NPZ captures) remains under
`prepared/numerov-solve-profile-public-replay-20261019/`.

All **14 final public payloads** fresh-verify. Local public-source tests pass
**89 with two wider-format skips**, exercised by remote producer/replay. The
[14:11:53UTC capture](progress-numerov-solve-profile-result-20261019.json)
reverifies29 collected hashes,71 sources,242 installed artifacts and prior report.
No scientific containers remain;60.67 GiB persistent disk is free. Public:
`https://data.qscat.org/co-numerov-solve-profile-20261019/`.
