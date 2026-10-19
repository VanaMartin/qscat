# Recurring mixed-solve timing attribution

The [dense scaling benchmark](../../projects/ukrmol_co/NUMEROV_DENSE_SCALING.md)
passes166 checks and measures extended recurring rates.1353/.2714 s per interval
for free/varying W. It establishes recurring mixed solves as the priority, but
does not distinguish residual evaluation from binary64 LU. Resolve that decision
before implementing a compiled primitive.

## Fixed scope and instrumentation

Reuse the exact N=32 free/varying synthetic298-channel inputs, both arithmetic
formats, from the frozen scaling report. Run four untraced rows in free/varying,
binary64/longdouble order, then four traced rows in the reverse order: **eight
propagations /256 intervals**. All eight must reproduce the original repetition0
decimal/binary64 K and complete numerical diagnostics exactly. The original
coupled implementation and scaling input builder remain unchanged.

During each solve only, Python line tracing measures the time from an executable
line to the next line/return in `refined_solve`. Nested calls are included in
their parent's interval, not timed twice. Handler bookkeeping is excluded from
the line intervals and remains in the encompassing solve timer. Source lines are
recorded literally and their timing categories resolved against current source:

- `residual = b - a @ x`: residual evaluation, including allocation/subtraction;
- LU factorization and input cast; initial solve and casts; correction and casts;
- denominator/error norms; all remaining lines.

This distinguishes source statements, not a unique elementary operation or a
numerical-error cause. Save nonnegative line durations/counts, disjoint categories,
solve roles, actual input/output shapes/dtypes and encompassing wall times.
Require partitions to fit inside their encompassing timer with1e-6-s slack;
both tracer and original numerical bindings are restored on success/failure.

## Reusable operands and qualification

At recurrence interval16, save A,B and accepted X for each occurring recurring
role (one free defect, varying potential/defect), in both dtypes: **six captures**.
Copy after the solve timer; save after propagation. Captures are diagnostic
operands, not additional solves. NPZ retains actual host longdouble values;
record shape, dtype, measured significand bits and SHA256. Padding bytes are not
portable numeric identity; replays compare loaded values and dtypes, not NPZ byte
identity. The captures exercise the synthetic regime, not arbitrary CO matrices.

In addition to exact archived K/diagnostics, require correct frozen solve count,
roles and shapes, precision preservation, at most three corrections and the
unchanged32-eps backward-error gates. The archived qualification remains valid
only within its original scope. New instrumentation cannot authorize CO runs.

## Prospective timing hypotheses and cost

Report each model/backend and each recurring role separately. On each extended
model, test the distinct timing hypotheses:

1. Residual evaluation is at least80% of traced recurring solve wall time.
2. Residual evaluation exceeds factorization+initial solve+correction by at least10×.
3. Traced/untraced recurring solve wall-time ratio is within[.8,1.2].

False timing hypotheses are retained independently of numerical qualification;
no repeats or new counts follow. Summarize the source-statement optimization
target only if supported. Any implementation and accuracy/performance
qualification require a separate contract retaining the Python differential oracle.

One worker, reserved physical cores0–3,4 GiB, single-threaded BLAS,0.2-s samples.
Forecast55 s (two copies of the previous four N=32 rows plus bounded tracing/
capture overhead); upper90 s; cap180 s. Preserve source/report/242 installed
artifact hashes, original exits and persistent dated `/home` evidence.
