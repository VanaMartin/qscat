# Preserve the zero-atom stop; correct only the source-equality control

The first5Z launch stops12:01:12UTC before any integral/SCF/CC call. Seven runtime
controls pass; the calculation-equality AST hash fails. Local Python3.12 and the
actual image's Python3.11 serialize different FunctionDef metadata. The actual
QZ and5Z `calculate` function source bytes are identical. This is a software
control-portability error,not an atomic convergence or physics rejection.

Preserve the original14-file capsule/archive,all original exits(control1,
command/container/owner/service/launcher2),and absence of atomic outputs. Do not
rerun or requalify that source/contract. Its archive SHA256 is
`69a94cd4509aa530dac03a89e50399cdc9497ba0e2dc6bd3533f2bc391763085`.

The separate corrected release changes **only the test/contract dispatch** to
require an exact `inspect.getsource(calculate)` byte hash. The expected hash is
derived from the executed QZ Git source,and current bytes compare equal:
`a07c02721c493ff6f289c6972fb411a0694f4850ed148ec317d7e9e77a146cd0`.
Source-byte equality is stronger than AST equality and independent of AST
metadata serialization. No numerical,spin,state,occupation,angular,stability,
memory,timeout or finite-basis diagnostic gate changes. No quantum job repeats.

Freeze/push new test/contract/source manifest with a distinct capsule/container.
Run all eight no-quantum controls **inside the actual pinned image** before
the separate science allocation. A failed control stops again with its original
output. If they pass,execute the same four-case5Z scope under unchanged bounds,
then independently reconstruct the actual result. This narrow implementation
repair serves the already user-approved endpoint control,not a new physical
branch,larger basis or fallback. No further repair ladder is released.
