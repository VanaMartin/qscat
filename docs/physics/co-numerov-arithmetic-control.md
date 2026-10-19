# Fixed-grid arithmetic-precision control

## Decision and unchanged physical problem

The [matrix attribution](../../projects/ukrmol_co/NUMEROV_PRECISION_ATTRIBUTION.md)
finds floating-point drift57–120 times larger than discrete truncation at the
failed compressed analytic200/400-bohr points. A channel basis change misses its
all-endpoint tenfold criterion. The next decision is whether more arithmetic
precision in **the same free constant-potential Numerov problem** removes that
drift, before changing spatial grids or returning to CO.

Keep the actual archived binary64 Q/threshold/energy inputs, sixth-order Taylor
initialisation/matching, identity-defect recurrence, zero angular momenta and
free tails. Internal coefficients use native Rydberg energies and bohr lengths.
The independent80/110-digit matrix discrete oracle is the finite-grid reference;
the continuous oracle remains a separate comparison. Retain every old failure.

## Precision-controlled implementation

Implement a small free-exterior diagnostic under`projects/ukrmol_co/`, with two
otherwise identical dtype executions: NumPy binary64 and NumPy longdouble.
SciPy LAPACK solves cannot retain longdouble precision, so use the same explicit
partial-pivot Gaussian solve in both executions. Require a longdouble significand
of at least64 bits, strictly wider than binary64's53; otherwise stop **unset**.
This is a precision comparison within the new typed implementation, not an
assertion that replacing SciPy alone repairs the old implementation.

For W=diag(threshold−E), B=I−h²W/12 is constant. Construct the Taylor **increment**
P(h)−I and derivative map D(h) directly, without subtracting identities. Initialise
`delta=(B change) [B(Q+change)]⁻¹` by a transposed solve. Each interval uses
`X=delta+diag(h² W/(1−h²W/12))`, then `delta=solve(I+X,X)`.
At the endpoint, match using
`A=B(P(−h)−I)+delta B`, `C=B D(−h)`, to exact free sine/cosine and decaying
columns. Column-equilibrate the final matching system; retain unsymmetrized K.
All arithmetic, trigonometric evaluation and matching stay in the chosen dtype.
The small solver is limited to at most three channels and is not a production
kernel or a general coupled-potential implementation.

## Finite prospective contract

Use six fixed identities from the previous compressed analytic study:

- Scalar h=.003125/400 and h=.0015625/800.
- Full rotated open/open/closed h=.003125/400 and h=.0015625/200/400/800.

Run both dtypes at all six: **twelve typed solves**. Reexecute the original SciPy
scalar finest800 and matrix finest400 once each, requiring exact stored K
identity: **fourteen new solves total**. No new spatial rung or endpoint.

Require80/110-digit discrete and continuous oracle drift<1e-20 before propagation.
Check the small solver with pivot-demanding systems and singular systems, and
short-grid scalar/open/closed recurrence against the independent discrete oracle.
For every typed row retain decimal K components and significand metadata without
rounding to binary64 before error assessment. Check finite values and physical
S unitarity/reciprocity at1e-6.

The extended-precision qualification gates are roundoff K<=1e-8 at every declared
row, and separate continuous K/S/phase<=1e-6 at every declared finest endpoint.
These are stricter than the original scalar roundoff1e-7 and preserve its1e-6
matrix/phase gates. Independently record the prospective hypothesis that extended
roundoff is at least100× smaller than the otherwise identical binary64 execution
at **all six** coordinates. A false hypothesis does not trigger another precision
or spacing rung. The coarser-row continuous errors are diagnostics.

Control failure/unset stops further submissions; accuracy gates are evaluated
after the fixed twelve typed rows. A qualification failure is preserved with its
original exit. No altered solver setting, tolerance, retry,
new CO invocation or full-grid successor follows automatically. The old SciPy
baseline and the typed binary64 implementation may differ in roundoff; both are
reported, never equated by assumption.

## Resource and evidence contract

One worker, reserved four physical cores0–3/4 GiB, single-threaded BLAS. The fixed
scope contains3957760 intervals. Explicit small solves add Python work, so use
a conservative **360-s forecast /600-s upper /1200-s cap** against the previous
61-s ten-solve measurement. Preserve0.2-s cgroup samples, original exits, dated
persistent `/home` evidence, pinned source/input bytes and before/after242 engine
artifact checks. Publish numerical outputs and independently verify saved
decisions; freeze any public numerical replay scope separately.
