# Saved-matrix diagnosis of the stretched real-axis stop

## Scope and question

This is the bounded source-and-saved-evidence follow-up to the
[original guard stop](../../projects/ukrmol_co/qualification-evidence/stretched-real-guard/README.md).
It performs no exterior ODE integration, native call, root search, target solve
or protected holdout analysis. The original six scientific/wrapper exits remain2;
the guard is not repaired or retrospectively qualified. Both own native exports
retain their separate scope. No physical basis comparison has yet been made.

The question is whether the failure is consistent with a badly conditioned final
linear solve, a units/metadata error, or loss of the physical amplitude identity
in a particular open–closed sector. The
[diagnostic reader](../../projects/ukrmol_co/diagnose_stretched_real_guard.py) first
re-authenticates the original33-file evidence, native K/surface sum and failed
gates through the independent reader, then measures only the saved Q,Y,A blocks
and asymptotic rates from the authenticated native thresholds.

## Observed localization

The open Jost amplitude matrix A has6 rows and400 boundary columns. All momenta
here use the native Rydberg radial equation; lengths are bohr. Q uses the own
Hartree poles at total energy E0+E_Ry/2, as already independently verified.

| Saved development index |60|80|
|---|---:|---:|
|Native radial energy /Ry|0.22785|0.30135|
|Condition number of A|1.07219|1.05238|
|Condition number of I−QY|9.88730|11.42059|
|Maximum amplitude magnitude|0.686316|0.737188|
|Open–open flux defect|6.79725e-10|6.29209e-8|
|Open–closed flux defect|2.36840e-10|**6.96992e-6**|
|Closed–closed flux defect|6.38625e-13|2.06123e-8|
|Wronskian residual|6.43222e-10|**1.49092e-5**|

The failed80 flux maximum is at zero-based indices(0,65). The Wronskian residual
maximum is at(65,0) at both energies. Column65 is native boundary channel66,
label[66,12,0,0]: target12 with radial angular labels l=0,m=0. Its target label
[12,1,1,0] and threshold **0.529019366Ry** are retained verbatim; no new electronic
state-identity or exact-Lambda claim is inferred from those labels.

Thus this is not evidence of a catastrophically conditioned final six-open
amplitude or I−QY matching solve. Nor is it large cancellation of huge saved
amplitudes: the largest saved A entry is below1. The source's pre-integration
amplitudes can have very different scaling, so these observations do **not**
bound conditioning along the entire contour or prove a unique cause.

The open-only part and final S are comparatively accurate while the open–closed
physical flux identity fails. That explains why S agreement, unitarity and
reciprocity can pass without satisfying all physical-amplitude prerequisites.
Do not drop the closed-channel columns, whiten A, enforce flux by projection,
or reinterpret the smaller open-only defect as the original full gate.

## A source-supported finite-tail contamination mechanism

The executed `outgoing_map.outgoing_map` initializes the finite inverse-power
tail F(z) and its logarithmic derivative Y=F′F⁻¹. Its open inverse-Jost rows carry
weights proportional to exp(−ikz), and are propagated by A′=−AY. This is distinct
from propagating only the outgoing subspace Y. The original numerical code is
unchanged, SHA256 `8860804fbc55ecee08d017937b5e5d116e7b52003d76cd42252ab13dc2923140`.

For an asymptotically diagonal open channel with k=sqrt(E) and a coupled closed
channel with κ=sqrt(T−E), consider a small unwanted **closed-dual component** in
an open inverse row at the finite endpoint. Along z=R+L exp(iθ), its magnitude
relative to the desired open dual, after backward transport, scales as

    exp[−L(κ cosθ − k sinθ)].

This follows directly from the diagonal dual equations A′o=−ik Ao and A′c=κ Ac:
their backward magnitudes shrink as exp(−kL sinθ) and exp(−κL cosθ), respectively.
It is an asymptotic sensitivity factor, **not a rigorous error estimate for the
coupled finite CO exterior**. Finite-tail coefficient errors, coupled evolution,
ODE discretization and roundoff are not separated by this observation alone.

At the originalθ=40degrees andL=120−18=102bohr:

|Quantity|60|80|
|---|---:|---:|
|k /bohr⁻¹|0.47733636|0.54895355|
|κ /bohr⁻¹|0.54878900|0.47714711|
|κ cosθ−k sinθ /bohr⁻¹|0.11357087|**0.01265535**|
|Critical angle atan(κ/k) /degrees|48.98327|**40.99697**|
|Relative damping factor|9.31180e-6|**0.275038**|

At80 the original contour lies only about1degree below this component's
asymptotic crossover, leaving very little attenuation of a finite-tail closed-
dual error. The relative damping factor grows by29536.5 between60 and80; the
observed Wronskian residual grows by23178.9. This is consistent in scale and
location with finite-tail closed-dual contamination. It is **supporting evidence,
not a fitted error bound, unique causal proof or a repaired guard**. A good
outgoing-subspace Y and a well-conditioned final S solve do not by themselves
qualify the full amplitude-row normalization.

## Decision and bounded next evidence

The saved-only diagnosis narrows the blocker to open–closed amplitude transport
and its finite-tail sensitivity. Retain the numerical stop; do not change the
CO angle, endpoint, tolerance, RHS budget, seed or channel model under the stopped
contract. No new CO map, complex pole, B2, fit or sweep is released here.

The next smallest method-level control is a separately frozen, exactly soluble
two-channel **dual-transport toy**, without CO input access: inject a prescribed
finite endpoint error into an inverse-Jost row and compare diagonal propagation
with the analytic exponential above on angles below/at/above the crossover.
Use that to distinguish the predicted transport sensitivity from an algebraic
sign/normalization error before proposing any CO numerical change. It does not
test the actual CO tail truncation or establish that another angle repairs it.
Design/test/source freeze and a separate execution decision must precede ODE
tests; no automatic numerical successor follows from this saved-only report.
