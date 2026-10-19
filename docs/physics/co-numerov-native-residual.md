# Experimental native extended residual

## Measured target and API

The [mixed-solve profile](../../projects/ukrmol_co/NUMEROV_SOLVE_PROFILE.md)
attributes94.6% of extended recurring solve time to `B - A @ X`,27.5× combined
LU/initial/correction work. The existing NumPy expression in the validated
coupled implementation is the permanent differential oracle. This primitive
does not change the Numerov equation, spatial grids, units or backward-error
gate. Matrix entries keep the caller's atomic-unit convention.

The experimental Python interfaces are `numpy_residual(a, b, x)` and
`native_residual(a, b, x)` under`projects/ukrmol_co/`. For A(m,k), X(k,n), B(m,n),
return an unsymmetrized real(m,n) residual, with no input mutation. Native mode
requires identical native-endian NumPy longdouble arrays, finite values and
positive compatible dimensions. It accepts transposes, Fortran order and
negative-stride views by copying numeric values into contiguous buffers.

## Precision boundary and algorithm

Stable Rust does not expose the host's80-bit extended scalar. A small C shim
under`native/qscat-extended/` owns only longdouble arithmetic, called by a
dimension/length-checked Rust/PyO3 buffer interface. Both layers and the compiler
are explicit provenance; this is a Rust/C kernel, not Rust binary128 arithmetic.

Native qualification is restricted to little-endian x86-64, sizeof(longdouble)=16,
LDBL_MANT_DIG=64 and LDBL_MAX_EXP=16384, matching the measured NumPy format.
The wrapper checks the NumPy descriptor and native capability before use.
Other hosts retain the CPU-runnable NumPy reference. No binary64 coercion occurs.

The shim copies A into aligned storage and transposes X into aligned contiguous
dot-product rows. For each output, accumulate products in increasing k order,
with64-bit-significand longdouble multiply/add; subtract from B only afterwards.
Four output columns are accumulated together for instruction/data reuse, with
the same per-output k order. Compiler flags `-fno-fast-math -ffp-contract=off`
forbid reassociation/contraction. Checked dimensions bound allocation sizes;
copies avoid byte-buffer alignment assumptions. Output stores the ten significant
little-endian bytes and zeroes padding. Python returns a newly owned array.

This is an exact-arithmetic reimplementation, not a promise of bitwise identity
to every NumPy reduction. Differential qualification uses explicit backward-scaled
budgets and independent high-precision dot products. Caller-level residual gates
must still be evaluated independently when integration is eventually qualified.

## Primitive qualification before integration

Freeze the build, finite numerical scope and costs before remote execution:

- Six archived A/B/X captures: three extended native comparisons and three
  binary64 NumPy-reference controls. Reuse only their actual numeric values.
- Independent deterministic rectangular dense and cancellation cases, including
  a2^-60 signal lost by binary64, plus contiguous/Fortran/transposed/reversed views.
-100-digit mpmath dot references; preserve the true error separately from NumPy
  agreement. Maximum-entry error divided by `||A||inf||X||inf+||B||inf` must be
  <=32 extended eps; a zero denominator requires a zero residual.
- Reject incompatible dtype/shape/nonfinite inputs; preserve all input values.
- Measure end-to-end native/NumPy residual cost, including packing/transposition/
  allocation, on the three298-channel extended captures: one warmup and two fixed
  measured repetitions in reversed backend order. Prospective performance
  hypothesis: each capture improves at least1.5×. A false hypothesis is retained
  independently of numerical correctness, with no automatic tuning ladder.
- Release build, Rust unit checks, clippy/fmt, differential Python tests and two
  Criterion dense-size benchmarks. Save exact compiler/toolchain flags, lockfiles,
  wheel/binary hashes, original exits and measured resources.

The coupled propagator is not modified by this primitive contract. Integration
requires a separately frozen comparison including analytic/coupled stress,
independent residual checks, all prior numerical gates and measured interval cost.
No native CO or inner solve is authorized. Synthetic captured conditioning is
not arbitrary298-channel CO conditioning.
