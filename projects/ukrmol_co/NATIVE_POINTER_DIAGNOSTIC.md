# CAS12 CSF pointer-identity diagnostic

The [first source-level diagnostic](progress-diagnostic-result-20261012T0701.json)
captures rank0 SIGSEGV in Class1 at `CSF_module.F90:506`, while reading
`this % csf % manager`. It occurs earlier than the release-engine post-Class4
rank4 failure. A common cause remains unproven.

## Source inspection and reduced probe

In the retained UKRmol-in 3.3.0 source, `CSFManager::create_csfs` assigns each
CSF its manager pointer, then assigns each determinant its parent CSF pointer.
Hamiltonian constructors retain polymorphic pointers to the CSFs, orbital
table, options and target CI. The main program's actual CSF/manager/target
objects lack explicit `TARGET` attributes, while the retaining dummy arguments
have them. This warrants investigation of lifetime, temporary arrays and
index provenance; it does not identify the observed crash's cause.

The [reduced probe](native_pointer_probe.f90), compiled in the pinned tools
image with the diagnostic flags, reports matching addresses for its retained
CSF array, manager, target array, parent backlink and manager backlink, both
with and without explicit `TARGET` on the actual objects. Its
[recorded result](pointer-interface-probe-20261012.json) preserves exit0,
compiler warnings and all output. Array `LOC` expressions warn about possible
temporaries. This small address probe does not reproduce the native object
graph, validate pointer lifetime generally or establish a `TARGET` correction.

The saved Class1 frame has `ido=62305` and `jdo=29198` against
`num_configurations=70746`; those printed loop indexes alone do not show an
out-of-range access. Other optimized/synthetic-pointer GDB values cannot
establish corruption. The native diagnostic must inspect the live backlinks
and determinant IDs before the failing dereference.

## Declared full-object audit

The [pointer-audit contract](cas12-pointer-diagnostic-contract.json) specifies:

- Audit every CSF's manager identity, parent backlinks, ID, orbital allocation/
  count and determinant bounds immediately after creation and before matrix
  construction. Check the retained CSF array and target/orbital/options context.
- Add flushed Class1 entry/exit markers for each target-symmetry index.
- Compare determinant parent addresses with the CSF passed to `slater_rules`
  immediately before each determinant dereference. On mismatch, print rank,
  CSF, determinant, determinant ID and actual/expected addresses, then
  `ERROR STOP 86`. Structural failures retain their source backtrace.
- Retain the earlier Class5/6/7 markers, GDB-per-rank capture, exact saved B1
  inputs, original libraries and `-O1 -g -fbacktrace -fcheck=bounds` flags.

`native_pointer_patch.py` requires four original source checksums and unique
anchors; `build_pointer_native.sh` selects this patch in the retained build
recipe. `native_pointer_owner.py` uses the frozen build/run helpers and requires
all twelve serial/MPI target-HF checks and three debugger smokes before one B1
invocation. The audit adds assertions and logging; it does not alter the
engine algorithm, model, target attributes or tolerances.

## Preserved instrumentation failures

1. Owner1179919 stops at **07:22:06 UTC** during compilation: the new
   Hamiltonian logging lacked a `const_gbl::stdout` import. Build/container
   exit2 and controller1 are preserved in the
   [failure capture](pointer-audit-build-failure-20261012.json).
   All50source/158input hashes reverify; no engine checks or B1 invocation occur.
2. The [import successor](cas12-pointer-diagnostic-build-successor-contract.json)
   stops at **07:27:50 UTC**, owner1182601. Its guarded patch refuses an ambiguous
   module-start anchor that also matches the end-module line. Build/container
   and controller exits1 are preserved in the
   [anchor-failure capture](pointer-audit-anchor-failure-20261012.json).
   All52source/167input hashes reverify; compilation, engine checks and B1 do
   not start.

The [separate anchor successor](cas12-pointer-diagnostic-anchor-successor-contract.json)
qualifies the module-start anchor with a preceding newline. All four source
pins, unique anchors, added Fortran line lengths and the single missing import
were checked against the retained source archive before launch. Each attempt
has a fresh source/build/output directory and keeps both earlier failures.

### Predicate isolated by the engine gate

Owner1182943 compiles successfully, then stops at **07:33:57 UTC** with9/12
reference cases passing. All three MPI target-HF cases fail the new parent
`ASSOCIATED(nonpolymorphic pointer, polymorphic array element)` check at
`after-create`, although scalar `LOC` prints identical actual/expected addresses.
The [preserved gate capture](pointer-audit-predicate-failure-20261012.json)
records CTest8, checks-container1 and controller1; no native B1 or debugger
smoke starts. All54source/174input hashes and242installed artifacts reverify.

The [independent reduced predicate probe](native_pointer_association_probe.f90)
reproduces false `ASSOCIATED` for all four valid polymorphic array elements
at **O0/O1/O2** in the pinned compiler. Its nonpolymorphic counterpart returns
true, while scalar-address checks accept the valid links and reject deliberately
redirected and null parents. All three runs exit0 without compiler warnings;
the [raw probe record](pointer-association-probe-20261012.json) preserves commands
and output. This diagnoses an audit-predicate flaw, not the native SIGSEGV.

The [separate address successor](cas12-pointer-diagnostic-address-successor-contract.json)
replaces only that parent predicate with scalar `LOC(parent) /= LOC(expected)`,
matching the immediate Slater-rule guard already declared. It retains the full
engine gate, exact physics inputs and earlier failed attempts.

## Allocation and interpretation

Native B1 retains eight MPI ranks on physical cores4–11, one library thread
per rank,80-GiB cgroup,8-GiB shared backing and16-GiB native `memp`. Native,
container and owner limits are24/25/30hours. The hourly observer is bounded
to30hours. `STOP` ends the active owner step; failures schedule no retry or
follow-on calculation.

Review the audit stage, CSF/determinant/index/address tuple and source stack
before declaring any engine correction or matched2048-root engine control.
Successful diagnostic construction would still require baseline reproduction
before resuming the scientific4096-root comparator. Pair stability, omitted
spectrum, model, continuity, extraction and neutral-reference qualification
remain distinct gates on the production sweep and potential fit.
