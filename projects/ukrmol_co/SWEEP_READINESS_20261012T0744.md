# CO sweep readiness — 12 October, 07:44 UTC

**The separately declared CSF address audit is running on Sadaharu.** Owner
1187176 starts at07:39:28UTC; native B1 starts at **07:42:49UTC** after fresh
compilation, all **12 serial/MPI target-HF reference checks** and all three
debugger tests pass. At **07:44:50UTC**, all eight ranks have passed the
after-creation and before-build manager/parent/ID/orbital-count/determinant-bound
audits of **990,990 CSFs /22,475,718 determinants**, plus retained target/
orbital/options context checks. Each rank completes Class1 target index1 and
enters index2. No pointer-failure marker or scientific comparison is recorded.

These checks show initial address/index consistency in this instrumented run.
They do not establish general lifetime validity or the cause of the prior
rank0/Class1 SIGSEGV, nor its relationship to the release-engine rank4 crash
after Class4. Immediate Slater-rule parent-address guards remain active.

The [source investigation](NATIVE_POINTER_DIAGNOSTIC.md) and
[address-audit contract](cas12-pointer-diagnostic-address-successor-contract.json)
preserve the reasoning and three failed instrumentation attempts:

- Missing logging import: build/container2, controller1; no native invocation.
- Ambiguous patch anchor: build/container1, controller1; no native invocation.
- An audit predicate rejects three MPI reference cases despite equal addresses:
 9/12tests pass, CTest8, checks-container/controller1; no native invocation.

The [reduced predicate probe](pointer-association-probe-20261012.json) isolates
false `ASSOCIATED` for nonpolymorphic pointers to polymorphic array elements
at O0/O1/O2 in the pinned compiler. Scalar-address guards accept valid links
and reject redirected/null links at each optimization, with exits0 and no
warnings. Replacing that audit predicate restores the full engine gate;
this is an instrumentation correction, not a native engine fix.

All **55 source /192 input-evidence hashes**, **242 installed artifacts** and
the five saved native inputs reverify. Numerical-library hashes match the
prior debug build. The [immutable launch capture](progress-pointer-audit-20261012T0744.json)
pins31remote records plus its frozen observation, including per-rank marker
logs. Frozen source is `prepared/cas12-pointer-address-source-20261012/`,
pinned to `d88fbc9`; original algorithms, deck, libraries and tolerances persist.

Native allocation remains8physical cores4–11,80-GiB cgroup,8-GiB shared backing,
16-GiB `memp` and one library thread per rank. Native/container/owner limits
are24/25/30hours. The current **799.22% CPU /11.33GiB memory** sample is not
final cost. `STOP` in the source directory ends the active owner step.

The original observer starts before the owner's first execution file and
fails1 without captures. Its record is preserved alongside two earlier
observer startup failures. The [observation-only successor](cas12-pointer-observer-successor-contract.json)
starts only after checking the live owner file: **1191504** is active with
one immutable capture, hourly readings and a30-hour bound. It schedules no
calculation or repair.

**Next:** inspect any parent/manager/context/index/address failure and its
source stack, or diagnostic completion, before separately declaring an engine
correction or matched2048-root engine control. A changed engine needs baseline
reproduction before scientific4096-root comparison. CAS12 pair stability and
omitted-spectrum verdicts remain unset. Electronic/channel convergence, state/
pole identity, extraction, neutral correlation/domain and a frozen restartable
production pilot continue to gate the sweep and held-out potential fitting.
