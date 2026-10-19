# Interior stress selection: analytic prerequisite rejects

**Status: false / original exits1/1; the proposed CO successor is not submitted.**
The earlier [14-identity corner qualification](NUMEROV_CORNER_CONTROL.md) remains
valid in its tested regime. A stronger analytic prerequisite exposes limits
before the new compressed/equilibrium peaks are propagated.

## Selection and intended decision

The [read-only inventory](numerov-interior-selection-20261019.json) examines all
594 native32 matrices. In each sector it selects the maximum absolute interior
K entry (first-index tie break), plus fixed midpoint index48. The compressed
peak is index87/4.45 requested eV (K entry207.915, eigenvalue≈313); equilibrium
is index53/2.75 eV (entry26.710, eigenvalue≈−51.446); stretched index18/1 eV is
already qualified and is reused. Selection does not observe candidate errors.

The proposed new scope was four peaks plus six midpoint identities,60 solves,
4190160 intervals. A leading free-channel dispersion diagnostic is103.81 times
larger at the compressed peak than at the qualified stretched peak. This
motivated a prospectively fixed finer grid, validated analytically before any
CO use. See the [derivation, budgets and stopping rules](../../docs/physics/co-numerov-interior-control.md).

## Completed analytic result

Source/design/selection frozen and pushed at **`b7f2278`** before execution.
Twenty-four scalar and rotated open/open/closed controls cover the two new
peak energies, tangent+320/−55, three spacings, three endpoints and tail6/8.
Continuous and scalar discrete Numerov oracles use80 and110 digits.

**64/68 checks pass. Four reject:**

| Compressed-model analytic control | Measured K error | Budget |
|---|---:|---:|
| Scalar, fine800, versus discrete oracle | 1.20892527e-7 | 1e-7 |
| Rotated, fine200, versus continuous oracle | 1.12543592e-6 | 1e-6 |
| Rotated, fine400/order8, versus continuous oracle | 1.11113104e-6 | 1e-6 |
| Rotated, fine400/order6, versus continuous oracle | 1.11113104e-6 | 1e-6 |

Here fine h=.0015625 bohr. The order6/order8 repeated400 result is retained as
two frozen checks, not counted as two independent physical failures. The rotated
800-bohr continuous K error8.46436e-7 passes, as does scalar800 continuous K
error5.76310e-8; these do not replace the rejected comparisons. Rotated400 error
worsens4.92619e-7→1.11113e-6 when halving h. All S/phase/invariant gates and
oracle-precision checks pass. Equilibrium-model analytic gates all pass.

The scalar discrete oracle separates its rejected floating-point drift from
discretisation. A continuous error alone does not establish the coupled problem's
individual roundoff and truncation contributions. No monotone refinement cure
or absolute native-error bound is inferred.

## Preservation and resource outcome

Measured **72.952831 s /1.21588 min**, sampled CPU72.491243 s, recorded kernel
high-water **47,329,280 bytes /0.044079 GiB**,363 samples at0.2 s. Four reserved
physical cores0–3/4 GiB,120-s forecast/600-s cap. All46 archived sources and242
installed engine artifacts verify before/after. Remote evidence is under
`/home/kooza/ukrmol/co-sa-20261006/prepared/numerov-interior-analytic-20261019/`;
the10-file collection preserves original scientific exits1/1.

The [independent saved-decision verifier](numerov_interior_replay.py) reconstructs
every oracle comparison and all68 decisions using separate scattering/error
expressions; it also reconstructs the native selection and coverage from pinned
reference/history bytes. [Verification](verification-numerov-interior-analytic-20261019.json)
preserves **false**, not an upgraded result.

**The preflight stopping rule is enforced.** No new CO contract is frozen or
submitted, no CO/native/inner solve occurs, and no new unique energy is added.
Candidate qualification remains14; historical coverage remains132/594 with462
uncomputed. Four potentially new identities were selected but are **not computed**.
Further work needs a separately frozen fixed-grid precision/algebra or independent
propagation attribution on this analytic regime. An unsubmitted orchestration
draft is archived outside the repository; it is not a completed experiment.

## Full public numerical replay and final capture

The separately frozen [replay contract](numerov-interior-public-replay-contract.json)
reexecutes **all24 analytic rows**, not a subset, from fresh public downloads.
Every stored K component reproduces exactly; all68 check values and the **false**
aggregate remain identical. Verification exits0/0 preserve scientific exits1/1.
The producer and numerical replay share the same frozen propagator/model/oracle;
the separately expressed saved-decision verifier checks their interpretation.

Replay finishes **2026-10-19T09:53:03.602698UTC**, measured72.218529 s /
71.695078 sampled CPU-seconds /44,679,168-byte recorded kernel high-water,
359 samples at0.2 s, four physical cores0–3/4 GiB. Its8-file collection preserves
owner, complete reexecuted report, exact-component comparisons, exits and samples.
All **16 final public payloads** fresh-verify;62 focused tests and saved decisions
pass. See [public evidence and instructions](qualification-evidence/numerov-interior/README.md)
and [public verification](verification-numerov-interior-public-20261019.json).

The [09:54:14UTC remote capture](progress-numerov-interior-result-20261019.json)
reverifies18 collected hashes,46 frozen sources,242 installed artifacts and all12
native references plus history. No scientific container runs;61.49 GiB remains.
