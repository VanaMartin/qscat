# A saved-input compressed A2 outgoing-pole contrast

## Decision served, not a production-model selection

The [finite compressed open-Pi pair](../../projects/ukrmol_co/PAIRED_POLE_CHARACTER.md)
supplies one CAS11/nine-A2 pole. It does not establish electronic-model accuracy.
The smallest matched **development-geometry** contrast already present is
CAS11/five versus nine retained A2 states per spin at R=1.9bohr. Both2048-root
inner calculations are complete. Change the observable from historical phase-fit
parameters to the **same coupled outgoing pole definition** on both models.

The production-relevant decision is which uncertainty to address before spending
more on electronic models: do retained A2 channels materially affect the pole
even at the compressed anchor where the historical phase gate passed? A resolved
rejection prioritizes channel-model evidence before a new basis/active-space
allocation; a pass returns priority to a separately justified basis/active-space
contrast. Neither outcome selects a production model or proves nine-state
convergence. A pass must **not** downgrade to five states across the domain:
the stretched phase rejection remains. This cheap contrast is triage for the
next expensive model allocation, not a universal sweep prerequisite.

Other choices are deferred:

- The completed same-engine CAS11→CAS12 scattering contrast is at equilibrium.
  Its old phase/fit results remain known, but **the equilibrium pole stays
  excluded** from model selection and scalar fitting.
- The inspected development-geometry CAS12/basis records document target imports,
  not a compatible completed N+1 outgoing-pole pair. A target import does not
  supply scattering eigenpairs. No global filesystem-absence claim is made.
- New50 is parked. Its lowest-five contrast mixes orbital/channel prescriptions;
  it is not an isolated averaging change. No new QC, properties or inner solve
  is justified merely to avoid reusing the completed A2 control.

## What the inventory actually establishes

[Remote inspection](../../projects/ukrmol_co/progress-model-contrast-inventory-20261025.json)
and [fresh public-byte/metadata verification](../../projects/ukrmol_co/verification-model-contrast-inventory-public-20261025.json)
authenticate the scoped native records. Five/nine use40/48 selected targets,
250/298 boundary channels and2048 exported poles each. Both have the same six
ground-target open labels, radius18bohr, two CF orders and **zero Buttle flag**.
There are244/292 closed channels; these are not scalar p-wave models.
Use each model's **own** exported eigenvalues and boundary amplitudes: the
contracted dimensions also differ (27546/27818). Slicing the nine-model channel
matrix is not the completed five-model calculation and cannot replace it.

The40-row target table is an exact energy/symmetry/spin subset of the48-row
table in both sectors. Recorded native input hashes for integrals, target data,
properties and sector CONGEN match; outer decks are byte-identical. This is
recorded shared-input provenance, not a fresh multi-GB integral rehash. Historical
qualification already audited those original shared inputs. Binary R is zero:
the physical1.9bohr anchor comes from the original hashed run/deck contract,
not an invented geometry decoded from that field.

### Exact missing deliverable: an explicitly mapped250-channel model view

The strict saved reader still rejects five-A2 channel labels. Its original six
unset cases and aggregate failure stay unchanged. Channels retain original
property IDs, while the selected target table is densely numbered. Compressed
selection in written row order is `[1,...,38,40,42]`;12 channels per sector need
remapping, six have IDs above40. An in-range-but-wrong dense interpretation can
change a threshold by0.015165426Ryd; clipping IDs is not a repair.

For sorted selected property IDs `s`, explicitly define `dense(s[j])=j+1` and
check, for **every** channel,

    stored_threshold_Ryd = 2 * (E_selected[dense(property_id)] - E_ground).

Inspection reproduces the target-table subset; the fresh public reading
authenticates those bytes and independently checks every mapped threshold.
That metadata result is **not a qualified new
decoder, real-S reconstruction or pole calculation**. The existing
`native_full_reconstruction.load_model` deliberately supports48/298 only.
The concrete next implementation is one private compressed-five mapped view,
retaining original labels as provenance and exposing dense lookup separately,
plus a narrow case/loader adapter to the existing pole runner. Do not silently
relax the strict reader or overwrite native files, and do not build a new solver.

### Why the raw IDs do not require an engine patch for this path

The authenticated source trace is specific to the saved fixed-nuclei K path:

- `Postprocessing_module.F90::get_channel_info`,748–778, keeps original property
  ID/energy in ICHL/ECHL; `::write_channel_info`,1070–1083, separately densifies
  selected target rows. `outerio.f::READTC`,1295–1305, reads both fields as stored.
- `rsolve_compak.f`,160–166, defaults `NGEOM=1`, `NPOLE=0`, `CALCAK=0`;
  the saved decks do not override them and explicitly have `CALCDIP(1)=0`.
- The ICHL-indexed target-energy dereference in `RESIDR_compak`,1616–1632,
  is guarded by `abs(IBUTTL)==1`. All inspected boundary headers haveIBUTTL0,
  unchanged by `newbut=0`. The pole sum instead uses total energy and W/D.
- The adiabatic-nuclei call is guarded by `NPOLE>0` at654–677, and wavefunction/
  dipole branches by their disabled flags at503–543 and845–873.
- Active asymptotic matching receives **stored EVCHL thresholds** at685–689,
  followed by propagation and open-open K output at719–745; it does not rebuild
  those thresholds from an invalid dense target pointer.

Source hashes/excerpts are in the inventory and published raw/writer archives.
This traces the relevant pointer-sensitive branches; it is not an independent
compiled-engine test or qualification of other consumers, Buttle corrections,
photoionization or vibrational paths. The future mapped view must still reproduce
two development full-S matrices before any complex evaluation. Failure stops;
no native repair, rebuild or fresh scattering solve follows automatically.

## Explicit implementation API and gates

`read_interface(..., selected_property_ids=...)` sorts a checked integer bijection
and adds separate `channel_target_rows`/`selected_property_ids` metadata. Raw
channel labels stay intact. With no keyword it retains strict dense-row behavior
and the same historical failure. Mapped thresholds must agree within1e-10Ryd;
selected written row IDs must be consecutive. No bytes are rewritten.

`compressed_five_view.load_compressed_five(source)` accepts only compressed/five
B1/B2, checks original pilot versus actual OUTER_INTERFACE selection and the
48-state target subset (1e-10Ha), then shares native assembly with the original
48/298 loader. Its W/D/CF come from the **five** files. The old loader remains
48/298-only by default; the shared assembly permits just these two declared layouts.

`compressed_a2_pole_contrast` supplies the narrow case dispatcher, proposal-to-
pole-contract translator and `pole_insensitivity` band/identity decision. There
is no automatic-execution CLI. The existing pole runner requires a checksum-
pinned real-guard report covering exactly indices60/80, recomputes its measured
error gates and only then permits a five-model complex use. A canary pass for
the nine model cannot substitute for these two real checks. All ODE, bordered,
contour and derivative numerical routines are unchanged.

`compressed_five_view_gate` is a separate **saved-only** CLI. After source/tests/
contract freeze it may authenticate the original inputs, preserve the default
strict failure, check own W/D bytes and nine-model identity-view equality. It
does not construct S, propagate, search or release a pole job. Source/tests and
an exact owner for a numerical contrast still need a separate execution freeze.

## Prospective staged test — not released for execution

The [proposal data](../../projects/ukrmol_co/compressed-a2-pole-contrast-proposal-20261025.json)
are fixed before any five-model pole evaluation. Implementation, controls, exact
source/inputs and owner must separately freeze/push before numerical application.

1. **Saved-only mapped-view gate:** verify selection IDs against both original
   pilot and actual IDTARG deck; require a unique bijection, target metadata/
   thresholds and exact250/40/2048 layout. Preserve all raw labels/bytes and the
   old strict-reader failure. Test unsorted selection, ID42→dense40, wrong
   in-range lookup, clipping, duplicate/missing IDs, threshold corruption and
   identity mapping on the48-target reference. No energy-nearest guess.
2. **Five/B1 first:** at development indices60 and80, reconstruct full S at the
   final outgoing setting. Retain the existing1e-6 entry/invariant and phase-sum
   budgets,1e-7 flux budget. No holdout index is used. Only after both pass,
   search one local box around the saved nine/B1 pole: +/-0.01Ha in each
   coordinate, original seed and +2e-5Ha. This allocation is not a physical
   admissible interval. Keep all250 channels and the uncorrected2048-root sum.
3. **Same original witnesses:** settings80/30/order8/rtol1e-10,
   120/30/order8/rtol1e-11 and120/40/order10/rtol1e-12; root caps16/12/12;
   fixed radius.001Ha /32 contour points, nested16; eight probes /seed314159;
   derivative steps1e-5 and5e-6Ha. Keep every original boundary-null, winding,
   interior-resolvent and visible rank-one open-residue gate. Use the existing
   phase-free Pi/support and four-witness stability policy unchanged.
4. **Identity caution:** label-match the six-channel intensity profiles with the
   nine-sector reference. Cross-model TV above.25 is **inconclusive identity**,
   not a sensitivity rejection; below it is only a necessary local identity
   screen, not continuous tracking. This conservative25% redistribution cutoff
   is a declared assignment policy, not a symmetry theorem. Include a synthetic
   wrong-Pi-pole/different-l-profile control. Report residue/norm changes, but do
   not require cross-model residues to match: the scalar target omits them.
   Ambiguous identity, box exit or multiple directions stops; no homotopy ladder.
5. **Only if B1 is a valid insensitivity pass:** evaluate five/B2 using the same
   real checks and the prior paired-sector allocation: seed the five/B1 pole,
   +/-0.002Ha box, root caps12/8/8, unchanged witnesses. Require five-model B1/B2
   position/full width agreement within0.1meV and the original paired character/
   canonical-profile gates. If B1 already gives a resolved sensitivity rejection,
   stop without spending on B2; report that sector rejection and leave the new
   five-model paired assignment unset.

No nine-model propagation is repeated; reuse its authenticated pole/residues.
No legacy fitted width or new rational continuation is substituted for a pole.
Use `E_total=E_ground+z_Ha`, exterior incident energy `2*z` Rydberg and canonical
`qscat.units` for eV. Keep stored CF scaling and all closed channels; no extra
half/radius/conjugation in the native boundary sum or double polar correction.

## Decision bands and uncertainty boundary

For each qualified sector compare `Re z` and full width `Gamma=-2 Im z` with
the **nine-model reference**. Position tolerance is.05eV, width tolerance
`max(.001eV, .05*Gamma_nine)` (about.10423eV here). Do not compare percentages
of half widths or select background/window subsets.

Retain.0001eV same-root position/full-width numerical budgets per model.
Use operational boundary bands.0002eV for the position contrast and.000205eV
for width (two model allowances plus5% of the reference-width allowance):
pass only below `tolerance-band`, reject only above `tolerance+band`, otherwise
inconclusive. These guard bands are a prospective decision policy, **not rigorous
error enclosures**, and do not cover electronic or omitted-spectrum uncertainty.
Test below/inside/above bands and a full-width/half-width negative control.

If both sectors qualify and pass, report only finite compressed A2 pole-target
insensitivity at the tested settings. Rejecting a valid B1 contrast is sufficient
to reject that insensitivity claim. Any failed prerequisite or unresolved pole
identity is inconclusive, not a physics rejection. Identical accepted seeds
are not independent iterates, exact degeneracy or zero uncertainty. Constructions
share D/W and the outgoing ODE; no independent exterior/absolute CO accuracy.

## Costs, effort cap and stopping

The original five/nine two-sector calculations cost1.651263/1.827292h and
6.011147/7.245762GiB. **Do not repeat them:** the inner eigenpairs are present.
The metadata inspection takes.0272s on the remote host; SSH/download/hash work
is separate. Two fresh public payloads independently authenticate26 scoped
input pins and the channel/target metadata, without S or pole construction.
One preliminary metadata-only path typo exits1 and is recorded, not relabelled
as a scientific failure or hidden by the successful corrected inspection.

Proposed mapped-view checks: four physical cores0–3, one thread,4GiB,60s hard cap.
Future propagation: at most **79 maps B1 +67 B2 =146 total**, eight physical
cores0–7 /one BLAS thread per worker /up to8 contour workers /8GiB. Prior298-channel
47/42-map applications cost181/173s, not a measured250-channel forecast.
Conditional target600s /upper1800s total, hard1800s per sector /3600s combined;
90s per map. No paid provisioning, image/dependency rebuild or background owner.

Implementation/diagnosis: **one active cycle, at most two hours** for the explicit
mapped view, narrow adapter and controls. If they cannot freeze within that bound,
stop at this exact prerequisite rather than extend an interface project. Any
input, normalization, cap, identity or witness failure preserves partial artifacts
and original exits; no retry, expanded region, successor, precision/channel/root
ladder or fresh native solve. Review immediately after the first gate/result.

This design cycle launches **zero** native calls, propagations or pole evaluations.
All114 route holdout S matrices and the equilibrium pole-target exclusion remain
protected. New50 stays parked. Model accuracy, full channel/basis/active-space
convergence, geometry continuity, neutral correlation, factory nodes/asymptote
and physical full-S embedding remain unresolved. No curve, fit or sweep release.
