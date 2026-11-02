# Stretched finite DZ/TZ local pole candidates are not physically insensitive

Both [own-DZ](../stretched-dz-pole/README.md) and [own-TZ](../stretched-tz-pole/README.md)
finite outgoing poles independently qualify. The authorized saved-only review
`fbf10be` reconstructs both complete public bundles, authenticates native open
labels and compares their saved high-rule residues. Four freshly downloaded
review payloads reproduce exactly. No new ODE/native/electronic calculation,
geometry or114-holdout analysis. Both original scientific successes and every
earlier failure retain their original exits; review0 is not a production pass.

|R=2.5bohr, own CAS12/64/400/B1/2048|DZ|TZ|TZ−DZ|
|---|---:|---:|---:|
|Position above each own neutral /eV|0.914117093|1.017487033|+0.103369940|
|Full width /eV|0.267092530|0.311264209|+0.044171679|

**Candidate position shift103.369940meV exceeds50meV. Full-width shift
44.171679meV (16.537969% of DZ) exceeds5%/1meV-floor=13.354626meV.** Using TZ
as reference also rejects:14.191056%>5%, width budget15.563210meV. Physical
comparison uses these budgets, never the same-model0.1meV numerical tolerance.
Even subtracting two declared numerical budgets leaves both shifts above their
physical budgets; this is a robustness diagnostic, not a rigorous error enclosure.

## What supports a local counterpart, and what remains unknown

Both models share the same40-component orbital-average objective/6/3/3/0 active
symmetry and eight roots per spin/sector, but use separately optimized DZ/TZ
frames. Equal counts do not prove target-state or wavefunction correspondence.
The six own-native open labels match exactly and represent the ground-target
channel with |m|=1/3 partial waves. Each own seed converges to a locally isolated
rank-one open-Pi visible pole, with unit winding and stable derivative/contour
residues. Each0.001Ha circle counts a local pole only; neither run counts the
whole union of search boxes or excludes other candidates globally.

The phase-free open-residue profiles are compatible with the same local
ground-channel Pi-dominated feature: first two components DZ0.781974/0.217521,
TZ0.755496/0.243882; total-variation change0.026477666. Normalized residue-magnitude
distance0.044440327, norm ratioTZ/DZ1.163826. These diagnostics are invariant to
channel signs and are **reported without an invented cross-model identity
cutoff**. Existing same-model witness-stability thresholds are not applied
between bases. Residue weights are not probabilities; matrix magnitude agreement
is not a coherent wavefunction overlap, global assignment or proof of identity.

**Inference:** necessary saved observable evidence shows no local-character
contradiction, and these compatible local candidates differ beyond both physical
budgets. This rules out using this comparison to claim basis-plus-frame
insensitivity. A strict same-electronic-state basis-convergence verdict remains
unset without cross-basis physical wavefunction/subspace correspondence. No
absolute CO error, larger-basis convergence, global geometry identity, Pi pair,
width curve, potential fit or production model is established.

## Campaign disposition: stop this insensitivity route; redirect model choice

Do not fit a production potential to DZ or promote TZ merely because it is
larger. Retain both poles as finite exploratory evidence and nine as the working
reference, with new50 parked,2048 baseline/equilibrium exclusion protected.
The missing production prerequisite is now an **explicit electronic-model and
remaining uncertainty decision**, not another outgoing-pole implementation or B2
rescue. Channels/omitted spectrum, cross-basis identity, geometry continuity,
neutral correlation/wider-domain reference and full-scattering embedding remain
open. A common-ground energy origin does not remove these uncertainties.

Next bounded deliverable is a saved-evidence electronic-route choice: at most
30minutes design, zero QC/CI/native/maps/holdouts, name the model uncertainty and
one finite next allocation only if it could change production selection. No
allocation, channel/root/basis ladder, geometry/fit/sweep is authorized here.
Any new numerical branch requires a separate explicit decision and freeze.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-dz-pole
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-tz-pole
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-pole-counterpart
uv run --no-sync python -m projects.ukrmol_co.stretched_pole_counterpart_review projects/ukrmol_co/qualification-evidence/stretched-dz-pole projects/ukrmol_co/qualification-evidence/stretched-tz-pole
```
