# Bounded spectral-resolution control for the combined polar exterior

## Decision and prospective scope

The [original combined exterior study](../../projects/ukrmol_co/COMBINED_EXTERIOR_BENCHMARK.md)
rejects its frozen gate: degree128→256 changes by 7.90991e-7 against a
1e-7 budget. Exactly three configurations fail, all ε=30, L=80, θ=20°,
one per dipole including d=0. Their finest independent-method differences are
8.36e-14–4.09e-11; the full matrix's finest discrepancy is 2.654e-9.
This supports a resolution experiment rather than changing the tail physics.
It does not establish that the degree256→512 change will pass: differentiation
roundoff and conditioning may grow with degree.

The separately frozen [resolution contract](../../projects/ukrmol_co/combined-exterior-resolution-contract.json)
uses **degrees128/256/512**, with the entire original 432-case matrix repeated.
The [model, outgoing sheets, equations and methods](co-combined-exterior-benchmark.md)
and implementation are identical. Every energy, coupling, domain, contour and
ODE tolerance is retained; no failed configuration is removed. Degree128 remains
as the coarse rung to expose the known discrepancy. The four equal-threshold
analytic controls remain at the unchanged runner's fixed degree256; the zero-
dipole unequal-threshold analytic limits are evaluated throughout the new family.

All sixteen checks and numerical budgets are unchanged. In particular the
maximum degree256→512 entry change must be ≤1e-7 and improve by a factor of ten
over the maximum degree128→256 change, subject to the existing 1e-9 roundoff
floor. Finest independent-method agreement, analytic-limit agreement, domain/
contour spread, real-axis transpose symmetry and flux positivity remain required.
The two lower-half-plane energies remain on the declared outgoing sheets.

## Cost, preservation and stopping

The parent costs 11.178 seconds with 0.073 GiB sampled kernel high-water on
Sadaharu. Doubling degree multiplies dense cubic work by approximately eight
and dense storage by approximately four; the conservative estimate is **under
90 seconds**, with a **300-second timeout and 4-GiB limit**. Reserve physical
cores0–3 for this small job and use one BLAS thread. Freeze design, contract and
the unchanged runner's source hashes before execution; preserve persistent
source, dependency wheel, logs, exits, container metadata and sampled resources.

One fixed run is authorized. False or unset stops without retries or another
grid. Preserve the parent's false verdict and exits1/1 separately. Passing this
family qualifies only the defined toy outgoing DtN problem. The next decision
would be a prospective independent saved-native open-open K reconstruction,
including radius/step/tail and Gailitis normalization controls; it does not
authorize native pole interpretation or a production sweep.
