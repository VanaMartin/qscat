# Retained A2 channel audit — 11 October 2026

The [public quarter replay](qualification-evidence/quarter-state-tracking/README.md)
independently reconstructs the matrices used here. The
[root-resolved audit](a2-channel-coverage-audit-20261011.json) checks the two
remaining retained routes in both physical and AO-following Lowdin metrics:

- Singlet A2 root4: **4→4→5→6→6**.
- Triplet A2 root5: **5→6→6→5→5**.

## Coverage result

| Physical link / incoming root | Projection squared into roots1–5 | Into roots1–6 | Into roots1–12 |
|---|---:|---:|---:|
| Singlet midpoint→upper quarter / root5 | 0.000432 | 0.550711 | 0.872660 |
| Singlet upper quarter→stretched / root6 | 0.000227 | 0.823709 | 0.871927 |
| Triplet equilibrium→lower quarter / root5 | 0.000390 | 0.766582 | 0.873334 |
| Triplet lower quarter→midpoint / root6 | 0.000221 | 0.712833 | 0.872916 |

Root6 carries **0.550279 / 0.823482 / 0.766192 / 0.712611** squared overlap
on these links. Strong local assignment is compatible with an almost absent
projection into the five-root retained sector. Some first-five overlap blocks
also have a singular value indistinguishable from zero. This is a measured
retained-target subspace gap,independent of the direct/two-leg path disagreement.

Six roots alone do not capture the complete observed manifold. On the singlet
midpoint→upper-quarter link,roots7–12 carry **0.321949** physical squared overlap:
root7 contributes **0.193921**,root9 **0.122257**. Eight roots capture **0.744632**
of the observed twelve-root projection **0.872660**. The separate AO-following
metric reaches the same coverage conclusion; its values remain separately
recorded. These finite-manifold projections are not a claim that all missing
Hilbert-space weight has been accounted for.

## Next control

After the equilibrium CAS12 cost/export pilot and a qualified spectral comparator,
design a **fixed-orbital explicit-A2-channel comparison**. Preserve the forty-
component orbital average and vary retained A2 channels; first compare five
versus six per spin,then include the specifically measured root7/root9 coverage
in a separately budgeted comparator if the observable or projection gate requires
it. Six A2 roots per spin would mean42 total target channels; nine would mean48,
with the other six spin/irrep sectors held at five. These counts are a design
consequence,not launched native models. Verify available native target eigenpairs,
DENPROP selection/multipoles and per-model target maps before preparing decks.

Use phases and identical fixed-window position/full-width fits to test the actual
scattering consequence. Increasing only A1 averaging roots does not fill this A2
gap. Geometry subdivision is not the next control supported by this audit.
Continuity and resonance-pole identity remain unset.

## Reconstruction

Fetch the public quarter archive and read
`state-averaged-evidence/diagnostics/cas11-quarter-interval-state-tracking-20261010/result.json`.
Its SHA256 is recorded as `producer_result_sha256` in the audit. For each
`composed_path_assignments` route and each corresponding `pairs[].comparisons[]`
matrix `O=root_overlap_matrix`,take the incoming one-based row `a` and evaluate
`sum(abs(O[a-1,:k])**2)` for k=5,6,8,12. Record each omitted column separately;
check the assigned column against `assignments`,and compute singular values of
`O[:5,:5]` and `O[:6,:6]`. Keep `physical` and `ao_transport` records distinct.
The independently reconstructed public report anchors the underlying matrix
accuracy; this audit performs only projections and subspace singular values.
