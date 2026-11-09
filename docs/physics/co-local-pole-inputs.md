# Boundary-only inputs for a local CO pole node

An outgoing pole needs the finite native spectrum, boundary amplitudes,
thresholds, angular labels and exterior multipoles—not the saved real-energy
K grid. A boundary-only loader avoids implicitly decoding protected scattering
samples. It changes no propagation, root, contour or residue algorithm.

`local_pole_inputs.load_boundary_model(source)` accepts only the declared neutral
48-target/298-channel/2048-root, two-multipole, 18-bohr reference. Thresholds
must match individual target totals; every closed channel is retained. `c2`
includes centrifugal l(l+1), and `c3` has no extra Rydberg factor. Native total
and ground energies are Hartree; thresholds/exterior coefficients follow the
Rydberg radial equation. Binary geometry fields are zero in these exports;
the archived target/deck provenance establishes R and electronic prescription.

Pole dispatch accepts stretched/nine/B1 or B2 only with the explicit
`local-cas11-nine-2048-v1` view. Its owner supplies hash-pinned successful
same-case real guards at development indices 16 and 20, not reserved every-fifth
samples. Existing compressed dispatch/guards are unchanged. The local view
authenticates four boundary/channel files, not unrelated raw scattering grids.

Synthetic controls check threshold conversion, closed-channel retention, packed
symmetric coefficients, layout/radius rejection and refusal of absent views,
wrong cases or inadequate guards. Every scientific row still needs unchanged
pole/residue/character checks and its own bounded external release. No multi-R
identity, omitted-spectrum or electronic-accuracy claim follows from loading.
