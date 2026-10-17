# Native CO surface amplitudes, energy units and boundary matching

## Scope and source chain

This is a source-level derivation for the saved **nine-A2, fixed-nuclei,
uncorrected 2048-pole, zero-Bloch** exports. It resolves how to assemble the
native boundary matrix. It is not a comparison with independently propagated
native K or proof of the missing-spectrum error. Source files are preserved
with their licences and SHA256 manifests in the combined-exterior evidence.

Read the actual pinned engine, not a generic R-matrix convention:

1. `gbtolib/source/gto_routines.f90::CGTO_amplitude`, lines 1001–1018:
   Gaussian primitive contraction is multiplied by `r**(l+1)`, yielding the
   **reduced radial** surface value. The harmonic normalization is included.
2. `atomic_basis_mod.f90::calculate_amplitudes`, lines 3926–4032, and
   `molecular_basis_mod.f90::calculate_amplitudes`, lines 4103–4155: AO surface
   values are contracted to MOs. `normalize_to_a=.false.` avoids renormalizing
   continuum orbitals at a newly requested evaluation radius; it does not
   divide the surface value by a.
3. `ukrmol_interface.f90::EVAL_AMPLITUDES`, lines 1039–1074, and
   `::UKP_READAMP`, lines 1108–1205: **line 1176 divides by sqrt(2)**. Its
   explanatory comment explicitly says the factor absorbs the half in the
   Hartree R-matrix formula. This is upstream of the CI contraction.
4. `Postprocessing_module.F90::get_boundary_data`, lines 899–969:
   `eval_amplitudes(r,.false.)`, `ukp_readamp`, then orbital-amplitude/CI-vector
   multiplication. The later **sqrt(2) multiplication is guarded by
   `ir<nfdm+1`**, so it applies only to interior RMT sampling points.
5. `::write_boundary_data`, lines 1243–1255: write total-Hartree poles
   `ei+core_energy` and the **surface** `wamp(nfdm+1)` in channel-fast, pole-slow
   Fortran order. No radius or extra half is introduced here.
6. `libouter/outerio.f::READRM` entry, lines 124–209: read poles and amplitudes
   unchanged. Its comment at line 154 names `RMATR*RMATRX`; the active assembler
   and matching algebra below determine which dimensional matrix is actually
   passed. Do not add a radius division merely from that comment.
7. `compak/rsolve_compak.f::RESIDR_compak`, lines 1527–1532 and 1565–1578:
   compute `WAMP/(EPOLE−ETOT)` and multiply by WAMP, **without another half**.
   `SQUARM_compak`, lines 2149–2166, only unpacks the symmetric triangle.
8. `::KMAT_compak`, lines 2204–2229: form `F−R(F'−BSTO F)`. Thus the active
   matrix maps a **radial derivative** to a value and carries length units.

## Derived matching equation

For reduced radial eigenfunctions normalized in the inner-region Hilbert space,
write physical surface overlaps as γ_in. Native stored amplitudes are
`w_in=γ_in/sqrt(2)`. The Hartree spectral Green identity gives

\[
 Q_{ij}(E)=\frac12\sum_n\frac{\gamma_{in}\gamma_{jn}}{E_n-E}
          =\sum_n\frac{w_{in}w_{jn}}{E_n-E},\qquad
 u(a)=Q(E)[u'(a)-\beta u(a)].
\]

Q has units bohr; the commonly dimensionless R is **Q/a**, paired with
`a u'−b u` and b=aβ. For the saved decks **β=BBLOCH=0**: the keyword is absent
and the native default is zero (`rsolve_compak.f`, line 164). A Python consumer
must either use `(Q,u')` or `(Q/a,a u')` consistently. It must not mix Q/a
with u', insert an additional half, or undo the already absorbed sqrt(2).
Nonzero native BBLOCH input semantics across all propagation paths are not
qualified by these zero-Bloch decks.

At complex energy use **transpose products**, not Hermitian conjugation:
Q(E) is complex symmetric and satisfies `Q(E*)=Q(E)*` for these real amplitudes.
Its energy derivative is `Σ wwᵀ/(E_n−E)²`, which is positive semidefinite on
real intervals away from poles. None of this makes a truncated pole set complete.
For a qualified outgoing exterior map Y, the matching singularity is
`det[I−Q(E)(Y(E)−βI)]=0`; this formula is not permission to search CO poles
before native matrix reconstruction and spectral/model controls pass.

## Energy zero and outer potential

`rsolve_compak.f`, lines 247–251 converts requested eV with **0.0735 Ryd/eV**;
line 368 selects e0=first target total Hartree energy; line 611 computes
`ETOT=e0+ENRYD/2`. `RESIDR` therefore receives total **Hartree** energy, not
incident Rydberg energy. Channel thresholds are `2(E_target−e0)` and the outer
radial k² is `ENRYD−threshold_Ryd` (`ASYM2_compak`, lines 1104–1110).
Keep the historical requested-eV conversion when reproducing those files;
it differs from conversion with the modern Hartree/eV constant by about
**1.84445×10⁻⁵ relatively**.

`Postprocessing_module.F90`, lines 1408–1420 already multiplies multipole
properties by **2** to convert the radial operator to Rydberg units.
`POTL_compak`, lines 2132–2142 then uses

\[
 u''=\left[\operatorname{diag}(\Delta_i-E_{\rm inc,Ryd})
 +\operatorname{diag}(l_i(l_i+1))/r^2-2Z I/r
 +\sum_{\lambda=1}^{2}CF_\lambda/r^{\lambda+1}\right]u.
\]

The saved CF matrices are already twice the Hartree potential coefficients;
**do not multiply them by two again**. They contain angular couplings and are
packed row-wise lower triangles per multipole order. Neutral CO has Z=0.
The two orders are inverse-square/inverse-cube terms, not evidence of a complete
inverse-fourth-power polarization model.

## Independent analytic normalization check

A free radial box, u(0)=0,u'(a)=0, has
`p_n=(n+1/2)π/a`, `E_n=p_n²/2`, normalized `γ_n=(-1)^n sqrt(2/a)`.
Therefore native `w_n=(-1)^n/sqrt(a)` and the exact sum is

`Σ w_n²/(E_n−E)=tan(ka)/k`, k²=2E.

This directly reproduces u(a)/u'(a), while the dimensionless convention yields
`tan(ka)/(ka)`. At finite N the omitted sum falls as `2a/(π²N)`; three-rung
tests must observe that truncation rather than claiming finite-N equality.
Multiplying Q by half/two or dividing by a gives persistent discrepancies.
This mathematical check anchors the derived formula independently of CO phases.

## Qualification boundary

The source trace resolves surface sqrt(2), total-Hartree denominators, the
length-valued zero-Bloch matrix, and the already converted CF coefficients.
Native full-K reproduction, Gailitis incoming/outgoing normalization and general
nonzero-Bloch propagation remain separate evidence requirements. The combined
exterior benchmark's rejected refinement gate is preserved. Five-A2 label
readings remain unset; no mapped decoder, engine correction, native rerun or
CO complex-pole result is implied here.
