# Compressed B1 finite-reference outgoing pole —23 October 2026

## Result and production consequence

**One numerically witnessed outgoing pole now exists for the saved compressed
CAS11/nine-A2/B1 finite reference**, R=1.9bohr, all298 boundary channels and
2048 exported inner eigenpairs:

| Quantity | Value |
|---|---:|
|Incident pole /Hartree|0.1381917721280088 −0.0383040828514810i|
|Position above native ground threshold /physical eV|**3.7604**|
|Full width Gamma=−2 Im z /physical eV|**2.0846**|

Raw reporting values are3.760389687392806 /2.0846143865399505eV. Extra digits
preserve provenance, not physical uncertainty. Native total energy uses the
saved ground target plus this incident Hartree value. Reporting uses canonical
`qscat.units`, not the historical0.0735Ryd/requested-eV grid conversion.

This is the first pole/residue datum for this **declared finite reference**,
not a qualified CO width curve, an electronically converged result, a scalar
potential, a production model choice or a full-scattering embedding. The
[target-map limitations](../../docs/physics/co-effective-resonance-target.md)
remain. B1/B2 consistency, other geometries and model uncertainty are still open;
do not feed one B1 pole into the factory as a qualified CO target curve.
B1 by itself is not a Lambda assignment: paired-sector degeneracy alone also
would not distinguish Pi from higher angular character in the same irreps.
Historical phase-fit plots retain their own conventions; subtracting those
candidates from this pole is not an isolated, matched electronic-model contrast.

## Frozen experiment and witnesses

[Design](../../docs/physics/co-compressed-pole-reference.md),
[contract](compressed-pole-reference-contract-20261023.json), source and tests
freeze/push at **`9b4f588`**, before CO root use. Search remains inside the
prescribed box with the original seed pair; **47/77 allowed maps** are used:
10 boundary-root evaluations, one complex Jost center, four derivative points
and32 contour points. No alternate seed, box expansion, retry or repair occurs.
All eleven controls pass in the pinned runtime. Original scientific/command/
container/owner exits **0/0/0/0**; no timeout/OOM.

- **Boundary root:** the first setting takes eight secant evaluations. Both
  prescribed tighter settings accept the same starting candidate in one
  evaluation each, with bordered |alpha|≈3e-14 and second-smallest singular
  value≈0.170. They do **not** independently iterate to new root estimates.
  Identical reported positions mean the starting candidate passes each residual
  gate, not that root uncertainty is zero. The same-root position/full-width
  budget remains0.1meV; the finite setting family is not an absolute error enclosure.
- **Interior resolvent:** nested16/32-node contour moments recover the same pole.
  Raw contour-to-root position/full-width differences are2.16e-13 /1.23e-13eV,
  with nonzero projected residue norm0.01801775 and rank ratio5.67e-16.
  These tiny construction agreements are **not** absolute CO error bounds.
- **Isolation:** discrete determinant winding is+1 for both contour rules;
  maximum phase jumps0.543 /0.272rad meet the frozen pi/2 gate. This is a
  finite resolution/isolation test, not a rigorous global pole count.
- **Visible scattering residue:** the6x6 open-channel residue is nonzero,
  Frobenius norm0.05689072Ha, rank ratio5.20e-15. It is not just a hidden
  interior-state pole. Its16→32 relative change is1.35e-12. No off-axis
  real/imaginary split or use of the real-axis S formula is involved.
- **Derivative normalization:** fixed1e-5 /5e-6Ha steps give residue differences
  from the contour of2.54e-7 /6.33e-8 relative, below1e-3. The center null and
  bilinear residue symmetry gates pass separately.

Bordered matching, interior resolvent moments and open-channel residue are
different constructions, but **share the exterior ODE and native D/W**. Small
controls have independent analytic/dense/BVP oracles; this is not a second
298-channel complex exterior engine. No omitted-spectrum error bound follows
from2048 roots or the earlier real-axis2048→4096 phase control. The accepted
economical baseline is retained, not automatically enlarged.

## Cost, preservation and independent saved reconstruction

Run finishes **14:45:31UTC**, container wall **180.990698s**, Python180.326268s.
Parent CPU110.490766s excludes contour workers; final cgroup CPU capture
**660.464411s** includes the worker/control processes. Cgroup peak
2,086,596,608bytes (**1.94329GiB**); parent peak338,056KiB (**.32240GiB**) is a
separate measure. Eight physical cores0–7, up to eight contour workers, one BLAS
thread each,8GiB allocation. All123 original inputs, four prerequisite pins and
25 frozen source entries verify before/after. No canary or pole job remains.

The earlier canary's original dependency stop2/2/2 is preserved; this experiment
does not relabel it or the older propagation/extraction failures. There are
**zero new engine, QC, property or inner solves, extra geometries, potential fits
or holdout S analyses**. The114 route holdout S values remain unanalyzed and the
equilibrium pole exclusion stays protected. No `ResonanceTarget` is exported.

[Fresh public saved reconstruction](verification-compressed-pole-reference-public-20261023.json)
verifies **7 payloads /86 archived files /25 source entries**, downloads three
prerequisite payloads and authenticates all123 original inputs. It reconstructs
Q for all47 authenticated map captures, ten bordered values,32 contour points
and both derivative witnesses without importing producer numerical modules or
propagating an ODE.
It uses independently decoded native records/Q contraction, opposite-order
**I−YQ** Schur/residue solves, NumPy bordered solves and einsum contour moments.
Maximum Q drift8.97e-15; relative interior/transition drifts6.23e-13 /6.42e-13,
moment drift8.71e-14 and derivative-residue drift3.07e-11. Every saved pole/
residue gate reconstructs, with original exits retained. Shared saved Y are
authenticated, **not independently repropagated**. Reproduction accuracy does
not turn the finite-model result into absolute CO accuracy.

[Fetch/reproduce](qualification-evidence/compressed-pole-reference/README.md).
Complete raw points, Q/Y/A matrices, projections, moments, root traces, source
capsule, owner, logs and resources persist on Sadaharu at
`prepared/compressed-pole-reference-source-20261023/`, with public mirror
`prepared/publication-compressed-pole-reference-20261023/`.

## Goal checkpoint and one next decision

**Evidence delta:** the pole-target route now has one coupled, open-channel-visible
finite-reference datum, rather than only real-axis fit parameters or a non-pole
complex map probe. This removes a local pole-identification blocker. It does
not select the electronic model, certify omitted spectrum, establish continuity
or deliver a multi-geometry width curve/potential.

**Disposition: continue with a separately bounded paired-sector/character consistency
proposal at this same compressed anchor**, before extending to other geometries
or fitting a curve. Test the B2 counterpart at the same-root budgets and the
open-channel angular character; failure would stop the proposed single-Pi target
assignment. Degeneracy alone is not that assignment.
Success should return effort to a decision-changing matched electronic-model
contrast, not another general exterior/precision campaign. No second sector,
anchor, root ladder or geometry job is launched by this outcome. Model, neutral,
continuity, factory training-node minimum and held-out physical scattering remain
separate gates; new50 stays parked. Review on that proposal or the next active
working-cycle boundary.
