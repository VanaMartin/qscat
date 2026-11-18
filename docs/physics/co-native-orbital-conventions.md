# Native CGTO point values and the saved-reader normalization trap

This is a source-resolved mathematical contract, not a validated native adapter.
It supports the [common-domain metric](co-common-domain-orbital-metric.md). The
existing orbital export contains **post-read reconstructed** shell data, not a
guaranteed literal copy of all normalization fields stored on the native file.
Its unchanged MO coefficients and stored-overlap column do not alone establish
that the exported AO point values use the coefficient generator's normalization.

## Point-value convention

For shell s, use its center A_s and Cartesian displacement X=x-A_s in bohr.
Let d_sp, v_sp, alpha_sp and K_s denote exported contractions, primitive-factor
array, exponents and shell factor. The native shell evaluator uses

    chi_s,m(x) = K_s R_l,m(X) sum_p d_sp v_sp exp(-alpha_sp |X|^2).

There is no additional radial division, spherical-Y normalization or independent
primitive renormalization in this point-value formula. The array v may already
include contraction and sphere factors; its field name is not a promise that it
contains only a bare primitive norm. Preserve the distinction between the live
post-read object and the original integral-generation object.

The solid harmonics R have R_0,0=1 and angular norm 4pi/(2l+1), not one. With
complex spherical harmonics Y_l,m including the Condon–Shortley phase,

    R_l,0 = sqrt(4pi/(2l+1)) r^l Y_l,0
    R_l,+m = sqrt(8pi/(2l+1)) (-1)^m r^l Re Y_l,m,  m>0
    R_l,-m = sqrt(8pi/(2l+1)) (-1)^m r^l Im Y_l,m,  m>0.

At the origin use their polynomial limits, not an undefined angular direction.
Within each shell, ordinal i=1,...,2l+1 means m=i-l-1. The first shells are:

| l | Values in increasing m order |
|---|---|
| 0 | 1 |
| 1 | y, z, x |
| 2 | sqrt(3)xy, sqrt(3)yz, (3z²-r²)/2, sqrt(3)xz, sqrt(3)(x²-y²)/2 |

Here x,y,z in the table are the displacement components, not absolute
coordinates. These positive p signs matter; a self Gram cannot determine them.
For a bare all-space-normalized primitive the factor multiplying R_l,m is

    N_l(alpha) = sqrt(2 (2alpha)^(l+3/2) / Gamma(l+3/2))
                 sqrt((2l+1)/(4pi)).

Do not multiply N_l into an exported v again. The construction path computes
primitive normalization, a contracted self normalization, then folds the latter
into the primitive-factor array and sets the separate shell factor to one.
Continuum normalization can subsequently fold another factor into that array.

## Reading is not a lossless normalization export

The pinned GBTO path is:

1. `molecular_basis_mod:read_basis_and_integrals` invokes the atomic basis reader.
2. `basis_data_generic_mod:read_bdg` reads each shell, then inserts it through
   `atomic_basis_mod:add_shell`.
3. Insertion calls the shell's normalization method. For a CGTO this recomputes
   bare primitive and all-space contraction factors from exponents/contractions;
   stored primitive/shell normalization fields are not simply retained.
4. Insertion folds the recomputed contraction factor into the primitive array.
5. Integral generation separately calls `normalize_continuum(a)` before
   evaluating the finite-radius integrals. The convenience **read** path does
   not itself repeat that restoration after reconstructing the atomic basis.

The source-resolved mathematical formula above therefore evaluates the exported
live object, not automatically the AO basis to which the saved C coefficients
refer. An identity-like **stored** MO overlap does not test this mismatch: that
column was read separately, not integrated from exported AO values.

This withdraws the stronger claim that the current convenience-reader export
preserves raw saved normalization factors. It does not withdraw byte decoding,
stored-Gram arithmetic or accepted resonance/bound endpoints. No numerical size
of a molecular normalization error follows from source inspection alone.

## Conditional recovery without another electronic solve

For a continuum shell the source requires A_s=0. If the original generator used
a positive radius a0, a candidate restoration factor is q_s=1/sqrt(M_s(a0)),
where M_s(a0) is the self integral of the **read-reconstructed** shell on that
sphere. For h=l+3/2, it is independent of m:

    M_s(a0) = |K_s|² 2pi/(2l+1)
       sum_pq (d_sp v_sp)* (d_sq v_sq)
       (alpha_sp+alpha_sq)^(-h)
       gamma(h, (alpha_sp+alpha_sq) a0²).

gamma here is the lower incomplete gamma. Target shells have q_s=1 for this
specific native reconstruction. This is an analytically justified **candidate**,
not an admitted patch: confirm the original generation path/radius and compare
independent absolute point values and restored self integrals before molecular
use. Alternatively, a separately admitted raw shell-record decoder can preserve
file fields without entering insertion; it still needs the generation-stage
interpretation. Neither route needs SCF/CI or a Hamiltonian eigensolve.

Let diagonal D repeat q_s over all AOs in the shell. With AO functions from the
read-reconstructed object, the candidate generator-consistent MO values use
C_eff=D C. Keep the original C and export immutable, with restoration recorded
as a separate derived stage. This is a declared basis-scale restoration, not
whitening, fitting a Gram, reorthogonalizing MOs or renormalizing sectors.
When the physical comparison domain changes, **hold a0 and D fixed**: otherwise
each radius would define new orbitals instead of testing domain sensitivity.

## Addressing, coordinates and integration domains

AO function ordinals are fixed by shell start/address tables and the m order
above. MO coefficients have AO rows and symmetry-concatenated MO columns; use
the exported mapping, not a guessed energy/irrep reorder. Point evaluation sums
AO values times the original column coefficients. No automatic molecular-axis
alignment, origin shift, orbital sign flip or phase matching is performed.

For positive native radius a0, Gaussian AO integrals start as all-space
integrals. Tails are subtracted only when **both AO shells** are marked
continuum (and therefore origin-centered); target-target and target-continuum
AO pairs retain their all-space integrals. This is a hybrid **AO** convention,
not a statement that every MO cc block is a pure sphere integral: continuum MOs
may mix target and continuum AOs.

Write G_in for the uniformly integrated ball AO Gram and G_all for all space,
both in one specified AO scale. The native hybrid AO matrix H equals G_in on
continuum-continuum AO entries and G_all on entries with at least one target AO.
The future comparison must distinguish

    stored-native prediction: C† D H D C
    uniform physical ball:    C† D G_in D C.

Their difference is the correspondingly transformed target-containing AO tail
matrix, not a normalization defect to erase. Reverse consistency and complete
many-electron self norms still use one uniform physical domain.

## Admission consequence

The unmodified export plus an assumed norm-preserving reader is **not admitted**
as a matched physical self-metric input. A next bounded stage must first validate
the read-versus-generator scale distinction with generated/analytic controls,
tag a source-justified reconstruction and predict the native-hybrid versus
uniform-domain difference. Only then consider one existing-export self-metric
sentinel, with absolute pointwise s/p/d signs, actual angular-range controls,
three-level quadrature, domain sensitivity, state-level error allocation and
hard work/memory caps. No additional native export, other-R input, molecular
integration or vector recovery is released by this source contract.

Source anchors are GBTO `gto_routines` (`eval_cgto`, `normalize_cgto`, `dngto`,
`cms_gto_norm`), `special_functions` (`cfp_solh`), `basis_data_generic`,
`atomic_basis` and `cgto_integrals` (`GG_shell_integrals`). Their exact source
hashes, caller locators and restricted-source verification belong to the
external audit record; raw native source is not redistributed here.
