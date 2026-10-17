# Recovering scattering matrices from saved CO outputs

## Bounded question and interfaces

The endpoint-isolated width audit shows that local phase fidelity alone does
not identify a pole. Before implementing another extraction method, inventory
the saved native information and verify the real-axis matrices already written
by UKRmol+. This Python-stage reader does not fit phases or continue data into
the complex plane. `saved_matrices.read_k(path)` returns native Rydberg energies,
real open-open K matrices and header metadata; `scattering_matrix(k)` returns
the corresponding dimensionless S matrix. The command verifies a byte manifest
before reading the fixed twelve cases in `saved-matrices-contract.json`.

## Conventions and prospective checks

For the real symmetric K matrix in the saved flux-normalized open-channel
basis, use

\[
 S=(I-iK)^{-1}(I+iK),\qquad
 \delta_\Sigma=\sum_j\arctan\kappa_j\pmod\pi.
\]

The matrix inverse is evaluated as a linear solve. A separate spectral route
`Q diag(exp(2i arctan κ)) Qᵀ` cross-checks the solve. Analytic tests use rotated
channels with specified eigenphases; unitarity and reciprocity use the ordinary
Hermitian product, not the ECS c-product. Native eigenphase sums are printed
to limited precision: require agreement modulo π within **1e-6 rad**, while
unitarity, reciprocity, determinant-phase and spectral/Cayley consistency have
**1e-10** absolute entrywise budgets. These are representation/reproduction
checks, not the 0.001-rad prediction or 0.05-rad model-comparison budgets.

The binary layout is checked against the pinned external engine's
`source/libouter/outerio.f`, `WRITKH/WRITKM`, `READTC` and `READRH/READRM`.
The supported format is little-endian double/int64 payloads with 4-byte Fortran
record markers, one dataset, fixed nuclei, no dissociation, no Buttle/partition
correction and no saved full CI vectors in the boundary file. Unsupported or
truncated records raise; header record counts alone are not completeness proof.
Actual 99 energy records and EOF must agree. The K payload stores **only one
triangle**; unpacking it is not an independent measurement of the original
unsymmetrized K matrix's reciprocity.

The fixed inventory covers five/nine-A2 at compressed/equilibrium/stretched,
both Pi sectors, using their actual deck unit numbers. Require finite ordered
energies, a 99-point 0.1–5 requested-eV grid (`0.0735` Ryd/requested-eV), six
open channels with labels consistent with `fort.10`, 2048 boundary poles and
250/298 total channels respectively, radius 18 bohr, and two complete finite
packed multipole matrices. Check target/channel/boundary geometry and threshold
consistency; do not substitute a single-p-wave asymptotic model. All cases and
failed checks remain recorded. Failed or malformed data stop the command without
a successor calculation. This does not alter earlier owner exits or fit verdicts.

## What the files can and cannot establish

`fort.10` supplies target energies (Hartree), channel thresholds (Rydberg),
target pointers and angular labels. `fort.21` supplies total inner eigenenergies
(Hartree), boundary amplitudes and the angularly coupled long-range multipole
coefficients. The selected runs use an **uncorrected truncated pole sum**;
calling its 2048 exported poles a complete physical spectrum would be wrong.
The external source's `READRM` describes its convention as
`a R_ij(E) = Σ_k w_ik w_jk / (E_k − E)` for this no-correction layout.
An independent propagator must derive its boundary normalization/Bloch convention
from the actual writer and assembler before using that formula.

The decks and actual `rsolve` implementation also matter: the active path is
COMPAK, with coupled multipole propagation from 18 to 100 bohr and Gailitis
asymptotic matching. Saved open-open K matrices already include native closed-
channel elimination. They suffice to reconstruct **the existing real-axis S
samples**; they do not supply an independently propagated S matrix at arbitrary
real/complex energy. Nor do they uniquely specify complex continuation from a
finite set of samples. Closed-channel sheet choices, threshold square roots,
outgoing solutions, multipole normalization and radius/truncation controls are
separate physical/numerical work.

## Next gate

If the inventory and K-to-S checks pass, reuse these inputs. No new inner solve
or export is justified by a missing scalar phase background. First freeze a
physical toy scattering problem with independently known real-axis scattering
and complex poles, validate boundary matching and sheet/outgoing conventions,
then a coupled polar-tail benchmark. Only after those gates should an
independent outer propagator consume the native boundary data and reproduce the
saved full K matrices over the declared grid. Native pole extraction must wait
for that reconstruction, radius/step/tail convergence and explicit sheet/pole
stability tests. Purely rational interpolation of saved K samples is insufficient
as independent pole evidence.
