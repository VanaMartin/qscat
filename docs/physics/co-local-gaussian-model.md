# A finite-domain Gaussian interaction candidate for local CO fitting

The delivered local curves supply neutral V0(R) and incident poles
z(R)=Er(R)-i Gamma(R)/2, not a unique electronic interaction. A small prospective
candidate is U(r,R)=-lambda(R) exp[-alpha(R) r²], with effective ell=1 and a
separately added centrifugal term ell(ell+1)/(2r²). Keep V0's molecular reference
unchanged; the fixed-R electronic fit targets incident z without adding V0 twice.
Coordinates are bohr, energies Hartree, alpha inverse-bohr² and lambda Hartree.

`local_gaussian_model` supplies only the six-parameter forward shape. On
x=(R-2.2)/0.3, log(lambda) and log(alpha) are degree-two Chebyshev polynomials,
three real coefficients each. Positivity is built into this parameterization,
not repaired by clipping a negative width or potential parameter. Parameters and
finite real R are explicit; R is restricted to1.9–2.5bohr. Electronic r may be
complex for a qualified electronic ECS calculation. Complex nuclear R, global
extrapolation and inferred asymptotes are deliberately unsupported.

The full trial surface would be V0(R)+U(r,R), with the centrifugal term added
once in the radial electronic Hamiltonian. This module does not fit coefficients,
choose a pole, solve nuclear dynamics, infer off-shell couplings or supply a
neutral/atomic energy zero. A saved V0 evaluator and a frozen fitting experiment
must establish those inputs independently. Pole equivalence does not establish
background/channel equivalence or VE/DA accuracy.

A local calibration should use actual development geometries only, preserve all
validation exclusions, retain failed attempts and assess parameter sensitivity.
At least three electronic resolutions and two ECS angles must establish the
forward-pole numerical error separately from model residuals. Any fit is a new
finite-model candidate, never a refit of the accepted curve coefficients or a
production/global QSCAT model.

## Dimensionless scale and parameter bounds

For the mass-one fixed-R radial problem, `s=sqrt(alpha) r` gives
`H=alpha[-(1/2)d²/ds² + ell(ell+1)/(2s²) - (lambda/alpha)exp(-s²)]`.
Thus pole energies scale as `z=alpha epsilon(lambda/alpha)`. The width-to-position
ratio constrains the dimensionless depth, while the energy scale fixes alpha.
Positivity alone does not imply a particular upper bound on lambda. A declared
parameter-box escape is a failed numerical prescription, not proof that this
Gaussian family lacks a solution. Preserve that failure and separately justify
any new parameter range; it must not relax pole fidelity, resolution or validation
gates. This scaling identity supplies no physical or off-shell validation.

## Evaluating an external candidate

Candidate coefficients and qualification sidecars are immutable experiment data,
not registry constants. Find their current status and full checksums in the
[experiment catalogue](../../projects/ukrmol_co/experiments.md). Authenticate an
explicitly downloaded candidate before using it; an archived development status
is historical and must be read with its separate validation qualification.

```python
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebval

from projects.ukrmol_co.local_gaussian_model import local_gaussian_interaction

candidate = json.loads(Path("/external/qualified-candidate.json").read_text())
R = np.array([1.9, 2.2, 2.5])[None, :]   # real nuclear coordinate, bohr
r = np.array([0.1, 1.0, 5.0])[:, None]  # electronic coordinate, bohr
U = local_gaussian_interaction(candidate["coefficients"], r, R)
V0 = chebval((R - 2.2) / 0.3, candidate["neutral_coefficients_Ha"])
surface = V0 + U                       # Hartree, original molecular reference
# The ell=1 radial electronic Hamiltonian adds 1/r**2 separately, once.
```

This evaluates a finite real-R surface, not a global factory model. A validation
pass at a declared fidelity budget does not imply the tighter interpolation
accuracy of the original accepted curves. Keep stricter rejected predictions
visible, and do not tune the frozen coefficients to the validation geometries.
Once scored, those data are not unused prospective tests for a later retuned model.
Channel/background information, off-shell coupling and nuclear/scattering
validation remain separate requirements before physical predictions.

## A real-R fit domain is not a nuclear boundary condition

A fitted interval does not automatically define a boundary-insensitive nuclear
Hamiltonian. For `H_N=-(1/2mu)d²/dR²+V0(R)` on a real Dirichlet box, shrinking
the interval changes the problem, not just its discretization. Resolve grid/order
error separately before interpreting wall sensitivity. The ordinary Hermitian
norm applies to this neutral problem; an electronic Gamow-state c-product is
not its probability norm.

Inward wall moves can falsify a stated sensitivity budget, but their energy
shifts are not an upper bound on the original box's error against an unknown
larger domain. Conversely, a failed wall screen does not reject a fixed-R pole
fit. An explicitly wall-dependent finite-box diagnostic is a distinct claim from
a physical molecular vibrational level or an outgoing nuclear resonance.

Do not extrapolate the coefficients or bypass the complex-nuclear-R rejection to
make a global or nuclear-ECS solver run. Broader-domain input and a qualified
continuation need a separate decision; more resolution cannot remove a change
in the prescribed boundary conditions.
