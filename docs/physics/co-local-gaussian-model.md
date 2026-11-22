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
