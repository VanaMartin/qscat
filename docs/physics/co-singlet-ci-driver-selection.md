# Selecting a singlet CI driver by spatial symmetry

`projects.ukrmol_co.target.build_target` accepts an optional
`target_singlet_irrep_drivers` mapping in its JSON configuration, for example:

```json
{"target_singlet_irrep_drivers": {"B1": "spin0", "B2": "spin0"}}
```

Only those singlet sectors use PySCF's `direct_spin0_symm`; all other sectors
retain `direct_spin1_symm`. Missing entries preserve the defaults. The legacy
`target_singlet_a1_driver` setting remains the fallback for singlet A1; an explicit
mapping entry takes precedence. Unknown irreps or driver names are rejected.
The selected driver is recorded by sector in the existing QC diagnostics.

This changes the numerical representation of a CI solve, not the intended
Hamiltonian, active space, spin, number of ensemble states or state weights.
The restricted driver uses alpha/beta-exchange-symmetric CI vectors. This admits
even-spin states, not only singlets: the spin penalty and independent physical
spin/action checks remain necessary. A converged restricted-driver flag alone
does not certify the desired physical state or lowest-root coverage.

Use full, unrestricted determinant-space Hamiltonian actions to check each
returned vector, including physical Rayleigh energies/residuals, and compare
root energies and energy-matched subspaces against a qualified reference.
Energies and action residuals are in Hartree. Preserve the orbital optimization,
fresh-CI, spin, Gram, Pi-partner and native-import gates; a fixed-orbital diagnostic
does not qualify the orbitals or authorize importing an unsuccessful checkpoint.
The prior [singlet-A1 method discussion](co-electronic-qualification.md#compressed-cas12-singlet-a1-driver-diagnosis)
describes the same distinction between a repaired CI objective and accepted QC.

Synthetic tests in `projects/ukrmol_co/test_singlet_irrep_driver.py` cover selector
isolation and invalid configurations, a five-singlet analytic spectrum in B1/B2,
and an actual mixed-spin orbital optimization compared with dense physical full
CI at the final orbitals. These controls test the implementation and its stated
small-system regimes, not electronic accuracy of the CO model. Changing a driver
in a CO calculation still requires its own frozen diagnostic rationale and release.
