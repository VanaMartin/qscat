# Admit pole candidates before ranking them

The potential factory's `ElectronicPair.pole` matches two ECS spectra on an
absolute energy window. Energies are Hartree. Its reported pole is the midpoint
of a nearest two-angle pair; its residual is their absolute separation.

The most angle-stable pair need not be an admissible resonance. A low-energy
continuum pair can have a smaller absolute separation than a genuine resonance
while failing the relative-width or energy-floor gate. Selecting that pair and
then returning `None` hides other admissible pairs in the same window. Recentered
Newton windows can expose this even when the bootstrap window did not.

Apply the existing absolute-residual, bound-state, positive-width,
relative-width and energy-floor gates to every nearest pair first. Return the
admitted pair with minimum residual, breaking exact ties by its original
angle-A index. No admitted pair means `None`. This changes candidate admission,
not the generic two-spectrum matcher's global-argmin contract or any tolerance.

This ranking is not unique physical state identification. Newton's separate
continuity guard remains necessary; bound/resonance bridging retains its own
rules. No candidate is accepted merely because it moves toward the target.

The midpoint's derivative and a single angle's Hellmann–Feynman derivative are
different finite-grid estimators. Correct candidate admission alone does not
establish their agreement or certify electronic resolution. Synthetic inverse
recovery, derivative consistency and convergence of the forward pole remain
distinct validation gates before molecular fitting.
