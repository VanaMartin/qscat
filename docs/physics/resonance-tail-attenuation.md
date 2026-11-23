# Tail attenuation for complex resonance momenta

The real-channel ECS tail tuner uses the amplitude decay rate `K sin(theta)`.
That formula does not apply unchanged to a resonance `k=a-ib`, `a>0, b>=0`.
For a contour endpoint `Z=R0+L exp(i theta)`,

    |exp(i k Z)| = exp(b R0 - (a sin(theta)-b cos(theta)) L).

The tail exposes the outgoing resonance only if
`mu=Im(k exp(i theta))=a sin(theta)-b cos(theta)>0`. Exposure alone does not
make the finite endpoint sufficiently absorbing: real-region growth `exp(b R0)`
also consumes part of the attenuation budget. This explains how a tail can be
well order-converged yet retain substantial angle-dependent boundary reflection.

The project-level `tune_resonance_tail` adapts the existing real-channel
`qscat.tuning.tune_ecs_tail` without changing that primitive's contract. It passes
`K_eff=mu/sin(theta)` and `D_tail=D exp(-b R0)`. The returned element lengths then
give total plane-wave amplitude `|exp(i k Z)|<=D`, up to floating-point roundoff.
The adapter rejects unexposed branches and unrepresentable attenuation budgets
rather than silently using `Re(k)`, `abs(k)` or the real-energy channel momentum.

This is an a-priori plane-wave criterion, not a certified pole-energy bound.
For a p wave, the centrifugal factor and boundary-induced root sensitivity remain;
non-free tails require their own asymptotic model. Check the actual pole against
an independent outgoing oracle and a controlled tail/resolution family. A safe
potential-continuation angle must be established separately, with the stricter
double-ECS condition when both coordinates are complex. A one-point synthetic
check does not certify a molecular parameter range, nuclear grid or 2-D observable.

The adapter stays under `projects/potential_factory`; it does not change default
grids or promote a new prescription to QSCAT. A generated config is an output
requiring a separately scoped adoption decision, not an automatic replacement
for an accepted production deck.
