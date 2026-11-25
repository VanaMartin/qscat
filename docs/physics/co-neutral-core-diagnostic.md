# Matched neutral core-correlation diagnostic

`projects.ukrmol_co.neutral_core_diagnostic` runs a **diagnostic-only**
frozen-core/all-electron CCSD(T) pair on one shared stable real RHF reference.
It neither modifies nor replaces the accepted `projects.ukrmol_co.neutral`
worker. Implementation/testing is not permission to acquire molecular data;
the current CO handoff and a fresh bounded release govern physical execution.

## Observable and limits

For two explicitly admitted geometries a and b, define the energy increment

    delta_X(B) = E_X(B,b) - E_X(B,a),

where X is frozen-core (FC) or all-electron correlation (AE), and B is the
one-particle basis. All energies are Hartree and separations are bohr.
One matched pair per geometry supplies a finite-method core contrast

    C = delta_AE(W) - delta_FC(W).

Comparing its FC increment with an already accepted FC increment in another
basis V supplies a separate basis-replacement contrast

    B = delta_FC(W) - delta_FC(V).
    delta_AE(W) - delta_FC(V) = B + C.

This is an exact bookkeeping identity along that finite-method path, not an
order-independent decomposition of the physical correlation error. No fitted
offset, fragment energy or equilibrium interpolation is needed. The same
reference cancels the RHF contribution in the AE-minus-FC contrast at each
geometry. A change of basis can change reference orbitals and contractions as
well as radial flexibility. Included perturbative triples are not an estimate
of the missing higher-order correlation error.

Neither contrast establishes a complete-basis limit, relativistic correction,
physical reference uncertainty, a required nuclear accuracy allocation or a
replacement curve. No result is added to the delivered 5Z coefficients. A
literature residual is not a matched-implementation test when its Hamiltonian,
core treatment and basis differ.

## Narrow worker contract

The worker currently restricts separations to **2.45 or 2.50 bohr** and the basis
to **aug-cc-pwCVQZ**. This is deliberately not an arbitrary new sweep interface.
The versioned PySCF 2.11.0 stock basis file is checked by SHA256. Expected
dimensions are 218 spherical AOs, seven occupied orbitals and 211 virtuals.
The underlying molecule is neutral 14-electron singlet CO, with no ECP,
nonrelativistic Hamiltonian, C2v symmetry and the existing molecular origin.

1. Construct and validate the molecule/basis/runtime and parent/child source
   identities. This alone does not run SCF.
2. Converge RHF once and require both internal and external stability. No UHF
   substitution or unstable-reference exception is provided.
3. Run conventional real restricted CCSD(T) with `frozen=2`, then `frozen=None`
   on that exact reference. Check masks, occupied counts, full virtuals, finite
   amplitudes/energies, lambda convergence and a 14-electron AO-density trace.
   T1 norms are normalized by the actual correlated electron count, ten or 14.
4. Fingerprint the reference orbital coefficients, occupations, orbital energies
   and RHF total before/between/after operations. Preserve raw stage exits,
   shared-reference identity and successful partial records on failure.

SCF/CC convergence controls retain the accepted worker's tolerances. No density
fitting, virtual cutoff, extra frozen orbitals, reference restart or automatic
geometry successor is enabled. A failed AE stage does not erase the successful
FC row; a missing-only resume after failure is **not implemented** and requires
a separate decision, not a repeat of the successful peer.

## API, CLI and storage

- `preflight(config)` loads/checks the declared molecule and installed basis,
  without calling SCF/CC. It still requires the appropriate runtime-admission
  scope. Synthetic backends test this API without importing PySCF.
- `calculate_pair(config)` is the single-pair worker API. Physical use goes
  through the bounded launcher, not a direct unowned function call.
- `launch_pair(workdir, config)` creates a fresh work directory, records exact
  source identities and profiles one worker process group. It does not create
  a container, choose physical CPU topology or admit the host allocation.

For a separately approved runtime preflight:

```bash
python3 -m projects.ukrmol_co.neutral_core_diagnostic --preflight --bond-length 2.45
```

Only inside a separately released compute allocation, an operator can use:

```bash
python3 -m projects.ukrmol_co.neutral_core_diagnostic \
  --execute-reviewed --bond-length 2.45 --workdir /out/fresh-pair
```

The flag records operator intent; it is not independent scientific or resource
authority. It never releases the other geometry automatically. The worker's
source path/hash must match the parent's actual source, avoiding an image's
bundled module silently shadowing mounted source.

`diagnostic-pair.json` and `diagnostic-resources.json` use distinct diagnostic
schemas/labels, always `eligible_for_delivered_curve=false`. No `neutral.json`
is produced. A zero child exit without complete diagnostic records is rejected.
Scientific operation exits and numeric gate rejections remain distinguishable.
Raw `run.log`, profiles, shared RHF and separate FC/AE CC checkpoints belong in
the external run capsule, not Git.

The launcher monitors reference and per-core-treatment wall caps, an overall
wall cap, allocated disk and container charge. It terminates only its own child
process group and retains the original child exit separately from watchdog or
schema rejection. Polling is a sampled guard, not a mathematical bound on
transient resource use. Outer admission must independently enforce the container
memory/no-swap/CPU limits, host RAM/disk floors, overall two-geometry deadline,
image identity and terminal/OOM checks. An unavailable cgroup counter is zero,
not proof of an admitted memory allocation.

## Qualification boundary

Generated backends and inert processes test reference reuse, correct masks,
failure propagation, source/basis identity, partial retention, resource stops
and diagnostic-schema ineligibility. They do not validate PySCF's electronic
physics, the installed image/basis parser, molecular convergence, numerical
uncertainty or physical accuracy. No molecular smoke test is implied. Reusable
source must be committed/pushed before any separately admitted physical use.
