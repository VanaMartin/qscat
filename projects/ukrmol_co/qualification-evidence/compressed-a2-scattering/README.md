# Compressed five-versus-nine-A2 scattering evidence

**The compressed comparison passes its declared observable budgets.** Both
sectors give a maximum phase change of **0.0115944 rad at 3.9 requested-eV**,
centers **13.241–18.210 meV lower**, and full widths **2.218–4.404 meV smaller**.
All four background orders pass in the prospectively declared 2–5 requested-eV
window. Budgets remain 0.05 rad / 0.05 eV / 5% with a 0.001-eV floor.

The [declaration](../../a2-compressed-scattering-contract.json) preserves the
original compressed CAS11/cc-pVDZ orbitals and forty-component averaging. Both
models share the qualified 48-state target/property representation, integrals,
CONGEN, audited eight-rank engine and 2048 scattering roots. Only retained A2
counts/sequences and corresponding interface map change. Owner completion is
16 October at 16:55:18.733 UTC; all native, resource, reader and owner exits are
zero. The observer ends normally after nine captures.

The [independent reading](../../verification-a2-compressed-scattering-20261017.json)
reconstructs **73 checksum-fetched files**: complete selected CI/boundary exports,
telemetry, grids, Pi symmetry, actual decks/maps and both matched fit gates.
Actual large shared native input equality is checksum-audited remotely. Four
dense input hashes match the frozen complete 27546-pole source. The output
archive omits the large integrals/CONGEN and complete dense boundary; the
[target bundle](../compressed-a2-target/README.md) holds all 48 target vectors.

Five/nine costs are **1.651263 / 1.827292 hours**, **6.011147 / 7.245762 GiB**;
dense outer replay costs **4.248 minutes / 0.446 GiB**. Dense→five baseline phase
difference is 2.2e-6 rad. Both selected models have zero Pi phase split and
native residuals at most 1.37e-10 Hartree. Host replay takes 0.815 seconds;
largest host/local fitted-parameter difference is 4.018e-9 eV. These shared
readers/fitter reconstruct outputs independently of the owner, but are not an
independent scattering engine.

## Reproduce

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-a2-scattering
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/compressed-a2-scattering"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.a2_scattering_reading \
  --directory "$REPLAY_ROOT/raw" --output "$REPLAY_ROOT/reconstructed.json"
```

Verify 58 source-file hashes from `reader-source-manifest.json` before replay.
The reader verifies all 73 raw manifest hashes and recomputes both gates.
Scientific runtime is frozen at `89a5275`; unchanged reader source is `537f5b5`.
Digests are content addresses. Equilibrium also passes, while stretched rejects
after its named refinement. Full channel convergence, absolute omission, model/
basis/ensemble convergence, global pole identity and width extraction remain
unqualified. The original compressed 24-case extraction rejection is preserved.
