# Equilibrium five-versus-nine-A2 scattering evidence

The [matched declaration](../../a2-equilibrium-scattering-contract.json) completes
with all native,resource,comparison and owner exits0. Both original orbitals and
forty-component averaging remain unchanged. The models share the independently
qualified48-state native target/property representation,integrals/CONGEN,
audited8-rankengine and2048scattering roots. Only retained A2counts/sequences
and corresponding interface map change.

**Practical observable-insensitivity passes at equilibrium:** maximum phase
change0.0163665rad at2.55requested-eV;centers9.427–10.924meV lower and full
widths2.672–3.645meV smaller,all within0.05rad/0.05eV/5%with0.001eVfloor.
Dense→five-A2baseline phase difference5.1e-6rad passes. This is a fixed-orbital
channel control at2048,not full channel convergence or absolute omission error.
The stretched rejection and width-extraction sensitivity remain preserved.

The [local reading](../../verification-a2-equilibrium-scattering-20261016.json)
reconstructs73checksum-fetched files:full selectedCI/boundary exports,
telemetry,grids,Pi,actual decks/maps and fits. Actual large target/integral/
CONGEN inputs are hashed remotely;the four dense input files match the frozen
complete27546-pole source. The output bundle does not contain those large
integrals,CONGEN files or complete dense boundary. The separate
[target bundle](../equilibrium-a2-target/README.md) contains full48-state vectors.

Costs:five **1.65608h /6.40784GiB**,nine **1.81748h /6.89606GiB**;dense outer
replay **4.198min /0.446GiB**. Host read costs1.167seconds;host/local largest
fitted parameter difference8.293e-10eV. These readers share the same fitter and
export routines;they are independent output reconstruction,not independent
scattering engines.

## Reproduce

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/equilibrium-a2-scattering
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/equilibrium-a2-scattering"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.a2_scattering_reading \
  --directory "$REPLAY_ROOT/raw" --output "$REPLAY_ROOT/reconstructed.json"
```

Verify58source-file hashes from`reader-source-manifest.json`before replay.
The reader checks all73raw manifest hashes and recomputes both gates. Scientific
runtime is frozen at`a5954fb`;reader source is`537f5b5`. Artifact digests are
content addresses. Production,pole identity,model/basis/ensemble convergence and
width-extraction qualification remain open.
