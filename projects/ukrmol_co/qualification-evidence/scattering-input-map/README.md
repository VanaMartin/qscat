# Saved scattering phase and property-ID map

[Outcome](../../SCATTERING_INPUT_MAP.md):317,988 phase-defining first determinants
reproduce saved B1/B2 target phases exactly. New lowest-five property selection
needs **25 changed idtarg entries**; energies agree with original CIDATA within
4.91639e-10Ha. Original679/new676 moment files remain unchanged. No native call.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/scattering-input-map
python projects/ukrmol_co/qualification-evidence/scattering-input-map/reconstruct-saved.py \
  --directory projects/ukrmol_co/qualification-evidence/scattering-input-map
```

Six public payloads include the results archive (1,477,010 bytes, SHA256
`ec8a60e1a2b5b437dd333d007142dbd1673bb6c6e42f1acad8b5583629e86dca`),
collection hashes, manifest, summary, a prior CIDATA archive pointer and an
independent stdlib saved-witness reconstruction. The exact archive hash is also
in `artifacts.json` and `collection-sha256.json`.

The50 archived files preserve both owners, freezes, original exits and partial
outputs, source manifests, cards/property files, compact first-term witnesses,
target-only phase slices and raw scattering header/phase records. Initial freeze
`1899892` fails1/1/1 in card parsing; parser-only successor `07a64e8` passes0/0/0.
No failed directory is overwritten. All19 original inputs plus4 predecessor pins
verify during the successful remote audit.

Public replay uses **permutation-cycle parity**, independently of the producer's
inversion-count formula. It authenticates all6 payloads /50 files and a freshly
downloaded prior archive containing original CIDATA, reproduces317,988 saved
phases and reconstructs property IDs/absolute energies. No native program or
producer reader is called. A fresh `verification-saved.json` is created; an
existing output is not overwritten.

**Limit:** the original1.4GB N+1 CSFs are NOT republished. Their complete hashes
were checked remotely; derived compact slices do not let public replay verify
those original hashes or all determinant terms. Do not call this a full raw-CSF
reconstruction or native scattering-consumer qualification. Persistent originals
remain on Sadaharu at the exact paths/pins in the archived contract. Actual
native scattering-destination reads, downstream property-consumer semantics and
the scientific contrast remain gated. Both physical A1 rejections persist.
