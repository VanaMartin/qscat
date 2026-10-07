# Original staged QC and fifty-component ensemble diagnostics

The original CAS(10,12) staged batch retains exits **[1,1,1,1,0]**. The first
four failures remain numerical rejections:

| Original attempt | Rejection | Wall / seconds |
|---|---|---:|
| R=1.9 DZ | Singlet A1 CI | 3551.82 |
| R=2.5 DZ | Singlet A2 and triplet B1/B2 CI | 7795.03 |
| Equilibrium projected TZ | Orbital gradient 1.05004e-5 and singlet A2 CI | 29192.74 |
| Equilibrium projected aug-TZ | Triplet A2 CI | 8676.26 |

The fifth entry is an equilibrium DZ **fifty-component** QC pilot, averaging
`[10,5,5,5]` roots per spin. Its CAS(10,11) counterpart also passes QC. Compare
like root sets: the table uses the **common forty roots evaluated at the new
orbitals**, rather than comparing the two differently weighted raw objectives.

| Diagnostic, fifty versus forty components | CAS(10,11) | CAS(10,12) |
|---|---:|---:|
| Reported ground energy change / eV | +0.456512 | +0.007930 |
| Ground z dipole change / a.u. | +0.034687 | +0.003604 |
| Maximum common-root energy change / eV | 4.281167 | 0.052079 |
| Common-forty objective change / eV | +0.170147 | +0.004204 |
| Wall / seconds | 3164.60 | 2450.85 |
| Peak container / GiB | 0.3679 | 0.8784 |

These are QC/ensemble diagnostics. Independent lowest-root coverage and native
all-root/dipole import remain unset; the large CAS11 root changes require those
checks before any scattering interpretation. Both initial fifty-component
checkpoints hash-identically to their corresponding tight forty-component
equilibrium checkpoint. The projected-TZ rejection remains excluded from
automatic retry.

## Fetch and reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/staged-targets
mkdir -p /tmp/co-staged-targets
tar -xzf projects/ukrmol_co/qualification-evidence/staged-targets/staged-targets.tar.gz \
  -C /tmp/co-staged-targets
uv run python /tmp/co-staged-targets/state-averaged-evidence/verify-staged-targets.py \
  /tmp/co-staged-targets/state-averaged-evidence
```

The verifier checks **545 payloads / 281438 profile samples**, completed owner
records and supervisor/manifest hashes, original batch exits, failure flags and
gradients, passing QC spin/orbital/Molden/ensemble consistency, initial checkpoint
identity and common-root diagnostics. Both executable target versions travel
with the data: the earlier source is resolved byte-exactly from commit `5fe2bda`,
and the later producer snapshot remains frozen with its original hash.
Original logs, profiles, diagnostics, checkpoints and staged source/commands
are retained. Sorted normalized tar/gzip repackaging is byte-identical.
