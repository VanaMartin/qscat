# Stretched aug-TZ supported-space coverage

All **48 fresh probes** at trial spaces **80/160/240** pass at the unchanged
CAS12 aug-TZ checkpoint: two spins, four irreps and 5/8 roots, 600 cycles and
the original residual/spin/energy/orthogonality gates.

| Independent reconstruction | Maximum |
|---|---:|
| Ensemble energy error / Hartree | 1.990e-13 |
| Physical residual / Hartree | 7.084e-10 |
| Spin-penalized residual / Hartree | 9.986e-10 |
| Trial-space energy error / Hartree | 2.274e-13 |
| Root-count energy error / Hartree | 1.421e-13 |
| Spin-square error | 8.438e-15 |

Cost is **871.86 seconds / 0.8383 GiB**, with **269 payloads / 4346 profiles**
independently reconstructed. The production CI trial space80 is covered by this
passing convergence family. The original space40 scan remains rejected.

The combined owner still exits **one** because an independent outer wrapper
rejects I_XSECS's valid completion text; the supported-space scan and verifier
each exit **zero**. The archive also preserves the earlier serial-under-MPI
stdin-EOF launch and the completion-wrapper failure. Overall controller failure
does not erase the independently completed coverage result.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/augtz-supported-spaces
tar -xzf projects/ukrmol_co/qualification-evidence/augtz-supported-spaces/augtz-supported-spaces.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/augtz-supported-spaces
PYTHONPATH=projects/ukrmol_co/qualification-evidence/augtz-supported-spaces/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/augtz-supported-spaces/state-averaged-evidence/verify-augtz-supported-bundle.py \
  projects/ukrmol_co/qualification-evidence/augtz-supported-spaces/state-averaged-evidence
```

Public fetch/reconstruction/byte-identical repackaging pass. The
[separate native-import release](../../cas12-augtz-supported-import-contract.json)
uses this proof without changing the original rejected scan or queued import.
