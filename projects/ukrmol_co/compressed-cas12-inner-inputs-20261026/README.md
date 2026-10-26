# Frozen compressed CAS12 B1 inner/export cards

`scatci.B1.inp` is the **exact** selected-root/native-boundary card resolved
from the qualified compressed preparation and its own40-state property map.
2048 roots,1e-12 solver tolerance,1000 iterations,16GiB native memory input,
in-memory vectors and uncorrected250-channel native export. Both the adapter
and actual native output must match these bytes. No native solver algorithm changes.

The three outer cards are byte-identical geometry-free native recipes from
the qualified CAS11 anchor pipeline. They read their thresholds/CF/amplitudes
from the **new CAS12** fort.10/fort.21 and never supply equilibrium physical
energies, moments or target IDs. RSOLVE computes99 native K references; EIGENP
and RESON are the existing pipeline diagnostics. Native phase-fit widths are
not substituted for the later coupled outgoing pole.

All input/source pins and units are in
`../compressed-cas12-inner-contract-20261026.json`. Only B1 is released by that
contract. B2/Pi pairing, complex propagation/pole/character/model comparison,
electronic accuracy and a production sweep remain separate gates.
