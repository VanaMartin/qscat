# Fixed-QZ diffuse control: neutral completes, anion SCF fails

The original attempt stops3November2026 at14:25:19UTC. Same atomic calculation
function as original QZ,same state/core/ROHF-reference UCCSD(T) conventions and
numerical gates;only five extra diffuse shells s–g,80→105AOs. All four actual-
image/runtime controls pass,including unchanged parent-shell/overlap checks.

- Neutral O:SCF/internal stability/CCSD/triples complete,
  total **−74.99518604197657Ha**,original worker exit0.
- O−:ROHF oscillates through100cycles;SCF convergence fails,worker exit1.
  No O− CCSD/triples call. Its final SCF iterate is not an endpoint.
- Science/container/owner/launcher original exits2. **No EA or diffuse shift.**
  No rerun,optimizer fallback,gate relaxation,timeout or OOM.

Owner15.860873s,worker CPU44.236479s,kernel cumulative peak5,896,859,648bytes
on four physical cores0–3 within16GiB. Immediate terminal observation precedes
the first600s periodic check. No campaign quantum job/monitor remains active.

Independent reconstruction authenticates26 raw files and15 executed Git sources
from`2b21b43a33d94a0fb4157dc7d8a4c4a1ace4da24`,checks native neutral energies,
the original100-cycle failure and absence of anion correlation/affinity outputs.
Both the successful neutral and original failure are preserved,not requalified
as a completed diffuse test. BSE version0 input exponents and source URL/digest
are in the frozen input ledger;the public companion preserves the API response.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/oxygen-affinity-diffuse-stop
uv run --no-sync python -m projects.ukrmol_co.verify_oxygen_diffuse \
  projects/ukrmol_co/qualification-evidence/oxygen-affinity-diffuse-stop/evidence.tar.gz \
  projects/ukrmol_co/qualification-evidence/oxygen-affinity-diffuse-stop/collection.json
```

**Decision:** stop this diffuse-SCF route. Its contribution remains unmeasured;
the failure does not prove that no stationary ROHF state exists. Select one
state/core-consistent correlation diagnostic before any separate release. See
[the diagnosis](../../../../docs/physics/oxygen-affinity-diagnosis.md).
