# Published artifacts

Some results in this repository are not committed. A converged sweep is a few
hundred kilobytes of CSV that takes minutes to hours of MUMPS solves to
produce — reproducible from a committed config, but expensive enough that
recomputing it to look at a number would be absurd. Those bytes live in public
object storage; what is committed is a small pointer beside them.

## Three tiers, by who reads the number

Format follows audience, not habit:

| tier | what | format |
|---|---|---|
| **machine** | sweeps, calibration runs, packet histories — read by plotting code and comparisons, never by a person | compressed `.npz`, values **float32**, energy axis float64 |
| **published** | short tables meant to be read: resonance positions, BO levels, anion levels | CSV |
| **quoted** | numbers a note or a model states: fitted parameters, fit reports | JSON, full precision |

float32 in the machine tier is measured, not assumed: on the O₂ deck it costs a
relative 5.9e-8 on σ, against the tightest tolerance anything here is held to
(1e-3), and takes a sweep from 390 kB of CSV to 74 kB. The **energy axis stays
float64** — values tolerate rounding independently, an axis does not, because
rounding two neighbouring mesh points onto the same float turns a curve into a
multivalued one. The O₂ mesh has 134× margin at its finest spacing and loses no
points, but level-aware meshes exist to resolve peaks a few meV wide, so that
margin is not worth spending for 13 kB.

Everything a claim *depends* on stays in git: the golden inputs tests read, the
fit reports that lock model constants, the figures the physics notes discuss.
The classification, the measurements behind it, and why the line is drawn there
are recorded in `docs/adr/0008-computed-artifacts-live-in-public-object-storage.md`
— a repository-only document, so it is named by path rather than linked: `adr/`
is excluded from this site.

## Fetching

```console
$ qscat-run fetch validation/factory/results/o2-ve
validation/factory/results/o2-ve: fetched 3 file(s)
```

No account, no credentials, no extra package — reads are anonymous HTTPS. To
see what a directory would pull without pulling it:

```console
$ qscat-run fetch validation/factory/results/o2-ve --list
validation/factory/results/o2-ve  (from 69742d8)
  https://data.qscat.org/o2-ve/cross_section.830cffb8a044.csv
  https://data.qscat.org/o2-ve/cross_section.927f86a1ff2f.npz
  https://data.qscat.org/o2-ve/cross_section.f8808b6e6355.png
```

Only outputs are on that list. Every published run directory also carries a
committed `config.resolved.yaml` and `manifest.json` — the input the sweep is a
function of, and the record of what produced it — so a clone with no network
can still say exactly what a published number came from and how to re-run it.
Both are kilobytes; putting them behind a download would make them unavailable
precisely when they are wanted.

Every file is checked against the sha256 recorded when it was published. A
mismatch is an error, and the bad bytes are not written — a truncated download
must never end up on disk looking like data. A file already present and correct
is skipped, so re-running costs nothing and an interrupted fetch resumes.

Directories without an `artifacts.json` keep their results in git and need no
fetching.

## Referencing an artifact directly

The URL is stable, so anything that can read a URL can use one — a plotting
notebook, a script in another language, a `curl` in a shell, a colleague who
has never cloned this repository:

```console
$ curl -O https://data.qscat.org/o2-ve/cross_section.<digest>.npz
```

```python
import io, urllib.request, numpy as np
url = "https://data.qscat.org/o2-ve/cross_section.<digest>.npz"
d = np.load(io.BytesIO(urllib.request.urlopen(url).read()))
E, sigma = d["energy"], d["ti:ve:v0->0"]
```

Take the digest from the run's `artifacts.json`, or from
`https://data.qscat.org/<experiment>/index.json`.

The path is the addressing scheme:

```
https://data.qscat.org/<experiment>/<name>.<sha256[:12]>.<ext>
```

The folder is a readable label; the **digest is the address**. That is what
makes such a URL safe to paste into a paper, an issue or a message: different
content hashes differently and therefore lives somewhere else, so a link
cannot quietly come to mean something other than what you cited. If a result
is corrected, the correction gets its own URL and the old one keeps the old
value for as long as it is kept, so the two remain comparable.

It also means a re-run that reproduces its numbers republishes to the same
address and changes nothing — which is the usual outcome here.

There is deliberately no "latest" alias; that would be the one URL whose
meaning could change. To find the current bytes for an experiment, read
`artifacts.json` in the repository at the commit you care about — or fetch
`https://data.qscat.org/<experiment>/index.json`, the one object addressed by
name rather than by content, which is why it may be replaced.

Nothing expires today. Blobs are shared between pointers once two runs produce
the same bytes, so deleting by age is unsafe; the eventual mechanism is
reachability — dropping blobs no pointer names — rather than a clock.

## Publishing

Publishing is maintainer-only and lives in the private `qscat-infra`
repository, which holds the Terraform for the bucket and the upload tool. It
needs an R2 token with write access to the artifacts bucket.

This is not obscurity: the bucket is bound to `data.qscat.org` as a read-only
hostname, and every write goes through a separate authenticated S3 endpoint.
There is no public write path to close.
