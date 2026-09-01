# Published artifacts

Some results in this repository are not committed. A converged sweep is a few
hundred kilobytes of CSV that takes minutes to hours of MUMPS solves to
produce — reproducible from a committed config, but expensive enough that
recomputing it to look at a number would be absurd. Those bytes live in public
object storage; what is committed is a small pointer beside them.

Everything a claim *depends* on stays in git: the golden inputs tests read, the
fit reports that lock model constants, the figures the physics notes discuss.
The classification, the measurements behind it, and why the line is drawn there
are recorded in `docs/adr/0008-computed-artifacts-live-in-public-object-storage.md`
— a repository-only document, so it is named by path rather than linked: `adr/`
is excluded from this site.

## Fetching

```console
$ qscat-run fetch validation/factory/results/o2-ve
validation/factory/results/o2-ve: fetched 4 file(s)
```

No account, no credentials, no extra package — reads are anonymous HTTPS. To
see what a directory would pull without pulling it:

```console
$ qscat-run fetch validation/factory/results/o2-ve --list
validation/factory/results/o2-ve  (from 69742d8)
  https://data.qscat.org/o2-ve/config.resolved.e3a2d2aa3203.yaml
  https://data.qscat.org/o2-ve/cross_section.830cffb8a044.csv
  https://data.qscat.org/o2-ve/cross_section.927f86a1ff2f.npz
  https://data.qscat.org/o2-ve/cross_section.f8808b6e6355.png
  https://data.qscat.org/o2-ve/manifest.edcf292d07f6.json
```

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
$ curl -O https://data.qscat.org/o2-ve/cross_section.830cffb8a044.csv
```

```python
import pandas as pd
df = pd.read_csv("https://data.qscat.org/o2-ve/cross_section.830cffb8a044.csv")
```

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
