# Published artifacts

Some results in this repository are not committed. A converged sweep is a few
hundred kilobytes of CSV that takes minutes to hours of MUMPS solves to
produce — reproducible from a committed config, but expensive enough that
recomputing it to look at a number would be absurd. Those bytes live in public
object storage; what is committed is a small pointer beside them.

Everything a claim *depends* on stays in git: the golden inputs tests read, the
fit reports that lock model constants, the figures the physics notes discuss.
See {doc}`adr/0008-computed-artifacts-live-in-public-object-storage` for the
classification and why it is drawn there.

## Fetching

```console
$ qscat-run fetch validation/factory/results/o2-ve
validation/factory/results/o2-ve: fetched 4 file(s)
```

No account, no credentials, no extra package — reads are anonymous HTTPS. To
see what a directory would pull without pulling it:

```console
$ qscat-run fetch validation/factory/results/o2-ve --list
validation/factory/results/o2-ve  (from c884f51)
  https://data.qscat.org/main/c884f51/o2-ve/config.resolved.yaml
  https://data.qscat.org/main/c884f51/o2-ve/cross_section.csv
  https://data.qscat.org/main/c884f51/o2-ve/cross_section.npz
  https://data.qscat.org/main/c884f51/o2-ve/manifest.json
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
$ curl -O https://data.qscat.org/main/c884f51/o2-ve/cross_section.csv
```

```python
import pandas as pd
df = pd.read_csv("https://data.qscat.org/main/c884f51/o2-ve/cross_section.csv")
```

The path is the addressing scheme:

```
https://data.qscat.org/<scope>/<sha7>/<experiment>/<file>
                       main | branch | tag
```

`<sha7>` is the commit the run was made from, so a URL names both the numbers
and the code that produced them. **Published paths are immutable**: a key is
never overwritten. If a result is later corrected, the correction is published
under the new commit's path and the old URL keeps the old values — so a note
that cites a number keeps citing the number it was written about, and you can
always compare the two.

That immutability is why these URLs are safe to paste into a paper, an issue,
or a message. What it does not give you is a "latest" alias; there is
deliberately no such thing. Follow the pointer in the repository at the commit
you care about.

Artifacts under `main/` and `tag/` are permanent. Artifacts under `branch/`
expire after 90 days — fine for sharing a work-in-progress number, not for
citing one.

## Publishing

Publishing is maintainer-only and lives in the private `qscat-infra`
repository, which holds the Terraform for the bucket and the upload tool. It
needs an R2 token with write access to the artifacts bucket.

This is not obscurity: the bucket is bound to `data.qscat.org` as a read-only
hostname, and every write goes through a separate authenticated S3 endpoint.
There is no public write path to close.
