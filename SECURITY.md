# Security Policy

## Supported versions

| Version | Supported |
|---|---|
| `main` (this repository) | yes |

qscat is not on PyPI; the only distribution is this repository, and only
the tip of `main` is supported. There are no maintained release branches.

## Reporting a vulnerability

Report privately through GitHub's advisory form:
<https://github.com/VanaMartin/qscat/security/advisories/new>
(repository **Security** tab → "Report a vulnerability"). Please do not
open a public issue for a suspected vulnerability. You can expect an
acknowledgement within one week.

Scope worth knowing: `qscat-run` executes YAML configuration files that
select solvers and write artifacts to paths named in the config. Treat
configs from untrusted sources as untrusted input.
