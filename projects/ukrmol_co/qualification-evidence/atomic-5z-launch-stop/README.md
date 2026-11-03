# Original5Z startup-control failure is preserved

The first separately frozen5Z release stops **3November12:01:12UTC**,before any
integral/SCF/CC call. Seven runtime controls pass; calculation-equality AST hashing
fails between local Python3.12 and image Python3.11 metadata serialization. The
calculation function source is actually byte-identical to the executed QZ source.
No atom directory,result,profile or checkpoint was produced by this attempt.

Original exits are **controls1;command/container/owner/service/launcher2**,scientific
exit unset. Owner0.386853s,no timeout/OOM,source integrity passes. Eleven freshly
fetched public payloads reconstruct14 full raw files and19 original Git sources.
The original stop is not retrospectively requalified by the later correction.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-5z-launch-stop
uv run --no-sync python -m projects.ukrmol_co.verify_atomic_5z projects/ukrmol_co/qualification-evidence/atomic-5z-launch-stop/evidence.tar.gz projects/ukrmol_co/qualification-evidence/atomic-5z-launch-stop/collection.json --startup-stop
```

The [separate correction](../../atomic-5z-source-compat-decision-20261103.md)
uses exact calculation-source bytes instead of version-dependent AST metadata,
then runs all eight controls in the actual image before quantum execution. All
physical/numerical/state/timeout/memory gates remain unchanged. The
[later successful5Z calculation](../atomic-5z-endpoints/README.md) has its own
source/contract/capsule and does not repeat any successful quantum calculation.
