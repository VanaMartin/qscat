# Saved compressed-five mapped view

**Qualified view only:** both B1/B2 retain their own2048 eigenpairs and all250
channels. Twelve channel labels per sector map from original property IDs to
dense target rows; raw labels/native bytes stay unchanged. The default reader
still raises `ValueError('Invalid channel label')`. The historical aggregate
failure and all six unset cases remain preserved.

[Result/resources](../../progress-compressed-five-view-result-20261025.json) ·
[Independent public reconstruction](../../verification-compressed-five-view-public-20261025.json) ·
[Scientific design](../../../../docs/physics/co-compressed-a2-pole-contrast.md)

Source/tests/contract freeze `9deb5a1` precedes the single saved-only application,
finished25 October22:59:00UTC. Original scientific/command/container/owner exits
are0/0/0/0, no timeout/OOM. Four physical cores0–3, one BLAS thread,4GiB /60s cap;
container2.073691s, gate.147262s. Cgroup peak.223804GiB and CPU1.486406s are sampled
at the final gate save and include preceding controls, not container teardown.
All36 runtime controls pass;39 focused tests passed locally before freeze.

Nine fresh public payloads authenticate53 archived files /39 source entries,
seven prerequisite pins and123 original inputs. Independent Fortran-record
walking, original-ID lookup and row-wise CF unpacking reproduce all six model
array hashes per sector exactly. The nine-model identity mapping is also exact
algebraically; no nine propagation is repeated. Frozen source hashes separately
match all39 Git blobs. A public-index fetch without the artifact User-Agent
returned403 before reconstruction started; standard artifact retrieval succeeded,
and the first actual reconstruction exited0. No scientific retry occurred.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-five-view
uv run --no-sync python projects/ukrmol_co/qualification-evidence/compressed-five-view/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/compressed-five-view/artifacts.json \
  verification-compressed-five-view-fresh.json
```

Replay uses NumPy and fresh public downloads; it imports no producer and runs
no controls, native executable, S construction, propagation or pole search.
Large source/result capsules live in public storage, not Git.

**Boundary:** no five-model S or pole has been evaluated. Numerical execution
still requires its exact real-guard/contrast owner and input/resource freeze.
No model-insensitivity, electronic accuracy, curve, fit or sweep verdict follows.
All114 holdout S matrices and equilibrium pole exclusion stay protected; new50
remains parked. The next deliverable is that execution freeze, not another
interface/solver development cycle or a repeat of the saved gate.
