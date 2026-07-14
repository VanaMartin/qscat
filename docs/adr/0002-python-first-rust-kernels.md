# 2. Python-first development with Rust kernels

Date: 2026-07-12
Status: accepted

## Context
Prior work (eMoScat) is C++/CUDA, judged increasingly costly to maintain. Research
velocity favors Python; heavy numerics need a compiled language. GPU (libXcuda) is
recovered but deferred.

## Decision
Python (uv) is the primary language. Validated methods live in `qscat`. Proven hot
paths move to Rust (PyO3/maturin) with the Python version kept as the differential
oracle. Everything runs on CPU locally and is containerizable. GPU/AWS deferred.

## Consequences
Fast iteration; safe optimization; reproducible builds. Two toolchains (Python + Rust)
must be maintained. See the `qm-method-lifecycle` skill.
