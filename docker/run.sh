#!/usr/bin/env bash
# Run a qscat-run experiment config end-to-end inside the CPU `test-deps`
# image (has MUMPS -- needed for large decks like H2P's full deck -- but does
# NOT run the test suite). Builds qmodeling-base then the `test-deps` app
# target (same pattern as docker/run-n2.sh / docker/build.sh), then runs
# `qscat-run run` inside the container, mounting CONFIG in read-only and
# OUTPUT_DIR out.
#
# Deliberately does NOT build the `test` target: that would run both test
# tiers -- including the production-scale one, measured at 13m31s on a
# 32-core x86 host and considerably longer on fewer cores -- before every
# compute invocation, and a single failing test would block compute work that
# has nothing to do with it. Run `docker/build.sh test` separately when you
# actually want the suite.
#
# Usage: docker/run.sh CONFIG [OUTPUT_DIR]
#   CONFIG      path to a qscat-run YAML config (e.g. apps/qscat-run/examples/f2-da.yaml)
#   OUTPUT_DIR  where artifacts are written (default: runs/<config-stem>)
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: docker/run.sh CONFIG [OUTPUT_DIR]" >&2
  exit 1
fi

config="$1"
config_stem="$(basename "$config")"
config_stem="${config_stem%.*}"
output_dir="${2:-runs/${config_stem}}"

mkdir -p "$output_dir"

export DOCKER_BUILDKIT=0
docker build -t qmodeling-base:latest -f docker/base.Dockerfile .
docker build --build-arg GIT_SHA="$(git rev-parse HEAD 2>/dev/null || echo unknown)" --target test-deps -t qmodeling:test-deps -f docker/Dockerfile .

# Mount the config at the SAME path it occupies in the repo, so that any
# RELATIVE path inside it (e.g. a `reference:` pointing at
# ../../../validation/n2/data/...) resolves in the container exactly as it does
# on the host. Mounting everything at /config.yaml made repo-relative configs
# resolve against / and fail. Falls back to /config.yaml for a config that
# lives outside the repo.
repo_root="$(git rev-parse --show-toplevel 2>/dev/null || echo "")"
config_abs="$(realpath "$config")"
if [[ -n "$repo_root" && "$config_abs" == "$repo_root"/* ]]; then
  in_image="/app/${config_abs#"$repo_root"/}"
else
  in_image="/config.yaml"
fi

docker run --rm \
  -v "$config_abs":"$in_image":ro \
  -v "$(realpath "$output_dir")":/out \
  qmodeling:test-deps \
  uv run --no-sync qscat-run run "$in_image" --output /out

echo "wrote artifacts to ${output_dir}"
