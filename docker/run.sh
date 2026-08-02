#!/usr/bin/env bash
# Run a qscat-run experiment config end-to-end inside the CPU `test` image
# (has MUMPS -- needed for large decks like H2P's full deck). Builds
# qmodeling-base then the `test` app target (same pattern as
# docker/run-n2.sh / docker/build.sh), then runs `qscat-run run` inside the
# container, mounting CONFIG in read-only and OUTPUT_DIR out.
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
docker build --target test -t qmodeling:test -f docker/Dockerfile .

docker run --rm \
  -v "$(realpath "$config")":/config.yaml:ro \
  -v "$(realpath "$output_dir")":/out \
  qmodeling:test \
  uv run --no-sync qscat-run run /config.yaml --output /out

echo "wrote artifacts to ${output_dir}"
