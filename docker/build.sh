#!/usr/bin/env bash
# Build the qModeling images in order: base (arch/vendor) then app (build/test/runtime).
# Usage: docker/build.sh [test|runtime]   (default: test)
set -euo pipefail
target="${1:-test}"
# Legacy builder avoids the slow Docker Hub BuildKit-frontend pull in this environment.
export DOCKER_BUILDKIT=0
docker build -t qmodeling-base:latest -f docker/base.Dockerfile .
docker build --target "$target" -t "qmodeling:${target}" -f docker/Dockerfile .
echo "built qmodeling-base:latest and qmodeling:${target}"
