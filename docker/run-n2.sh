#!/usr/bin/env bash
# Run the N2 LCP benchmark harness inside the CPU runtime image.
set -euo pipefail
export DOCKER_BUILDKIT=0
docker build -t qmodeling-base:latest -f docker/base.Dockerfile . >/dev/null
docker build --target runtime -t qmodeling:runtime -f docker/Dockerfile . >/dev/null
docker run --rm qmodeling:runtime python validation/n2/experiment.py
