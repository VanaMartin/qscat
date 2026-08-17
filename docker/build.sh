#!/usr/bin/env bash
# Build the qModeling images: a base (arch/BLAS-vendor/GPU layer) then the app
# (build/test-deps/test/runtime) FROM that base. `test-deps` installs the
# test/compute extras (mumps, plot) without running anything; `test` is FROM
# `test-deps` and adds the pytest run -- kept as a distinct, explicitly-opt-in
# target so building compute/test dependencies never silently pays for (or is
# blocked by) a full suite run. See docker/run.sh, which builds `test-deps`
# only.
#
# Usage: docker/build.sh [test|test-deps|runtime] [cpu|cpu-mkl|gpu]
#   arg1  app target   (default: test)
#     test-deps  install extras (mumps, plot) only -- no test execution
#     test       test-deps, then run the full pytest suite (~38 min)
#     runtime    the production image (no mumps)
#   arg2  base variant (default: cpu)
#     cpu      portable OpenBLAS + FFTW3 + system MUMPS (ARM/Graviton-friendly)
#     cpu-mkl  x86-64 Intel MKL variant  (SCAFFOLD — see base-cpu-mkl.Dockerfile)
#     gpu      CUDA base for future GPU kernels (SCAFFOLD — see base-gpu.Dockerfile)
set -euo pipefail

target="${1:-test}"
base="${2:-cpu}"

case "$base" in
  cpu)     base_dockerfile="docker/base.Dockerfile" ;;
  cpu-mkl) base_dockerfile="docker/base-cpu-mkl.Dockerfile" ;;
  gpu)     base_dockerfile="docker/base-gpu.Dockerfile" ;;
  *) echo "unknown base variant '$base' (expected cpu|cpu-mkl|gpu)" >&2; exit 2 ;;
esac

base_tag="qmodeling-base:${base}"

# Legacy builder avoids the slow Docker Hub BuildKit-frontend pull in this environment.
export DOCKER_BUILDKIT=0

docker build -t "$base_tag" -f "$base_dockerfile" .
# The portable CPU base also gets the :latest tag (the app Dockerfile's default).
if [ "$base" = "cpu" ]; then
  docker tag "$base_tag" qmodeling-base:latest
fi

docker build --build-arg "BASE_IMAGE=${base_tag}" \
  --build-arg "GIT_SHA=$(git rev-parse HEAD 2>/dev/null || echo unknown)" \
  --target "$target" -t "qmodeling:${target}-${base}" -f docker/Dockerfile .

echo "built ${base_tag} and qmodeling:${target}-${base}"
