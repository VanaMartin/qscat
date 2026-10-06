#!/usr/bin/env bash
# Compile the checksum-pinned UKRmol-in/out engine and bundled GBTOlib.
set -euo pipefail

root=/home/ukrmol/source
mkdir -p "$root"
cd "$root"
curl -fL --retry 3 --max-time 180 \
  'https://zenodo.org/records/18340340/files/ukrmol-in-3.3.0.tar.gz?download=1' \
  -o ukrmol-in-3.3.0.tar.gz
curl -fL --retry 3 --max-time 180 \
  'https://zenodo.org/records/18538198/files/ukrmol-out-3.3.0.1.tar.gz?download=1' \
  -o ukrmol-out-3.3.0.1.tar.gz
sha256sum --check <<'CHECKSUMS'
1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2  ukrmol-in-3.3.0.tar.gz
6fecd651d1b953004926bcd2485d18a1645f2a13fca93e176458b0d48f556a9b  ukrmol-out-3.3.0.1.tar.gz
CHECKSUMS
mkdir inner outer
tar -xzf ukrmol-in-3.3.0.tar.gz --strip-components=1 -C inner
tar -xzf ukrmol-out-3.3.0.1.tar.gz --strip-components=1 -C outer
test -f inner/source/gbtolib/CMakeLists.txt

cmake -S inner -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER="$(command -v mpifort)" \
  -DCMAKE_Fortran_FLAGS='-fdefault-integer-8' \
  -DBLAS_LIBRARIES=/opt/ukrmolp/lib/libopenblas.so \
  -DLAPACK_LIBRARIES=/opt/ukrmolp/lib/libopenblas.so \
  -DSCALAPACK_LIBRARIES=/opt/ukrmolp/lib/libscalapack.so \
  -DUKRMOL_OUT_DIR="$root/outer" \
  -DCMAKE_INSTALL_PREFIX="$root/install" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  '-DCMAKE_INSTALL_RPATH=/opt/ukrmolp/lib.double;/opt/ukrmolp/lib' \
  -DMPIEXEC_MAX_NUMPROCS=2 \
  -DWITH_MPI=ON -DWITH_GIT=OFF -DBUILD_DOC=OFF \
  -DBUILD_TESTING=ON -DWITH_PSI4=ON -DWITH_NUMDIFF=ON
cmake --build build -j"${1:-4}"
cmake --install build

python3 - <<'PY'
import hashlib
import json
import subprocess
from pathlib import Path

root = Path('/home/ukrmol/source')
record = {
    'ukrmol_in': '3.3.0', 'ukrmol_out': '3.3.0.1',
    'gbtolib': 'Sources bundled inside the pinned UKRmol-in archive',
    'precision': 'double', 'integer_bytes': 8,
    'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob('*.tar.gz')},
    'compiler': subprocess.check_output(['gfortran', '--version'], text=True).splitlines()[0],
    'cmake': subprocess.check_output(['cmake', '--version'], text=True).splitlines()[0],
    'toolchain_image': 'zdenekmasin/ukrmol_plus@sha256:1cef2ee6aeaca2a4d5b813ff31e7d28bef4436a4282398c8933a1c047a8ec76c',
}
(root / 'provenance.json').write_text(json.dumps(record, indent=2) + '\n')
PY
