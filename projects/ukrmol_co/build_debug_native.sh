#!/usr/bin/env bash
# Persistent diagnostic build over the retained source and unchanged numerical libraries.
set -euo pipefail
root=/home/ukrmol/debug-source
mkdir "$root/inner" "$root/outer"
cd "$root"
sha256sum --check <<'CHECKSUMS'
1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2  /work/upstream/ukrmol-in-3.3.0.tar.gz
6fecd651d1b953004926bcd2485d18a1645f2a13fca93e176458b0d48f556a9b  /work/upstream/ukrmol-out-3.3.0.1.tar.gz
CHECKSUMS
tar -xzf /work/upstream/ukrmol-in-3.3.0.tar.gz --strip-components=1 -C inner
tar -xzf /work/upstream/ukrmol-out-3.3.0.1.tar.gz --strip-components=1 -C outer
python3 -m projects.ukrmol_co.native_debug_patch --source "$root/inner" --output "$root/patch.json"
export OMPI_FC=/opt/compiler/bin/gfortran OMPI_CC=/opt/compiler/bin/gcc
cmake -S inner -B build \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_Fortran_COMPILER="$(command -v mpifort)" \
  -DCMAKE_Fortran_FLAGS='-fdefault-integer-8' \
  -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO='-O1 -g -fbacktrace -fcheck=bounds' \
  -DBLAS_LIBRARIES=/opt/ukrmolp/lib/libopenblas.so \
  -DLAPACK_LIBRARIES=/opt/ukrmolp/lib/libopenblas.so \
  -DSCALAPACK_LIBRARIES=/opt/ukrmolp/lib/libscalapack.so \
  '-DSLEPC_LIBRARIES=/opt/ukrmolp/lib/libslepc.so;/opt/ukrmolp/lib/libpetsc.so' \
  '-DSLEPC_INCLUDE_DIRS=/opt/ukrmolp/include;/opt/ukrmolp/include/petsc;/opt/ukrmolp/include/slepc' \
  -DUKRMOL_OUT_DIR="$root/outer" \
  -DCMAKE_INSTALL_PREFIX="$root/install" -DCMAKE_INSTALL_LIBDIR=lib \
  '-DCMAKE_INSTALL_RPATH=/opt/ukrmolp/lib.double;/opt/ukrmolp/lib' \
  -DMPIEXEC_MAX_NUMPROCS=2 \
  -DWITH_MPI=ON -DWITH_GIT=OFF -DBUILD_DOC=OFF \
  -DBUILD_TESTING=ON -DWITH_PSI4=ON -DWITH_NUMDIFF=ON
cmake --build build -j4
cmake --install build
python3 - <<'PY'
import hashlib,json,subprocess
from pathlib import Path
root=Path('/home/ukrmol/debug-source')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
binary=root/'install/bin/mpi-scatci'
sections=subprocess.check_output(['readelf','-S',str(binary)],text=True)
assert '.debug_info' in sections and '.debug_line' in sections
record={
 'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
 'patch':json.loads((root/'patch.json').read_text()),
 'cmake_cache_sha256':sha(root/'build/CMakeCache.txt'),
 'debug_info_and_line_sections_present':True,
 'installed_sha256':{str(p.relative_to(root/'install')):sha(p) for p in sorted((root/'install').rglob('*')) if p.is_file()},
 'numerical_library_sha256':{str(p):sha(p) for p in [Path('/opt/ukrmolp/lib')/n for n in ('libslepc.so','libpetsc.so','libopenblas.so','libscalapack.so')]},
 'packages':subprocess.check_output(['dpkg-query','-W'],text=True),
}
(root/'build-provenance.json').write_text(json.dumps(record,indent=2)+'\n')
PY
