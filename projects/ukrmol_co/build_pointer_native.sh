#!/usr/bin/env bash
# Reuse the exact prior build recipe, selecting only the newly declared assertion patch.
set -euo pipefail
python3 - <<'PY'
from pathlib import Path
source = Path('/opt/qmodeling/projects/ukrmol_co/build_debug_native.sh').read_text()
old = 'python3 -m projects.ukrmol_co.native_debug_patch'
assert source.count(old) == 1
new = source.replace(old, 'python3 -m projects.ukrmol_co.native_pointer_patch')
Path('/home/ukrmol/debug-source/build-pointer-driver.sh').write_text(new)
PY
bash /home/ukrmol/debug-source/build-pointer-driver.sh
