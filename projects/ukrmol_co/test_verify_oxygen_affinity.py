"""Evidence readers reject damaged/unpinned archives before any source execution."""

import pytest

from projects.ukrmol_co.verify_oxygen_affinity import verify as orbital_verify
from projects.ukrmol_co.verify_oxygen_diffuse import verify as diffuse_verify


@pytest.mark.parametrize("reader", [orbital_verify, diffuse_verify])
def test_unpinned_archive_rejected_before_unpack_or_quantum(tmp_path, reader):
    archive = tmp_path / "damaged.tar.gz"
    archive.write_bytes(b"not original evidence")
    args = (archive, {"archive": {"sha256": "0" * 64}})
    if reader is orbital_verify:
        args += (tmp_path / "missing-5z.tar.gz",)
    with pytest.raises(ValueError, match="archive"):
        reader(*args)
