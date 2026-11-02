"""Saved aug-TZ feasibility readers fail closed on packaging and changed bytes."""

import io
import tarfile

import pytest

from projects.ukrmol_co.stretched_augtz_feasibility import checked_files, supported
from projects.ukrmol_co.verify_stretched_properties import sha


def archive(entries):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as t:
        for name, data in entries:
            m = tarfile.TarInfo(name)
            m.size = len(data)
            t.addfile(m, io.BytesIO(data))
    return output.getvalue()


def test_positive_hash_pinned_saved_inventory_and_wrong_payload_reject():
    data = archive([("input", b"immutable")])
    collection = dict(
        archive=dict(bytes=len(data), sha256=sha(data)),
        files={"input": dict(bytes=9, sha256=sha(b"immutable"))},
    )
    assert checked_files(data, collection) == {"input": b"immutable"}
    for broken in [data + b"changed", data[:-1]]:
        with pytest.raises(ValueError, match="archive"):
            checked_files(broken, collection)
    collection["files"]["input"]["sha256"] = sha(b"changed")
    with pytest.raises(ValueError, match="Changed saved input"):
        checked_files(data, collection)


def test_unindexed_or_duplicate_native_records_do_not_qualify():
    data = archive([("input", b"saved"), ("input", b"saved")])
    collection = dict(archive=dict(bytes=len(data), sha256=sha(data)), files={})
    with pytest.raises(ValueError, match="duplicate"):
        checked_files(data, collection)


def test_unmaterialized_symlink_failure_remains_an_error():
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as t:
        m = tarfile.TarInfo("fort.16")
        m.type = tarfile.SYMTYPE
        m.linkname = "moints"
        t.addfile(m)
    data = output.getvalue()
    collection = dict(archive=dict(bytes=len(data), sha256=sha(data)), files={})
    with pytest.raises(ValueError, match="Nonregular"):
        checked_files(data, collection)


def test_a_different_or_relabelled_supported_proof_is_not_a_qualification():
    with pytest.raises(ValueError, match="Published supported-space proof"):
        supported(archive([("summary", b"all roots pass")]))
