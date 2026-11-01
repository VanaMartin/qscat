"""Independent native framing/finite-array/triangular known truths; no scientific calls."""

import io
import struct
import tarfile

import numpy as np
import pytest

from projects.ukrmol_co import verify_stretched_inner as reader


def record(data):
    return struct.pack("<I", len(data)) + data + struct.pack("<I", len(data))


def test_record_framing_and_complete_eof():
    data = record(b"first") + record(b"second")
    assert list(map(bytes, reader.records(data))) == [b"first", b"second"]
    for broken in [data[:-1], data + b"x", record(b"first")[:-4] + b"\0" * 4]:
        with pytest.raises(ValueError):
            reader.records(broken)


def test_native_array_size_and_nonfinite_values_are_rejected():
    a = np.asarray([1.0, 2.0, 3.0], dtype="<f8")
    np.testing.assert_allclose(reader.numbers(a.tobytes(), "<f8", 3), a, atol=0, rtol=0)
    for data, count in [
        (a.tobytes(), 2),
        (np.asarray([np.nan]).tobytes(), 1),
        (np.asarray([np.inf]).tobytes(), 1),
    ]:
        with pytest.raises(ValueError):
            reader.numbers(data, "<f8", count)


def test_packed_cf_and_k_symmetric_arithmetic_known_truth():
    packed = np.asarray([1, 2, 3, 4, 5, 6.0])
    expected = np.asarray([[1, 2, 4], [2, 3, 5], [4, 5, 6.0]])
    np.testing.assert_allclose(reader.symmetric(packed, 3), expected, atol=0, rtol=0)
    batch = reader.symmetric(np.stack([packed, packed + 10]), 3)
    np.testing.assert_allclose(batch, np.stack([expected, expected + 10]), atol=0, rtol=0)
    with pytest.raises(ValueError):
        reader.symmetric(packed[:-1], 3)


def test_stream_retains_telemetry_jsonl_and_requested_native_record(tmp_path):
    path = tmp_path / "evidence.tar.gz"
    payloads = {"eps-telemetry.jsonl": b'{"requested":2048}\n', "fort.10": record(b"header")}
    with tarfile.open(path, "w:gz") as archive:
        for name, data in payloads.items():
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
    collection = dict(
        archive=dict(bytes=path.stat().st_size, sha256=reader.sha(path.read_bytes())),
        files={n: dict(bytes=len(data), sha256=reader.sha(data)) for n, data in payloads.items()},
    )
    kept, _, _ = reader.stream_files(path, collection, capture_names=["fort.10"])
    assert kept == payloads


def test_post_csv_final_disk_measurement_excludes_its_later_summary_file():
    allocation = dict(
        files={
            "native/fort.21": dict(allocated_bytes=1024),
            "resources.csv": dict(allocated_bytes=512),
            "resources.json": dict(allocated_bytes=512),
        }
    )
    assert reader.disk_peak(1024, allocation) == (1536, 1536)
    assert reader.disk_peak(2048, allocation) == (2048, 1536)


@pytest.mark.parametrize("header_length", [208, 336])
def test_ci_header_extension_keeps_base_fields_and_exact_vector_records(header_length):
    header = bytearray(header_length)
    header[16:136] = b"CO CI calculation - doublet B1".ljust(120)
    struct.pack_into("<4q", header, 136, 3, 87198, 2048, 1)
    struct.pack_into("<ddq", header, 168, 0.5, 0.5, 15)
    struct.pack_into("<d", header, 192, 1.23)
    energies = np.arange(2048, dtype="<f8") * 0.01
    mixed = (
        np.ones(87198, dtype="<i8").tobytes()
        + energies.tobytes()
        + np.zeros(87198, dtype="<i8").tobytes()
    )
    data = record(b"********                CIDATA  ") + record(header)
    for center in [(0.0, 0.0, 0.0, 6.0), (0.0, 0.0, 2.5, 8.0), (0.0, 0.0, 0.0, 0.0)]:
        data += record(struct.pack("<q4d", 0, *center))
    data += record(mixed)
    zero = np.zeros(2256, dtype="<f8").tobytes()
    data += b"".join(record(struct.pack("<q", n) + zero) for n in range(1, 2049))
    np.testing.assert_allclose(reader.ci_export(data), energies + 1.23, atol=1e-12, rtol=0)
    with pytest.raises(ValueError):
        reader.ci_export(data[:-1])
