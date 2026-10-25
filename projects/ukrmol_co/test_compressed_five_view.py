"""Synthetic native records: explicit target maps never change stored channels/W/D."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.io import FortranFile

from projects.ukrmol_co.compressed_five_view import load_compressed_five, selected_ids
from projects.ukrmol_co.native_full_reconstruction import _model_from_interface, load_model
from projects.ukrmol_co.saved_matrices import read_interface

SELECTION = [*range(1, 39), 40, 42]


def write_interface(source, ids, labels, energies, channels, *, roots=4, amplitude_value=0.02):
    """Known native record layout; no production writer or physics solver."""
    source.mkdir(parents=True, exist_ok=True)
    header = [np.array([1, 2, 2, 0], dtype="<i8"), np.array([0, 1], dtype="<f8")]
    with FortranFile(source / "fort.10", "w", header_dtype="<u4") as stream:
        stream.write_record(np.array([10, 1, 0, 0, 0], dtype="<i8"))
        stream.write_record(np.array([b"synthetic".ljust(80)], dtype="S80"))
        stream.write_record(np.array([len(ids), 0, 0, len(channels)], dtype="<i8"))
        stream.write_record(*header)
        for i, (label, energy) in enumerate(zip(labels, energies, strict=True), 1):
            stream.write_record(np.array([i, *label], dtype="<i8"), np.array([energy], dtype="<f8"))
        lookup = dict(zip(ids, energies, strict=True))
        for i, (target, angular, projection) in enumerate(channels, 1):
            stream.write_record(
                np.array([i, target, angular, projection], dtype="<i8"),
                np.array([2 * (lookup[target] - energies[0])], dtype="<f8"),
            )
    cf = np.zeros((2, len(channels) * (len(channels) + 1) // 2))
    cf[0, 0] = 0.03
    cf[1, 0] = -0.04
    cf[0, 1] = 0.12
    cf[1, 1] = 0.17
    poles = -99.7 + amplitude_value + np.arange(roots) * 0.002
    w = np.full((roots, len(channels)), amplitude_value)
    with FortranFile(source / "fort.21", "w", header_dtype="<u4") as stream:
        stream.write_record(np.array([11, 1, 0, 0, 0], dtype="<i8"))
        stream.write_record(np.array([b"synthetic".ljust(80)], dtype="S80"))
        stream.write_record(np.array([len(ids), 0, 0, len(channels)], dtype="<i8"))
        stream.write_record(*header)
        stream.write_record(np.array([2, roots, 0, 0], dtype="<i8"), np.array([18.0], dtype="<f8"))
        stream.write_record(cf.astype("<f8").ravel())
        stream.write_record(poles.astype("<f8"))
        stream.write_record(w.astype("<f8").ravel())
    return poles, w


def small(source, ids=(1, 3, 7)):
    return write_interface(
        source,
        ids,
        [(0, 1, 0)] * 3,
        [-10, -9.9, -9.7],
        [(ids[0], 1, 1), (ids[1], 3, 1), (ids[2], 3, 3)],
    )


def write_k(source, targets, channels):
    (source / "rsolve.inp").write_text("&rslvin lukmt=921, calcdip=0,1,0, /\n")
    with FortranFile(source / "fort.921", "w", header_dtype="<u4") as stream:
        stream.write_record(np.array([11, 1, 0, 0, 0], dtype="<i8"))
        stream.write_record(np.array([b"synthetic K".ljust(80)], dtype="S80"))
        stream.write_record(np.array([1, 2, 2, 0], dtype="<i8"), np.array([0, 1], dtype="<f8"))
        stream.write_record(np.array([targets, 0, 0, channels, 1], dtype="<i8"))
        stream.write_record(
            np.array([1, 99], dtype="<i8"), np.array([0.00735, 0.003675], dtype="<f8")
        )
        for index in range(99):
            stream.write_record(
                np.array([6, 0, 21], dtype="<i8"),
                np.array([(0.1 + 0.05 * index) * 0.0735, *([0.0] * 21)], dtype="<f8"),
            )


def full_fixture(root):
    five = root / "raw/compressed/five/B1"
    nine = root / "raw/compressed/nine/B1"
    groups = [(irrep, spin, 0) for spin in [1, 3] for irrep in range(4) for _ in range(5)]
    by_id = dict(zip(SELECTION, groups, strict=True))
    by_id.update(
        {
            original: (3, 1 if i < 4 else 3, 0)
            for i, original in enumerate([39, 41, 43, 44, 45, 46, 47, 48])
        }
    )
    expected = {}
    for source, ids, count, value in [
        (five, SELECTION, 250, 0.02),
        (nine, list(range(1, 49)), 298, 0.07),
    ]:
        opened = [
            (1, angular, projection)
            for angular, projection in [(1, 1), (2, 1), (3, 1), (3, 3), (4, 1), (4, 3)]
        ]
        closed = [(ids[1 + i % (len(ids) - 1)], 1, 1) for i in range(count - 6)]
        expected[source.name + str(count)] = write_interface(
            source,
            ids,
            [by_id[i] for i in ids],
            [-100 + 0.3 * (i - 1) for i in ids],
            opened + closed,
            roots=2048,
            amplitude_value=value,
        )
        write_k(source, len(ids), count)
    order = SELECTION[::2] + SELECTION[1::2]
    (five.parent / "pilot.json").write_text(json.dumps(dict(model="five", target_map=order)))
    deck = root / "writer-evidence/compressed/five/B1/scatci.inp"
    deck.parent.mkdir(parents=True)
    deck.write_text(
        "&OUTER_INTERFACE ntarg=40, idtarg=" + ",".join(map(str, order)) + ", rmatr=18, /\n"
    )
    return five, nine, deck, expected["B1250"]


def test_explicit_mapping_preserves_raw_ids_thresholds_and_bytes():
    with TemporaryDirectory() as temporary:
        source = Path(temporary)
        small(source)
        before = {p.name: p.read_bytes() for p in source.iterdir()}
        with unittest.TestCase().assertRaisesRegex(ValueError, "Invalid channel label"):
            read_interface(source / "fort.10", source / "fort.21")
        actual = read_interface(
            source / "fort.10", source / "fort.21", selected_property_ids=[7, 1, 3]
        )
        assert [row[1] for row in actual["channel_labels"]] == [1, 3, 7]
        assert actual["channel_target_rows"] == [1, 2, 3]
        assert actual["selected_property_ids"] == [1, 3, 7]
        np.testing.assert_allclose(
            actual["channel_thresholds_rydberg"], [0, 0.2, 0.6], atol=1e-14, rtol=0
        )
        assert before == {p.name: p.read_bytes() for p in source.iterdir()}


def test_identity_mapping_matches_strict_metadata():
    with TemporaryDirectory() as temporary:
        source = Path(temporary)
        small(source, (1, 2, 3))
        strict = read_interface(source / "fort.10", source / "fort.21")
        mapped = read_interface(
            source / "fort.10", source / "fort.21", selected_property_ids=[3, 1, 2]
        )
        assert {key: mapped[key] for key in strict} == strict


def test_duplicate_missing_float_and_wrong_in_range_ids_fail():
    with TemporaryDirectory() as temporary:
        source = Path(temporary)
        small(source)
        for ids in [[1, 3, 3], [1, 3], [1.0, 3.0, 7.0], [1, 2, 7], [0, 3, 7]]:
            with unittest.TestCase().assertRaises(ValueError):
                read_interface(source / "fort.10", source / "fort.21", selected_property_ids=ids)


def test_threshold_corruption_and_clipped_ids_fail():
    with TemporaryDirectory() as temporary:
        source = Path(temporary)
        small(source)
        from projects.ukrmol_co.saved_matrices import records

        rows = records(source / "fort.10")
        label = np.frombuffer(rows[-1][:32], "<i8").copy()
        label[1] = 3  # Clipping the out-of-range original7 to dense3 is wrong.
        rows[-1] = label.tobytes() + rows[-1][32:]
        with FortranFile(source / "fort.10", "w", header_dtype="<u4") as stream:
            for row in rows:
                stream.write_record(np.frombuffer(row, "u1"))
        with unittest.TestCase().assertRaisesRegex(ValueError, "threshold"):
            read_interface(source / "fort.10", source / "fort.21", selected_property_ids=[1, 3, 7])


def test_corrupt_threshold_and_nonconsecutive_written_target_rows_fail():
    from projects.ukrmol_co.saved_matrices import records

    for mode in ["threshold", "row_order"]:
        with TemporaryDirectory() as temporary:
            source = Path(temporary)
            small(source)
            rows = records(source / "fort.10")
            if mode == "threshold":
                energy = np.frombuffer(rows[-1][32:], "<f8").copy() + 1e-4
                rows[-1] = rows[-1][:32] + energy.tobytes()
            else:
                label = np.frombuffer(rows[5][:32], "<i8").copy()
                label[0] = 3
                rows[5] = label.tobytes() + rows[5][32:]
            with FortranFile(source / "fort.10", "w", header_dtype="<u4") as stream:
                for row in rows:
                    stream.write_record(np.frombuffer(row, "u1"))
            with unittest.TestCase().assertRaises(ValueError):
                read_interface(
                    source / "fort.10", source / "fort.21", selected_property_ids=[1, 3, 7]
                )


def test_full_view_uses_own_eigenpairs_and_shared_assembly_not_nine_slice():
    with TemporaryDirectory() as temporary:
        five, nine, _, (poles, w) = full_fixture(Path(temporary))
        actual = load_compressed_five(five)
        np.testing.assert_array_equal(actual["poles"], poles)
        np.testing.assert_array_equal(actual["amplitudes"], w)
        assert actual["amplitudes"].shape == (2048, 250)
        for label, dense in zip(
            actual["interface"]["channel_labels"],
            actual["interface"]["channel_target_rows"],
            strict=True,
        ):
            if label[1] == 42:
                assert dense == 40
            if label[1] == 40:
                assert dense == 39
        np.testing.assert_allclose(actual["c2"][0, 0], 2.03, atol=1e-15, rtol=0)
        np.testing.assert_allclose(actual["c3"][0, 0], -0.04, atol=1e-15, rtol=0)
        np.testing.assert_allclose(actual["c2"][[0, 1], [1, 0]], [0.12, 0.12], atol=1e-15, rtol=0)
        np.testing.assert_allclose(actual["c3"][[0, 1], [1, 0]], [0.17, 0.17], atol=1e-15, rtol=0)
        strict = load_model(nine)
        mapped = _model_from_interface(
            nine,
            read_interface(
                nine / "fort.10", nine / "fort.21", selected_property_ids=list(range(1, 49))
            ),
            (48, 298),
        )
        for key in strict:
            if isinstance(strict[key], np.ndarray):
                np.testing.assert_array_equal(strict[key], mapped[key])
            else:
                assert strict[key] == mapped[key]
        assert not np.array_equal(actual["amplitudes"], strict["amplitudes"][:, :250])
        assert not np.array_equal(actual["poles"], strict["poles"])


def test_pilot_deck_selection_and_full_reference_metadata_must_agree():
    with TemporaryDirectory() as temporary:
        five, nine, deck, _ = full_fixture(Path(temporary))
        original = deck.read_text()
        deck.write_text(original.replace("ntarg=40", "ntarg=39"))
        with unittest.TestCase().assertRaisesRegex(ValueError, "selection"):
            selected_ids(five)
        deck.write_text(original)
        text = (five.parent / "pilot.json").read_text()
        pilot = json.loads(text)
        pilot["target_map"][-1] = 41
        (five.parent / "pilot.json").write_text(json.dumps(pilot))
        with unittest.TestCase().assertRaisesRegex(ValueError, "selection"):
            load_compressed_five(five)
        (five.parent / "pilot.json").write_text(text)
        from projects.ukrmol_co.saved_matrices import records

        rows = records(nine / "fort.10")
        energy = np.frombuffer(rows[4][32:], "<f8").copy() + 1e-4
        rows[4] = rows[4][:32] + energy.tobytes()
        with FortranFile(nine / "fort.10", "w", header_dtype="<u4") as stream:
            for row in rows:
                stream.write_record(np.frombuffer(row, "u1"))
        with unittest.TestCase().assertRaisesRegex(ValueError, "subset"):
            load_compressed_five(five)


def test_scope_cannot_select_equilibrium_or_another_channel_model():
    with unittest.TestCase().assertRaisesRegex(ValueError, "compressed/five"):
        load_compressed_five(Path("raw/equilibrium/five/B1"))


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
