"""Guard the selected-root input contract; raw engine gates are separate."""

import pytest

from projects.ukrmol_co.analyze import target_diagonalizer_properties
from projects.ukrmol_co.run import select_target_diagonalizer


def test_davidson_requires_serial_dispatch_and_preserves_output_contract(tmp_path):
    template = tmp_path / "target.scatci.inp"
    original = (
        "&input\n  megul = >>>MEGUL<<<,\n/\n&cinorn\n"
        "  nstat = >>>NSTAT<<<, nftw = >>>LUCITGT<<<, nciset = >>>NCISET<<<,\n"
        ">>>MPI_SCATCI<<<  vecstore = >>>VECSTORE<<<,\n"
        ">>>MPI_SCATCI<<<  memp = 2.5,\n>>>IGHT<<<igh = >>>IGH<<<,\n/\n"
    )
    template.write_text(original)
    select_target_diagonalizer(template, "davidson-serial", 1e-13, 800)
    selected = template.read_text()
    assert "igh = 0, forse = 1, crite = 1e-13, maxiter = 800," in selected
    assert selected.split("  igh = 0")[0] == original.split(">>>IGHT<<<")[0]
    for field in ("NSTAT", "LUCITGT", "NCISET", "VECSTORE", "MEGUL"):
        assert f">>>{field}<<<" in selected


@pytest.mark.parametrize("matches", [0, 2])
def test_reject_ambiguous_or_changed_upstream_layout(tmp_path, matches):
    template = tmp_path / "target.scatci.inp"
    original = "&cinorn\n" + ">>>IGHT<<<igh = >>>IGH<<<,\n" * matches + "/\n"
    template.write_text(original)
    with pytest.raises(ValueError, match="No unique target diagonalizer"):
        select_target_diagonalizer(template, "davidson-serial", 1e-12, 500)
    assert template.read_text() == original


def test_slepc_retains_distributed_dispatch(tmp_path):
    template = tmp_path / "target.scatci.inp"
    template.write_text("&cinorn\n>>>IGHT<<<igh = >>>IGH<<<,\n/\n")
    select_target_diagonalizer(template, "slepc", 1e-12, 500)
    assert "igh = -1, forse = 0, crite = 1e-12, maxiter = 500," in template.read_text()


@pytest.mark.parametrize("defect", [None, "dense-fallback", "roots", "incomplete", "failed"])
def test_selected_solver_gate_rejects_fallback_and_failed_iterations(tmp_path, defect):
    outputs = tmp_path / "output/CO/geom1/outputs"
    outputs.mkdir(parents=True)
    text = (
        " Sequential diagonalizations: T\nDiagonalization done with Davidson\n"
        "Requested # of eigenpairs       5\nDavidson diagonalisation completed:\n"
        " 35 iterations        39 matrix vector multiplies IERR =   0\n"
    )
    if defect == "dense-fallback":
        text = text.replace("Diagonalization done with Davidson", "ScaLAPACK chosen")
    elif defect == "roots":
        text = text.replace("eigenpairs       5", "eigenpairs       4")
    elif defect == "incomplete":
        text = text.split("Davidson diagonalisation completed:")[0]
    elif defect == "failed":
        text = text.replace("IERR =   0", "IERR =   1")
    (outputs / "target.scatci.singlet.A1.out").write_text(text)
    config = {
        "target_diagonalizer": "davidson-serial",
        "target_roots": 5,
        "target_singlet_roots": [5, 0, 0, 0],
        "target_triplet_roots": [0, 0, 0, 0],
    }
    if defect:
        with pytest.raises(ValueError, match="Incomplete or unexpected"):
            target_diagonalizer_properties(tmp_path, config)
    else:
        result = target_diagonalizer_properties(tmp_path, config)
        assert result["target_solver_sectors"] == {
            "singlet.A1": {"requested_roots": 5, "iterations": 35, "matrix_vector_multiplies": 39}
        }


@pytest.mark.parametrize("defect", [None, "dense-fallback", "incomplete", "failed"])
def test_slepc_gate_rejects_fallback_and_incomplete_output(tmp_path, defect):
    outputs = tmp_path / "output/CO/geom1/outputs"
    outputs.mkdir(parents=True)
    text = (
        "Optimized SLEPC Matrix Format chosen\nKRYLOVSCHUR used as Diagonalizer\n"
        "Requested # of eigenpairs       5\n"
        "/ Stopping condition: tol= 1.0000E-12, maxit=500, stopped at it=12\n"
        "EIGEN-ENERGIES\n"
    )
    if defect == "dense-fallback":
        text = text.replace("KRYLOVSCHUR used as Diagonalizer", "ScaLAPACK chosen")
    elif defect == "incomplete":
        text = text.replace("EIGEN-ENERGIES", "")
    elif defect == "failed":
        text += "Not all requested eigenpairs have converged!!!!"
    (outputs / "target.scatci.singlet.A1.out").write_text(text)
    config = {
        "target_diagonalizer": "slepc",
        "target_roots": 5,
        "target_singlet_roots": [5, 0, 0, 0],
        "target_triplet_roots": [0, 0, 0, 0],
    }
    if defect:
        with pytest.raises(ValueError, match="Incomplete or unexpected"):
            target_diagonalizer_properties(tmp_path, config)
    else:
        result = target_diagonalizer_properties(tmp_path, config)
        assert result["target_solver_sectors"] == {
            "singlet.A1": {"requested_roots": 5, "iterations": 12}
        }
