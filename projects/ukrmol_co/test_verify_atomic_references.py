"""Raw convergence/energy parsing must not infer a pass from a producer success flag."""

import pytest

from projects.ukrmol_co.verify_atomic_references import log_energies

LOG = """converged SCF energy = -37.6
<class 'ROHF'> wavefunction is stable in the internal stability analysis
UCCSD converged
E(UCCSD) = -37.7  E_corr = -.1
UCCSD(T) correction = -0.02
"""


def test_known_raw_log_energy_sum_and_convergence():
    assert log_energies(LOG)["total_energy_hartree"] == pytest.approx(-37.72, abs=1e-12, rel=0)


@pytest.mark.parametrize(
    "text",
    [
        LOG.replace("UCCSD converged", "UCCSD not converged"),
        LOG.replace("stable in the internal", "unstable in the internal"),
        LOG + "UCCSD(T) correction = -0.03\n",
        LOG.replace("converged SCF energy", "SCF failed"),
    ],
)
def test_missing_failed_or_duplicate_native_records_reject(text):
    with pytest.raises(ValueError):
        log_energies(text)
