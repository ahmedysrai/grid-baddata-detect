import pytest

from gridbad.chi2 import chi2_threshold


def test_known_value():
    # chi-square table, 13 degrees of freedom at 95%: 22.362
    assert chi2_threshold(40, 27) == pytest.approx(22.362, abs=1e-3)


def test_needs_redundancy():
    with pytest.raises(ValueError):
        chi2_threshold(27, 27)
