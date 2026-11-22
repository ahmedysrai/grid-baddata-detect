import copy

import numpy as np
import pandapower.networks as pn
import pytest

from gridbad.attack import (
    ac_attack,
    apply_attack,
    linear_attack,
    meters_touched,
    state_position,
)
from gridbad.chi2 import chi2_threshold
from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import normalised_residuals
from gridbad.wls import estimate

CASE14 = pn.case14()


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def angle_shift(res, bus, degrees):
    c = np.zeros(res.n_states)
    c[state_position(res, bus, "va")] = np.radians(degrees)
    return c


def test_ac_attack_moves_noise_free_estimate_by_c(case14_pf):
    meas = full_measurement_set(case14_pf)
    res = run(meas)
    c = angle_shift(res, 13, 10)
    c[state_position(res, 8, "vm")] = 0.02
    attacked = run(apply_attack(meas, res, ac_attack(res, c), case14_pf.sn_mva))
    np.testing.assert_allclose(attacked.x, res.x + c, atol=1e-8)
    assert attacked.J < 1e-10


def test_ac_attack_leaves_residual_tests_unchanged(case14_pf, rng):
    meas = add_noise(full_measurement_set(case14_pf), rng)
    res = run(meas)
    c = angle_shift(res, 13, 20)
    attacked = run(apply_attack(meas, res, ac_attack(res, c), case14_pf.sn_mva))
    # r is the same at x_hat + c, but H differs there, so J moves slightly.
    # stealthy_attack.py (seed 0) saw 58.0 -> 57.9 at 20 degrees.
    assert attacked.J == pytest.approx(res.J, rel=0.01)
    assert np.nanmax(normalised_residuals(attacked)) == pytest.approx(
        np.nanmax(normalised_residuals(res)), rel=0.02
    )
    j = state_position(res, 13, "va")
    assert attacked.x[j] - res.x[j] == pytest.approx(c[j], rel=0.01)


def test_linear_attack_is_caught_once_the_shift_is_large(case14_pf, rng):
    meas = add_noise(full_measurement_set(case14_pf), rng)
    res = run(meas)
    thr = chi2_threshold(res.n_meas, res.n_states)
    small = angle_shift(res, 13, 1)
    large = angle_shift(res, 13, 20)
    a_small = linear_attack(res, small)
    a_large = linear_attack(res, large)
    assert run(apply_attack(meas, res, a_small, case14_pf.sn_mva)).J < thr
    assert run(apply_attack(meas, res, a_large, case14_pf.sn_mva)).J > 10 * thr


def test_bus_14_angle_touches_only_its_own_meters(case14_pf):
    meas = full_measurement_set(case14_pf)
    res = run(meas)
    a = ac_attack(res, angle_shift(res, 13, 5))
    touched = meas.loc[res.meas_index[meters_touched(a)]]
    lines = case14_pf.line
    to_14 = set(lines.index[(lines.from_bus == 13) | (lines.to_bus == 13)])
    # P and Q injection at bus 14 and its neighbours 9 and 13, and P and Q
    # on the two lines into bus 14
    assert len(touched) == 10
    assert set(touched[touched.element_type == "bus"].element) == {8, 12, 13}
    assert set(touched[touched.element_type == "line"].element) == to_14
    assert set(touched.meas_type) == {"p", "q"}


def test_slack_angle_is_not_a_state(case14_pf):
    res = run(full_measurement_set(case14_pf))
    with pytest.raises(ValueError):
        state_position(res, 0, "va")
    assert state_position(res, 0, "vm") == 13
