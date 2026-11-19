import copy

import numpy as np
import pandapower.networks as pn
import pytest
from pandapower.estimation.state_estimation import StateEstimation

from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import lnr_identify, normalised_residuals, sensitivity_diagonal
from gridbad.wls import estimate

CASE14 = pn.case14()


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def row(meas, meas_type, element_type, element):
    mask = (
        (meas.meas_type == meas_type)
        & (meas.element_type == element_type)
        & (meas.element == element)
    )
    return meas.index[mask][0]


def test_sensitivity_trace_is_degrees_of_freedom(case14_pf, rng):
    res = run(add_noise(full_measurement_set(case14_pf), rng))
    s_ii = sensitivity_diagonal(res)
    assert s_ii.sum() == pytest.approx(res.n_meas - res.n_states, abs=1e-8)
    assert np.all((s_ii > 0) & (s_ii < 1))


def test_meas_index_follows_dropped_rows(case14_pf):
    meas = full_measurement_set(case14_pf).drop([5, 40])
    res = run(meas)
    assert sorted(res.meas_index) == sorted(meas.index)


def test_no_redundancy_means_every_measurement_is_critical(case14_pf, rng):
    # V everywhere plus P injection at the 13 non-slack buses: m = n = 27
    meas = full_measurement_set(case14_pf)
    keep = (meas.meas_type == "v") | (
        (meas.meas_type == "p") & (meas.element_type == "bus") & (meas.element != 0)
    )
    meas = add_noise(meas[keep], rng)
    res = run(meas)
    assert res.n_meas == res.n_states
    assert np.isnan(normalised_residuals(res)).all()
    removed, _ = lnr_identify(copy.deepcopy(CASE14), meas)
    assert removed == []


@pytest.mark.parametrize(
    "target", [("p", "bus", 3), ("q", "line", 2), ("v", "bus", 13)]
)
def test_single_error_without_noise(case14_pf, target):
    # With e = k s_i on measurement i only, r_i = S_ii e_i in the linear
    # model, so its normalised residual is k sqrt(S_ii) and it is the largest.
    meas = full_measurement_set(case14_pf)
    i = row(meas, *target)
    k = 10
    s_ii = sensitivity_diagonal(run(meas))
    meas.loc[i, "value"] += k * meas.loc[i, "std_dev"]
    res = run(meas)
    rn = normalised_residuals(res)
    pos = np.flatnonzero(res.meas_index == i)[0]
    assert np.nanargmax(rn) == pos
    assert rn[pos] == pytest.approx(k * np.sqrt(s_ii[pos]), rel=0.05)


def test_lnr_removes_the_bad_measurement_like_pandapower(case14_pf, rng):
    meas = add_noise(full_measurement_set(case14_pf), rng)
    i = row(meas, "q", "line", 2)
    meas.loc[i, "value"] += 20 * meas.loc[i, "std_dev"]

    removed, res = lnr_identify(copy.deepcopy(CASE14), meas)
    assert removed == [i]
    assert np.nanmax(normalised_residuals(res)) <= 3.0

    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    se = StateEstimation(net, tolerance=1e-8, algorithm="wls")
    assert se.perform_rn_max_test()
    assert set(meas.index) - set(net.measurement.index) == {i}
