import copy

import numpy as np
import pandapower.networks as pn
import pandas as pd
import pytest

from gridbad.measurements import (
    add_noise,
    full_measurement_set,
    load_measurements,
    pmu_angle_measurements,
)
from gridbad.wls import estimate

CASE14 = pn.case14()  # building it takes about 0.5 s, a deep copy 10 ms


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def test_counts_for_case14(case14_pf):
    res = run(full_measurement_set(case14_pf))
    # 14 buses x (V, P, Q) + 20 branches x (P, Q)
    assert res.n_meas == 82
    assert res.n_states == 27


def test_noise_free_measurements_give_power_flow_state(case14_pf):
    res = run(full_measurement_set(case14_pf))
    np.testing.assert_allclose(res.vm_pu, case14_pf.res_bus.vm_pu, atol=1e-8)
    np.testing.assert_allclose(res.va_degree, case14_pf.res_bus.va_degree, atol=1e-6)
    assert res.J < 1e-10


def test_pmu_angle_adds_one_row_and_fits_exactly(case14_pf):
    meas = pd.concat(
        [full_measurement_set(case14_pf), pmu_angle_measurements(case14_pf, [13])],
        ignore_index=True,
    )
    res = run(meas)
    assert res.n_meas == 83
    pos = np.flatnonzero(res.meas_index == meas.index[-1]).item()
    assert res.h(res.x)[pos] == pytest.approx(
        np.radians(case14_pf.res_bus.va_degree[13])
    )
    assert res.J < 1e-10


def test_shunt_bus_injection_excludes_shunt(case14_pf):
    # case14 has a capacitor at bus 9 (index 8). With the shunt's Q left in
    # the "measured" injection, the noise-free estimate was off by 0.01 pu.
    meas = full_measurement_set(case14_pf)
    q9 = meas[
        (meas.meas_type == "q") & (meas.element_type == "bus") & (meas.element == 8)
    ].value.item()
    assert abs(q9 - case14_pf.load.q_mvar[case14_pf.load.bus == 8].sum()) < 1e-9


def test_mean_J_is_close_to_degrees_of_freedom(case14_pf, rng):
    truth = full_measurement_set(case14_pf)
    n_runs = 60
    J = [run(add_noise(truth, rng)).J for _ in range(n_runs)]
    # J ~ chi2(55): mean 55, sd sqrt(110); sd of the mean over 60 runs is 1.35
    assert abs(np.mean(J) - 55) < 4.5


def test_gross_error_raises_J(case14_pf, rng):
    meas = add_noise(full_measurement_set(case14_pf), rng)
    clean = run(meas).J
    i = meas.index[(meas.meas_type == "v") & (meas.element == 13)][0]
    meas.loc[i, "value"] += 20 * meas.loc[i, "std_dev"]
    # J grows by roughly k^2 times the residual sensitivity of the measurement.
    # experiments/gross_error.py (seed 1) took J from 44 to 373 at k = 20.
    assert run(meas).J > clean + 200
