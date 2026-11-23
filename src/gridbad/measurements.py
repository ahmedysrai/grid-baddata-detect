"""Measurement sets for state estimation on a pandapower network.

Values are taken from a solved power flow, so the true state is known and
every error in z is one we put there ourselves.
"""

import pandapower as pp
import pandas as pd

# Standard deviations per measurement type. These are assumed round numbers
# for a first run, not taken from any utility's metering spec.
STD_VM_PU = 0.004
STD_POWER = 1.0  # MW or Mvar, i.e. 0.01 pu on the 100 MVA base of case14
# PMU voltage angle, degrees. IEEE C37.118.1-2011: a phase error of 0.57 deg
# (0.01 rad) alone gives the 1% total vector error limit. Taking that limit
# as about 3 sd is an assumption.
STD_VA_DEG = 0.2


def full_measurement_set(net):
    """Voltage magnitude and P, Q injection at every bus, plus P, Q at the
    from end of every line and the hv end of every transformer.

    `net` must hold a converged power flow (`net.res_bus` etc.). Returns a
    DataFrame with one row per measurement and the true value in `value`.

    Bus injections follow pandapower's load convention (consumption is
    positive). Shunt power is left out of them: the estimator treats shunts
    as part of the admittance matrix, so a "measured" injection that
    includes the shunt is wrong by the shunt's output.
    """
    rows = []
    shunt_p = net.res_shunt.p_mw.groupby(net.shunt.bus).sum()
    shunt_q = net.res_shunt.q_mvar.groupby(net.shunt.bus).sum()
    for b in net.bus.index:
        p = net.res_bus.p_mw[b] - shunt_p.get(b, 0.0)
        q = net.res_bus.q_mvar[b] - shunt_q.get(b, 0.0)
        rows.append(("v", "bus", b, None, net.res_bus.vm_pu[b], STD_VM_PU))
        rows.append(("p", "bus", b, None, p, STD_POWER))
        rows.append(("q", "bus", b, None, q, STD_POWER))
    for i in net.line.index:
        rows.append(("p", "line", i, "from", net.res_line.p_from_mw[i], STD_POWER))
        rows.append(("q", "line", i, "from", net.res_line.q_from_mvar[i], STD_POWER))
    for i in net.trafo.index:
        rows.append(("p", "trafo", i, "hv", net.res_trafo.p_hv_mw[i], STD_POWER))
        rows.append(("q", "trafo", i, "hv", net.res_trafo.q_hv_mvar[i], STD_POWER))
    cols = ["meas_type", "element_type", "element", "side", "value", "std_dev"]
    return pd.DataFrame(rows, columns=cols)


def pmu_angle_measurements(net, buses, std_deg=STD_VA_DEG):
    """Voltage angle at each bus in `buses`, in degrees, as a PMU would
    report it.

    The angle is taken relative to the slack bus, which the estimator holds
    at its power flow value. A real PMU measures against a GPS time
    reference, so this assumes the slack angle is known in that reference
    too. Same columns as `full_measurement_set`, so the two can be
    concatenated.
    """
    rows = [("va", "bus", b, None, net.res_bus.va_degree[b], std_deg) for b in buses]
    cols = ["meas_type", "element_type", "element", "side", "value", "std_dev"]
    return pd.DataFrame(rows, columns=cols)


def add_noise(meas, rng):
    """Copy of `meas` with independent Gaussian errors e ~ N(0, std_dev^2)."""
    out = meas.copy()
    out["value"] = out.value + rng.standard_normal(len(out)) * out.std_dev
    return out


def load_measurements(net, meas):
    """Replace the measurement table of `net` with the rows of `meas`.

    Each row keeps its index from `meas`, so estimator results can be
    traced back to it after rows have been dropped.
    """
    net.measurement = net.measurement.iloc[0:0]
    for row in meas.itertuples():
        pp.create_measurement(
            net,
            row.meas_type,
            row.element_type,
            float(row.value),
            float(row.std_dev),
            int(row.element),
            row.side,
            index=row.Index,
        )
