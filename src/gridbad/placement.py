"""Which states a set of protected PMUs pins down, and where to put them.

An attacker who can write every SCADA meter but none of the PMUs can add
any a = H c that is zero on the PMU rows. To first order that leaves the
residual unchanged (see attack.py), so the states the attacker can still
move are those with a nonzero entry in some c in the null space of H_pmu,
the PMU rows of the Jacobian. Bobba et al. (Workshop on Secure Control
Systems, 2010) showed for the DC model that protecting a basic measurement
set is necessary and sufficient to detect these attacks; here that is
H_pmu having full column rank, applied to the AC Jacobian at one point.

Three PMU models are compared:

- "angle": the voltage angle at the bus, as in `pmu_angle_measurements`.
- "phasor": voltage magnitude and angle at the bus.
- "currents": the phasor plus the current phasor at the bus end of every
  branch that touches it. With the branch admittances known, that fixes
  the voltage at each neighbour too.

Current phasors go into H_pmu as real and imaginary parts. The PMU reports
magnitude and angle; the two forms have the same rank wherever I is
nonzero, and the rectangular one has simpler derivatives.

The state order is the estimator's: non-slack angles (rad) by bus, then
every |V| by bus, so a mask from here lines up with `WlsResult.x`.
"""

from itertools import combinations

import numpy as np
import pandapower as pp
from pandapower.pypower.makeYbus import makeYbus
from scipy.linalg import null_space

PMU_MODELS = ("angle", "phasor", "currents")

# a state counts as pinned when its row of the null space basis is below this
PIN_TOL = 1e-8


def _ppc_with_voltages(net):
    if net.res_bus.empty:
        pp.runpp(net, numba=False)
    ppc = net._ppc
    if len(ppc["bus"]) != len(net.bus):
        # three-winding transformers and similar add internal buses
        raise ValueError("network has buses the PMU model does not map")
    return ppc


def bus_voltages(net):
    """Complex bus voltages of the power flow solution, in bus order."""
    va = np.radians(net.res_bus.va_degree.to_numpy())
    return net.res_bus.vm_pu.to_numpy() * np.exp(1j * va)


def voltage_derivatives(V, slack):
    """dV/dx for x = (non-slack angles, all magnitudes): an N x (2N - 1)
    complex matrix with dV_i/dtheta_i = j V_i and dV_i/d|V_i| = V_i / |V_i|.
    """
    n_buses = len(V)
    angle_buses = [b for b in range(n_buses) if b != slack]
    dV = np.zeros((n_buses, 2 * n_buses - 1), dtype=complex)
    for col, b in enumerate(angle_buses):
        dV[b, col] = 1j * V[b]
    for b in range(n_buses):
        dV[b, n_buses - 1 + b] = V[b] / abs(V[b])
    return dV


def pmu_jacobian(net, buses, model="currents"):
    """Rows of the Jacobian, at the power flow state of `net`, for PMUs on
    `buses` (0-based).

    For "angle" each PMU adds d(theta_k)/dx. "phasor" adds d|V_k|/dx as
    well. "currents" adds, for each branch end at bus k, the derivatives of
    Re(I) and Im(I) with I = Y_f V (or Y_t V at the to end), so
    dI/dx = Y_f dV/dx.
    """
    if model not in PMU_MODELS:
        raise ValueError(f"model must be one of {PMU_MODELS}, got {model!r}")
    ppc = _ppc_with_voltages(net)
    slack = int(net.ext_grid.bus.iloc[0])
    V = bus_voltages(net)
    dV = voltage_derivatives(V, slack)
    rows = []
    if model == "currents":
        _, Yf, Yt = makeYbus(ppc["baseMVA"], ppc["bus"], ppc["branch"])
        f_bus = ppc["branch"][:, 0].real.astype(int)
        t_bus = ppc["branch"][:, 1].real.astype(int)
        dIf = Yf @ dV
        dIt = Yt @ dV
    for k in buses:
        rows.append(np.imag(dV[k] / V[k]))
        if model in ("phasor", "currents"):
            rows.append(np.real(np.conj(V[k]) / abs(V[k]) * dV[k]))
        if model == "currents":
            for dI in [*dIf[f_bus == k], *dIt[t_bus == k]]:
                rows.extend([np.real(dI), np.imag(dI)])
    return np.array(rows).reshape(-1, len(dV[0]))


def pinned_states(H_pmu, n_states):
    """Boolean mask over the states: True where every c with H_pmu c = 0
    has c_i = 0, so no attack that leaves the PMUs alone can move state i.
    """
    if len(H_pmu) == 0:
        return np.zeros(n_states, dtype=bool)
    N = null_space(H_pmu)
    if N.shape[1] == 0:
        return np.ones(n_states, dtype=bool)
    return np.linalg.norm(N, axis=1) < PIN_TOL


def best_placements(net, n_pmus, model="currents"):
    """All placements of `n_pmus` PMUs that pin the most states.

    Brute force over every subset of buses, which is 1001 subsets for four
    PMUs on case14. Returns (number of states pinned, list of bus tuples).
    """
    n_states = 2 * len(net.bus) - 1
    best, found = -1, []
    for buses in combinations(net.bus.index, n_pmus):
        count = int(pinned_states(pmu_jacobian(net, buses, model), n_states).sum())
        if count > best:
            best, found = count, [buses]
        elif count == best:
            found.append(buses)
    return best, found
