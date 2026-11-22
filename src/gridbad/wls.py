"""Weighted least squares state estimation, using pandapower's solver.

Notation: z = h(x) + e, R = cov(e) diagonal, J(x) = sum((z_i - h_i(x))^2 / s_i^2).
States are the bus voltage magnitudes and angles, less the slack angle, so
a network with N buses has n = 2N - 1 states.
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from pandapower.estimation.algorithm.matrix_base import BaseAlgebra
from pandapower.estimation.state_estimation import StateEstimation


@dataclass
class WlsResult:
    vm_pu: np.ndarray
    va_degree: np.ndarray
    residual: np.ndarray  # z - h(x_hat), per unit, in pandapower's internal order
    std: np.ndarray  # measurement standard deviations, per unit, same order
    H: np.ndarray  # Jacobian dh/dx at x_hat, n_meas x n_states
    meas_index: np.ndarray  # row of net.measurement behind each entry of z
    x: np.ndarray  # state at the estimate: non-slack angles (rad), then all |V|
    h: Callable[[np.ndarray], np.ndarray]  # h(x) for any state, same order as z
    J: float
    n_meas: int
    n_states: int


def estimate(net, tolerance=1e-8):
    """Run WLS on the measurements in `net.measurement` from a flat start.

    No zero-injection constraints are added, so every row in the
    measurement table is one row of z and nothing else is.
    """
    se = StateEstimation(net, tolerance=tolerance, algorithm="wls")
    out = se.estimate(zero_injection=None)
    if not out["success"]:
        raise RuntimeError("WLS did not converge")
    eppci = se.solver.eppci
    r = eppci.z - se.solver.hx
    s = eppci.r_cov
    algebra = BaseAlgebra(eppci)
    x = eppci.E.copy()
    H = algebra.create_hx_jacobian(x).toarray()
    n_buses = len(net.bus)
    return WlsResult(
        vm_pu=net.res_bus_est.vm_pu.to_numpy(),
        va_degree=net.res_bus_est.va_degree.to_numpy(),
        residual=r,
        std=s,
        H=H,
        meas_index=np.asarray(se.solver.pp_meas_indices),
        x=x,
        h=algebra.create_hx,
        J=float(np.sum((r / s) ** 2)),
        n_meas=len(r),
        n_states=2 * n_buses - 1,
    )
