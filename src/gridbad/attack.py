"""Stealthy false data injection (FDI) against the WLS estimator.

Liu, Ning and Reiter (ACM CCS 2009, doi:10.1145/1653662.1653666) showed
for the DC model z = H x + e that adding a = H c to the measurements moves
the estimate from x_hat to x_hat + c and leaves the residual unchanged:

    r_a = z + H c - H (x_hat + c) = r.

So J(x) and every normalised residual are the same as without the attack,
and no residual test can tell the two apart. The attacker needs H (the
network model and which meters exist) and write access to every meter
whose row of H c is nonzero.

On the AC model used here, H is only the Jacobian at the estimate, so
a = H c is stealthy to first order. The exact version is
a = h(x_hat + c) - h(x_hat), which makes the residual at x_hat + c equal
to r; the new WLS solution is then x_hat + c up to the change in H between
the two points.

Attack vectors are in z's order and per unit, like `WlsResult.residual`.
"""

import numpy as np

# below this an entry of a counts as zero (per unit)
ZERO_TOL = 1e-12


def state_position(res, bus, quantity):
    """Index into `res.x` of the angle ("va") or magnitude ("vm") of `bus`.

    `bus` is the position of the bus in the network (0-based, as in
    pandapower's case14). The slack angle is not a state, so asking for it
    raises.
    """
    n_buses = (res.n_states + 1) // 2
    if quantity == "vm":
        return n_buses - 1 + bus
    if quantity != "va":
        raise ValueError(f"quantity must be 'va' or 'vm', got {quantity!r}")
    if bus == 0:
        raise ValueError("bus 0 is the slack; its angle is not a state")
    return bus - 1


def linear_attack(res, c):
    """a = H c, with H the Jacobian at the estimate."""
    return res.H @ c


def ac_attack(res, c):
    """a = h(x_hat + c) - h(x_hat)."""
    return res.h(res.x + c) - res.h(res.x)


def meters_touched(a):
    """Positions in z of the measurements the attack has to change."""
    return np.flatnonzero(np.abs(a) > ZERO_TOL)


def apply_attack(meas, res, a, sn_mva):
    """Copy of `meas` with the attack `a` (z order, per unit) added.

    pandapower stores bus injections in the load convention and converts
    them to injections on the way into z, so their sign flips here. Power
    goes from per unit to MW or Mvar on the base `sn_mva`; voltage
    magnitudes stay per unit.
    """
    rows = meas.loc[res.meas_index]
    power = rows.meas_type.isin(["p", "q"]).to_numpy()
    injection = power & (rows.element_type == "bus").to_numpy()
    scale = np.where(power, sn_mva, 1.0) * np.where(injection, -1.0, 1.0)
    out = meas.copy()
    out.loc[res.meas_index, "value"] = rows.value.to_numpy() + scale * a
    return out
