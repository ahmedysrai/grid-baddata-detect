"""Normalised residuals and the largest normalised residual (LNR) test.

Linearised around the estimate, the residual depends on the measurement
error through r = S e, with the residual sensitivity matrix

    S = I - H G^-1 H^T R^-1,    G = H^T R^-1 H.

Its covariance is Omega = S R = R - H G^-1 H^T. The normalised residual
r_i / sqrt(Omega_ii) is N(0, 1) when there is no bad data, so a threshold
of 3 has the same meaning for every measurement, which J(x) does not give.

A measurement with Omega_ii = 0 is critical: removing it makes the network
unobservable, its residual is always zero, and no residual test can see an
error on it.
"""

import numpy as np
from scipy.stats import norm

from gridbad.measurements import load_measurements
from gridbad.wls import estimate

# S_ii below this counts as zero. With m = n on case14, where every
# measurement is critical, rounding left all S_ii below 2e-14.
CRITICAL_TOL = 1e-10


def sensitivity_diagonal(res):
    """S_ii = Omega_ii / s_i^2 for every measurement, in z's order.

    S is idempotent with trace m - n, so the diagonal sums to the
    degrees of freedom of the chi-square test.
    """
    w = 1.0 / res.std**2
    G = res.H.T @ (w[:, None] * res.H)
    # diag(H G^-1 H^T) without forming the m x m product
    hgh = np.einsum("ij,ji->i", res.H, np.linalg.solve(G, res.H.T))
    return 1.0 - w * hgh


def normalised_residuals(res):
    """|r_i| / sqrt(Omega_ii), NaN for critical measurements."""
    s_ii = sensitivity_diagonal(res)
    rn = np.full(res.n_meas, np.nan)
    ok = s_ii > CRITICAL_TOL
    rn[ok] = np.abs(res.residual[ok]) / (res.std[ok] * np.sqrt(s_ii[ok]))
    return rn


def sidak_lnr_threshold(n_meas, alpha=0.05):
    """LNR threshold giving overall false alarm rate `alpha` over `n_meas`.

    If the m normalised residuals were independent N(0, 1), the largest of
    them stays below t with probability (2 Phi(t) - 1)^m. Setting that to
    1 - alpha gives a per-measurement two-sided level of
    1 - (1 - alpha)^(1/m). They are correlated through S, so the true
    overall rate at this t is a little below `alpha`.
    """
    per_meas = 1.0 - (1.0 - alpha) ** (1.0 / n_meas)
    return norm.isf(per_meas / 2)


def lnr_identify(net, meas, threshold=3.0, max_removals=10):
    """Repeat estimate, find the largest normalised residual, drop it.

    Stops when the largest one is at or below `threshold`. Returns the
    index labels of the removed rows of `meas`, in the order removed, and
    the last estimate. `net` gets the remaining measurements loaded.
    Critical measurements are never removed, because they are never
    flagged.
    """
    removed = []
    current = meas
    for _ in range(max_removals + 1):
        load_measurements(net, current)
        res = estimate(net)
        rn = normalised_residuals(res)
        if np.isnan(rn).all():
            return removed, res
        i = np.nanargmax(rn)
        if rn[i] <= threshold:
            return removed, res
        if len(removed) == max_removals:
            break
        label = res.meas_index[i]
        removed.append(label)
        current = current.drop(label)
    raise RuntimeError(f"still above {threshold} after {max_removals} removals")
