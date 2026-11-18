from scipy.stats import chi2


def chi2_threshold(n_meas, n_states, confidence=0.95):
    """Threshold for the chi-square bad data test on the WLS objective J(x).

    With Gaussian measurement errors and no bad data, J(x) follows a
    chi-square distribution with m - n degrees of freedom.
    """
    dof = n_meas - n_states
    if dof <= 0:
        raise ValueError(f"no redundancy: {n_meas} measurements, {n_states} states")
    return chi2.ppf(confidence, dof)
