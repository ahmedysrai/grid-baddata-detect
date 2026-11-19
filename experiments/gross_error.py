"""Add one gross error to a measurement on case14 and watch J(x) move.

    python experiments/gross_error.py [--seed 1]

Prints J for the clean (noisy but unbiased) measurement set, then for the
same noise plus an error of k standard deviations on one measurement,
against the 95% chi-square threshold. A second table gives the residual
sensitivity S_ii of each target and its normalised residual r_N, against
pandapower's default LNR threshold of 3.
"""

import argparse
import copy
import warnings

import numpy as np
import pandapower as pp
import pandapower.networks as pn
import pandas as pd

from gridbad.chi2 import chi2_threshold
from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import normalised_residuals, sensitivity_diagonal
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
# pandapower 3.5 triggers these inside its own estimation code
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)


CASE14 = pn.case14()


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    meas = add_noise(full_measurement_set(truth), np.random.default_rng(args.seed))

    clean = run(meas)
    thr = chi2_threshold(clean.n_meas, clean.n_states)
    print(f"m = {clean.n_meas}, n = {clean.n_states}, threshold (95%) = {thr:.1f}")
    print(f"clean: J = {clean.J:.1f}")

    # three different kinds of measurement, picked by hand
    targets = {
        "P injection, bus 4": (meas.meas_type == "p")
        & (meas.element_type == "bus")
        & (meas.element == 3),
        "Q flow, line 2-3": (meas.meas_type == "q")
        & (meas.element_type == "line")
        & (meas.element == 2),
        "V magnitude, bus 14": (meas.meas_type == "v") & (meas.element == 13),
    }
    ks = (5, 10, 20)
    s_ii = sensitivity_diagonal(clean)
    rn_rows = []
    for name, mask in targets.items():
        idx = meas.index[mask][0]
        pos = np.flatnonzero(clean.meas_index == idx)[0]
        line = []
        rn_line = []
        for k in ks:
            bad = meas.copy()
            bad.loc[idx, "value"] += k * bad.loc[idx, "std_dev"]
            res = run(bad)
            flag = "*" if res.J > thr else " "
            line.append(f"{k:>2} sd: J = {res.J:7.1f}{flag}")
            rn = normalised_residuals(res)
            top = "*" if np.nanargmax(rn) == pos and rn[pos] > 3 else " "
            rn_line.append(f"{rn[pos]:5.1f}{top}")
        print(f"{name:<20} " + "   ".join(line))
        rn_rows.append(f"{name:<20} {s_ii[pos]:5.3f}   " + "  ".join(rn_line))
    print("* above threshold")
    print()
    print(f"{'':<20}  S_ii   r_N at k = " + ", ".join(str(k) for k in ks))
    for r in rn_rows:
        print(r)
    print("* largest r_N in the set and above 3")
    print(f"S_ii over all {clean.n_meas}: min {s_ii.min():.3f}, max {s_ii.max():.3f}")


if __name__ == "__main__":
    main()
