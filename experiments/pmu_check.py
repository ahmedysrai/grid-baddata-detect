"""Does one PMU angle catch the stealthy bus 14 attack?

    python experiments/pmu_check.py [--trials 50] [--seed 0]

Same exact AC attack as stealthy_attack.py (shift the estimated angle of
bus 14 by c degrees), but the measurement set now has one PMU voltage angle
on top of the 82 SCADA meters, and the attacker cannot write to it: its
entry of a is set to zero. Three placements are compared: no PMU, a PMU at
bus 9 (a neighbour of bus 14) and a PMU at bus 14 itself. The first draw is
also run with the attacker owning the bus 14 PMU, as a check that the
attack is then invisible again.

Thresholds and the convergence rule are those of stealthy_attack.py.
Prints the fraction of --trials noise draws flagged by each test and
writes figures/pmu_check.png. About 4 minutes with the defaults.
"""

import argparse
import copy
import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandapower as pp
import pandapower.networks as pn
import pandas as pd

from gridbad.attack import ac_attack, apply_attack, meters_touched, state_position
from gridbad.chi2 import chi2_threshold
from gridbad.measurements import (
    add_noise,
    full_measurement_set,
    load_measurements,
    pmu_angle_measurements,
)
from gridbad.residuals import normalised_residuals, sidak_lnr_threshold
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

CASE14 = pn.case14()
TARGET_BUS = 13  # bus 14 in the IEEE numbering
# label -> PMU buses (0-based)
PLACEMENTS = {"no PMU": [], "PMU at bus 9": [8], "PMU at bus 14": [13]}
SHIFTS_DEG = (0, 0.25, 0.5, 1, 2, 5)
FIG = pathlib.Path(__file__).resolve().parent.parent / "figures" / "pmu_check.png"


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def flags(res):
    """(chi-square fired, LNR fired) at the 5% thresholds."""
    j_thr = chi2_threshold(res.n_meas, res.n_states)
    rn_thr = sidak_lnr_threshold(res.n_meas)
    return res.J > j_thr, np.nanmax(normalised_residuals(res)) > rn_thr


def attack(meas, clean, shift, pmu_rows, sn_mva, owns_pmu=False):
    c = np.zeros(clean.n_states)
    c[state_position(clean, TARGET_BUS, "va")] = np.radians(shift)
    a = ac_attack(clean, c)
    if not owns_pmu:
        a[np.isin(clean.meas_index, pmu_rows)] = 0.0
    try:
        return run(apply_attack(meas, clean, a, sn_mva)), a
    except RuntimeError:
        return None, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    scada = full_measurement_set(truth)
    rng = np.random.default_rng(args.seed)

    flagged = {(p, s): np.zeros(2) for p in PLACEMENTS for s in SHIFTS_DEG}
    failures = 0
    for trial in range(args.trials):
        # one noise draw per trial, shared by all placements; the PMU rows
        # get their own noise after it
        noisy = add_noise(scada, rng)
        for label, buses in PLACEMENTS.items():
            pmu = add_noise(pmu_angle_measurements(truth, buses), rng)
            meas = pd.concat([noisy, pmu], ignore_index=True)
            pmu_rows = meas.index[len(noisy) :]
            clean = run(meas)
            for shift in SHIFTS_DEG:
                if shift == 0:
                    res = clean
                else:
                    res, _ = attack(meas, clean, shift, pmu_rows, truth.sn_mva)
                if res is None:
                    failures += 1
                    flagged[label, shift] += 1
                else:
                    flagged[label, shift] += flags(res)
            if trial == 0 and buses == [TARGET_BUS]:
                print(f"first draw, PMU at bus 14, clean J = {clean.J:.1f}")
                print(" shift  PMU owned  meters      J   max r_N")
                for owns in (False, True):
                    res, a = attack(meas, clean, 2, pmu_rows, truth.sn_mva, owns)
                    rn = np.nanmax(normalised_residuals(res))
                    n = len(meters_touched(a))
                    print(
                        f"{2:>5}  {'yes' if owns else 'no':<9} {n:>6} {res.J:>7.1f} {rn:>8.2f}"
                    )

    print(f"\nfraction flagged over {args.trials} draws (chi2 / LNR)")
    print("shift  " + "  ".join(f"{p:>14}" for p in PLACEMENTS))
    for shift in SHIFTS_DEG:
        cells = [
            "{:.2f} / {:.2f}".format(*flagged[p, shift] / args.trials)
            for p in PLACEMENTS
        ]
        print(f"{shift:>5}  " + "  ".join(f"{c:>14}" for c in cells))
    print(f"estimator failures: {failures}")

    fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharey=True)
    for ax, t, name in zip(axes, (0, 1), ("chi-square", "LNR")):
        for label in PLACEMENTS:
            rate = [flagged[label, s][t] / args.trials for s in SHIFTS_DEG]
            ax.plot(SHIFTS_DEG, rate, marker="o", label=label)
        ax.axhline(0.05, color="grey", lw=0.8, ls=":")
        ax.set_title(name)
        ax.set_xlabel("angle shift at bus 14 (deg)")
    axes[0].set_ylabel("fraction of draws flagged")
    axes[0].legend()
    fig.tight_layout()
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=150)


if __name__ == "__main__":
    main()
