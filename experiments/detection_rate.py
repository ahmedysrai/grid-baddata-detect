"""Monte Carlo detection rate of the chi-square and LNR tests against k.

    python experiments/detection_rate.py [--trials 200] [--seed 0]

For each trial a fresh noise draw is made. The same draw is used for k = 0
(no bad data, which gives the false alarm rate) and for an error of k
standard deviations on each of the three targets from gross_error.py, so
the two tests and the three targets are compared on common noise.

A trial counts as detected by chi-square when J(x) is above the 95%
threshold, and by LNR when the largest normalised residual is above 3.
LNR also counts as identified when that largest one is the measurement
carrying the error. Writes figures/detection_rate.png and prints the
table behind it. About 10 minutes for 200 trials.
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

from gridbad.chi2 import chi2_threshold
from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import normalised_residuals
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

CASE14 = pn.case14()
LNR_THRESHOLD = 3.0
KS = (2, 4, 6, 8, 10, 12, 14)
# (meas_type, element_type, element) as in gross_error.py
TARGETS = {
    "P injection, bus 4": ("p", "bus", 3),
    "Q flow, line 2-3": ("q", "line", 2),
    "V magnitude, bus 14": ("v", "bus", 13),
}
FIG = pathlib.Path(__file__).resolve().parent.parent / "figures" / "detection_rate.png"


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def outcome(res, thr, target_label=None):
    """(chi-square flag, LNR flag, LNR picked the target) for one estimate."""
    rn = normalised_residuals(res)
    i = np.nanargmax(rn)
    lnr = rn[i] > LNR_THRESHOLD
    picked = lnr and target_label is not None and res.meas_index[i] == target_label
    return res.J > thr, lnr, picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    base = full_measurement_set(truth)
    labels = {}
    for name, (mt, et, el) in TARGETS.items():
        mask = (base.meas_type == mt) & (base.element_type == et) & (base.element == el)
        labels[name] = base.index[mask][0]
    thr = chi2_threshold(len(base), 2 * len(truth.bus) - 1)

    rng = np.random.default_rng(args.seed)
    clean = np.zeros((args.trials, 3), dtype=bool)
    hits = np.zeros((len(TARGETS), len(KS), args.trials, 3), dtype=bool)
    for t in range(args.trials):
        meas = add_noise(base, rng)
        clean[t] = outcome(run(meas), thr)
        for a, label in enumerate(labels.values()):
            for b, k in enumerate(KS):
                bad = meas.copy()
                bad.loc[label, "value"] += k * bad.loc[label, "std_dev"]
                hits[a, b, t] = outcome(run(bad), thr, label)

    fa = clean.mean(axis=0)
    rate = hits.mean(axis=2)
    print(f"{args.trials} trials, seed {args.seed}, chi-square threshold {thr:.1f}")
    print(f"false alarms at k = 0: chi-square {fa[0]:.3f}, LNR {fa[1]:.3f}")
    print()
    print(f"{'':<20} {'k':>3}  chi2   LNR  LNR right")
    for a, name in enumerate(TARGETS):
        for b, k in enumerate(KS):
            c, n, p = rate[a, b]
            print(f"{name if b == 0 else '':<20} {k:>3}  {c:4.2f}  {n:4.2f}  {p:4.2f}")

    ks = (0,) + KS
    fig, axes = plt.subplots(1, len(TARGETS), figsize=(6, 4), sharey=True)
    for a, (ax, name) in enumerate(zip(axes, TARGETS)):
        ax.plot(ks, [fa[0], *rate[a, :, 0]], "o-", ms=3, label="chi-square")
        ax.plot(ks, [fa[1], *rate[a, :, 1]], "s-", ms=3, label="LNR > 3")
        ax.plot(ks, [0.0, *rate[a, :, 2]], "--", label="LNR, right one")
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("error (sd)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("fraction of trials")
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle(f"IEEE 14-bus, m = {len(base)}, {args.trials} trials", fontsize=9)
    fig.tight_layout()
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=150)
    print(f"wrote {FIG.relative_to(FIG.parent.parent)}")


if __name__ == "__main__":
    main()
