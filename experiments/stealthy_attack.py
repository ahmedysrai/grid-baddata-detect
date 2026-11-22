"""Stealthy FDI on case14: shift the bus 14 angle without tripping either test.

    python experiments/stealthy_attack.py [--trials 50] [--seed 0]

The attacker adds c degrees to the estimated voltage angle of bus 14 and
writes the matching change to every meter it touches, in three ways:

- linear: a = H c with H at the operator's estimate (the DC-style attack)
- ac: a = h(x_hat + c) - h(x_hat), exact for the AC model
- random: same meters and the same weighted size sum((a_i / s_i)^2) as
  the ac attack, but a random direction, so it is plain bad data

The attacker is given the operator's own estimate x_hat, which is the
most it could know. Thresholds are the 95% chi-square point and the Sidak
LNR value, which NOTES.md 2026-11-21 found close to their empirical 5%
points. If the estimator fails to converge on attacked data, that counts
as flagged by both tests, since the operator sees the failure. Prints one
draw in detail, then the fraction of --trials noise draws flagged by each
test, and writes figures/stealthy_attack.png. About 4 minutes with the
defaults.
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

from gridbad.attack import (
    ac_attack,
    apply_attack,
    linear_attack,
    meters_touched,
    state_position,
)
from gridbad.chi2 import chi2_threshold
from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import normalised_residuals, sidak_lnr_threshold
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

CASE14 = pn.case14()
TARGET_BUS = 13  # bus 14 in the IEEE numbering
SHIFTS_DEG = (0, 1, 2, 5, 10, 20, 30)
KINDS = ("linear", "ac", "random")
FIG = pathlib.Path(__file__).resolve().parent.parent / "figures" / "stealthy_attack.png"


def run(meas):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net)


def make_attack(kind, res, c, rng):
    if kind == "linear":
        return linear_attack(res, c)
    a = ac_attack(res, c)
    if kind == "ac":
        return a
    touched = meters_touched(a)
    if len(touched) == 0:
        return a
    b = np.zeros_like(a)
    b[touched] = rng.standard_normal(len(touched)) * res.std[touched]
    return b * np.linalg.norm(a / res.std) / np.linalg.norm(b / res.std)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    exact = full_measurement_set(truth)
    rng = np.random.default_rng(args.seed)

    detail = None
    flagged = {(k, s): [0, 0] for k in KINDS for s in SHIFTS_DEG}
    diverged = {k: 0 for k in KINDS}
    for trial in range(args.trials):
        meas = add_noise(exact, rng)
        clean = run(meas)
        j = state_position(clean, TARGET_BUS, "va")
        j_thr = chi2_threshold(clean.n_meas, clean.n_states)
        rn_thr = sidak_lnr_threshold(clean.n_meas)
        rows = []
        for shift in SHIFTS_DEG:
            c = np.zeros(clean.n_states)
            c[j] = np.radians(shift)
            for kind in KINDS:
                a = make_attack(kind, clean, c, rng)
                if shift == 0:
                    res = clean
                else:
                    try:
                        res = run(apply_attack(meas, clean, a, truth.sn_mva))
                    except RuntimeError:
                        diverged[kind] += 1
                        flagged[kind, shift][0] += 1
                        flagged[kind, shift][1] += 1
                        rows.append((shift, kind, len(meters_touched(a))))
                        continue
                rn = np.nanmax(normalised_residuals(res))
                flagged[kind, shift][0] += res.J > j_thr
                flagged[kind, shift][1] += rn > rn_thr
                moved = np.degrees(res.x[j] - clean.x[j])
                rows.append((shift, kind, len(meters_touched(a)), res.J, rn, moved))
        if trial == 0:
            detail = (clean, j_thr, rn_thr, rows)

    clean, j_thr, rn_thr, rows = detail
    print(
        f"first draw: clean J = {clean.J:.1f}, max r_N = "
        f"{np.nanmax(normalised_residuals(clean)):.2f}"
    )
    print(f"thresholds: J {j_thr:.1f}, max r_N {rn_thr:.2f}")
    print(" shift  attack  meters      J   max r_N  angle moved (deg)")
    for shift, kind, n, *out in rows:
        if shift == 0:
            continue
        if not out:
            print(f"{shift:>5}  {kind:<7} {n:>5}  estimator did not converge")
            continue
        J, rn, moved = out
        print(f"{shift:>5}  {kind:<7} {n:>5} {J:>8.1f} {rn:>8.2f} {moved:>10.2f}")

    print(f"\nfraction flagged over {args.trials} draws (chi2 / LNR)")
    print("shift  " + "  ".join(f"{k:>11}" for k in KINDS))
    for shift in SHIFTS_DEG:
        cells = [
            f"{flagged[k, shift][0] / args.trials:.2f} / "
            f"{flagged[k, shift][1] / args.trials:.2f}"
            for k in KINDS
        ]
        print(f"{shift:>5}  " + "  ".join(cells))
    print("estimator failures: " + ", ".join(f"{k} {n}" for k, n in diverged.items()))

    fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharey=True)
    for ax, t, name in zip(axes, (0, 1), ("chi-square", "LNR")):
        for kind in KINDS:
            rate = [flagged[kind, s][t] / args.trials for s in SHIFTS_DEG]
            ax.plot(SHIFTS_DEG, rate, marker="o", label=kind)
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
