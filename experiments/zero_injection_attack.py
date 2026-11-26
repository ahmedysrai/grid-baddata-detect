"""Bus 8 attack with and without the bus 7 zero injection in the estimator.

    python experiments/zero_injection_attack.py [--trials 50] [--seed 0]

pmu_placement.py left one direction open with current PMUs on buses
{2, 6, 9}: the angle and |V| of bus 8, which hangs off bus 7 alone. To
first order those PMUs pin every other state, so the attacker can only
move bus 8. This runs that attack through the AC estimator (exact a =
h(x_hat + c) - h(x_hat), x_hat known to the attacker), once with every
meter as a plain row of z and once with P = Q = 0 at bus 7 as virtual
measurements, which the attacker cannot write. The PMUs themselves are not
in z: pandapower's estimator has no current angle rows.

Thresholds are the 95% chi-square point for m - n = 55 and the Sidak LNR
value for 82 residuals; the virtual rows replace the two bus 7 meters, so
m does not change. Failed estimates count as flagged. Prints one draw in
detail, the fraction of draws flagged, and writes
figures/zero_injection_attack.png. About 4 minutes with the defaults.
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
from gridbad.measurements import add_noise, full_measurement_set, load_measurements
from gridbad.residuals import normalised_residuals, sidak_lnr_threshold
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

CASE14 = pn.case14()
TARGET_BUS = 7  # bus 8 in the IEEE numbering
ZI_BUS = 6  # bus 7
SHIFTS = {
    "va": (0.25, 0.5, 1.0, 2.0, 5.0),  # degrees
    "vm": (0.005, 0.01, 0.02, 0.04),  # per unit
}
SETUPS = {"plain": None, "bus 7 ZI": [ZI_BUS]}
FIG = (
    pathlib.Path(__file__).resolve().parent.parent
    / "figures"
    / "zero_injection_attack.png"
)


def run(meas, zi):
    net = copy.deepcopy(CASE14)
    load_measurements(net, meas)
    return estimate(net, zero_injection=zi)


def shift_vector(res, quantity, size):
    c = np.zeros(res.n_states)
    c[state_position(res, TARGET_BUS, quantity)] = (
        np.radians(size) if quantity == "va" else size
    )
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    exact = full_measurement_set(truth)
    rng = np.random.default_rng(args.seed)

    cases = [(q, s) for q, sizes in SHIFTS.items() for s in sizes]
    flagged = {(name, q, s): [0, 0] for name in SETUPS for q, s in cases}
    clean_flagged = {name: [0, 0] for name in SETUPS}
    failures = 0
    detail = []
    for trial in range(args.trials):
        meas = add_noise(exact, rng)
        for name, zi in SETUPS.items():
            clean = run(meas, zi)
            j_thr = chi2_threshold(clean.n_meas, clean.n_states)
            rn_thr = sidak_lnr_threshold(clean.n_meas)
            clean_flagged[name][0] += clean.J > j_thr
            clean_flagged[name][1] += np.nanmax(normalised_residuals(clean)) > rn_thr
            for q, s in cases:
                c = shift_vector(clean, q, s)
                a = ac_attack(clean, c)
                try:
                    res = run(apply_attack(meas, clean, a, truth.sn_mva), zi)
                except RuntimeError:
                    failures += 1
                    flagged[name, q, s][0] += 1
                    flagged[name, q, s][1] += 1
                    continue
                rn = normalised_residuals(res)
                flagged[name, q, s][0] += res.J > j_thr
                flagged[name, q, s][1] += np.nanmax(rn) > rn_thr
                if trial == 0:
                    k = state_position(clean, TARGET_BUS, q)
                    moved = res.x[k] - clean.x[k]
                    moved = np.degrees(moved) if q == "va" else moved
                    top = res.meas_index[np.nanargmax(rn)]
                    detail.append(
                        (
                            name,
                            q,
                            s,
                            len(meters_touched(a)),
                            res.J,
                            np.nanmax(rn),
                            top,
                            moved,
                        )
                    )
            if trial == 0:
                print(
                    f"{name}: clean J = {clean.J:.1f}, thresholds J {j_thr:.1f}, "
                    f"max r_N {rn_thr:.2f}"
                )

    print("\nfirst draw (largest r_N row: -1 is a virtual ZI row)")
    print("setup      shift          meters       J   max r_N  row  moved")
    for name, q, s, n, J, rn, top, moved in detail:
        unit = "deg" if q == "va" else "pu"
        print(
            f"{name:<9}  {q} {s:<6} {unit:<3} {n:>6} {J:>8.1f} {rn:>8.2f} "
            f"{top:>4}  {moved:.4g}"
        )

    t = args.trials
    print(f"\nfraction flagged over {t} draws (chi2 / LNR)")
    print("shift          " + "  ".join(f"{n:>11}" for n in SETUPS))
    print(
        "none           "
        + "  ".join(
            f"{clean_flagged[n][0] / t:.2f} / {clean_flagged[n][1] / t:.2f}"
            for n in SETUPS
        )
    )
    for q, s in cases:
        unit = "deg" if q == "va" else "pu"
        cells = [
            f"{flagged[n, q, s][0] / t:.2f} / {flagged[n, q, s][1] / t:.2f}"
            for n in SETUPS
        ]
        print(f"{q} {s:<6} {unit:<4} " + "  ".join(cells))
    print(f"estimator failures: {failures}")

    fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharey=True)
    for ax, q, label in zip(
        axes, SHIFTS, ("bus 8 angle shift (deg)", "bus 8 |V| shift (pu)")
    ):
        for name, colour in zip(SETUPS, ("C0", "C1")):
            for test, style in ((0, "-"), (1, "--")):
                rate = [flagged[name, q, s][test] / t for s in SHIFTS[q]]
                tname = "chi-square" if test == 0 else "LNR"
                ax.plot(
                    SHIFTS[q],
                    rate,
                    style,
                    color=colour,
                    marker="o",
                    label=f"{name}, {tname}",
                )
        ax.axhline(0.05, color="grey", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xticks(SHIFTS[q], [f"{s:g}" for s in SHIFTS[q]])
        ax.minorticks_off()
        ax.set_xlabel(label)
    axes[0].set_ylabel("fraction of draws flagged")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=150)


if __name__ == "__main__":
    main()
