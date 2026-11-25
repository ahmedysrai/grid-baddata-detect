"""How many protected PMUs does case14 need before no SCADA-only attack
is left?

    python experiments/pmu_placement.py [--max-pmus 6]

For each PMU model in gridbad.placement and each number of PMUs, tries
every placement and prints the most states any of them pins, with the
placements that reach it (IEEE bus numbers, 1-based). Then takes the best
three-PMU placement with branch currents and shows what is still open:
the buses whose states can move, and how many SCADA meters an attack in
that direction has to write, and repeats the count with branch currents
plus the zero injection constraint at bus 7. Writes
figures/pmu_placement.png. About 10 seconds with the defaults.
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
from scipy.linalg import null_space

from gridbad.attack import meters_touched
from gridbad.measurements import full_measurement_set, load_measurements
from gridbad.placement import (
    PMU_MODELS,
    best_placements,
    pinned_states,
    pmu_jacobian,
    zero_injection_buses,
    zero_injection_jacobian,
)
from gridbad.wls import estimate

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

FIG = pathlib.Path(__file__).resolve().parent.parent / "figures" / "pmu_placement.png"


def ieee(buses):
    return "{" + ", ".join(str(b + 1) for b in buses) + "}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pmus", type=int, default=6)
    args = ap.parse_args()

    truth = pn.case14()
    pp.runpp(truth, numba=False)
    n_states = 2 * len(truth.bus) - 1
    # (label, PMU model, add zero injection rows)
    cases = [(m, m, False) for m in PMU_MODELS] + [("currents+ZI", "currents", True)]
    counts = {label: [] for label, _, _ in cases}
    best3 = None

    print(f"zero injection buses: {ieee(zero_injection_buses(truth))}")
    print(f"states pinned (of {n_states}) by the best placement")
    for label, model, zi in cases:
        for k in range(1, args.max_pmus + 1):
            count, placements = best_placements(truth, k, model, zero_injection=zi)
            counts[label].append(count)
            shown = " ".join(ieee(p) for p in placements[:3])
            more = f" and {len(placements) - 3} more" if len(placements) > 3 else ""
            print(f"{label:<11} {k} PMUs: {count:>2}  {shown}{more}")
            if label == "currents" and k == 3:
                best3 = placements[0]
            if count == n_states:
                break

    # what the best three are missing, and what an attack there costs
    H_pmu = pmu_jacobian(truth, best3, "currents")
    open_states = np.flatnonzero(~pinned_states(H_pmu, n_states))
    n_buses = len(truth.bus)
    # state i belongs to this bus: angles skip the slack (bus 0), then |V|
    bus_of_state = [*range(1, n_buses), *range(n_buses)]
    open_buses = sorted({bus_of_state[i] for i in open_states})
    c = null_space(H_pmu)[:, 0]
    c *= np.radians(1.0) / np.abs(c[: n_buses - 1]).max()

    net = copy.deepcopy(truth)
    load_measurements(net, full_measurement_set(truth))
    res = estimate(net)
    a = res.H @ c
    print(f"\nbest 3 with currents, {ieee(best3)}: open states {open_states.tolist()}")
    print(f"  buses that can still move: {ieee(open_buses)}")
    print(f"  |H_pmu c| = {np.abs(H_pmu @ c).max():.1e} for a 1 deg null direction")
    print(f"  SCADA meters it has to write: {len(meters_touched(a))} of {res.n_meas}")
    dS = zero_injection_jacobian(truth, zero_injection_buses(truth)) @ c
    dS_mw = dS * truth.sn_mva
    print(f"  change in bus 7 P, Q injection: {dS_mw[0]:.2f} MW, {dS_mw[1]:.2f} Mvar")

    fig, ax = plt.subplots(figsize=(6, 4))
    for label, ys in counts.items():
        ls = "--" if label.endswith("ZI") else "-"
        ax.plot(range(1, len(ys) + 1), ys, marker="o", ls=ls, label=label)
    ax.axhline(n_states, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("number of protected PMUs")
    ax.set_ylabel(f"states pinned (of {n_states})")
    ax.legend(title="PMU measures")
    fig.tight_layout()
    FIG.parent.mkdir(exist_ok=True)
    fig.savefig(FIG, dpi=150)


if __name__ == "__main__":
    main()
