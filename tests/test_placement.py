from itertools import combinations

import numpy as np
import pytest

from gridbad.placement import best_placements, pinned_states, pmu_jacobian

N_STATES = 27  # 14 magnitudes, 13 angles


def state_index(bus, quantity):
    # same layout as attack.state_position with the slack at bus 0
    return bus - 1 if quantity == "va" else 13 + bus


def states_of(buses):
    out = {state_index(b, "vm") for b in buses}
    return out | {state_index(b, "va") for b in buses if b != 0}


def pinned(net, buses, model):
    H = pmu_jacobian(net, buses, model)
    return set(np.flatnonzero(pinned_states(H, N_STATES)))


def neighbours(net):
    adj = {b: {b} for b in net.bus.index}
    pairs = [*net.line[["from_bus", "to_bus"]].values]
    pairs += [*net.trafo[["hv_bus", "lv_bus"]].values]
    for a, b in pairs:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def test_angle_pmu_pins_only_its_own_angle(case14_pf):
    assert pinned(case14_pf, [13], "angle") == {state_index(13, "va")}


def test_phasor_pmu_pins_angle_and_magnitude(case14_pf):
    assert pinned(case14_pf, [13], "phasor") == states_of([13])


def test_current_pmu_on_slack_pins_neighbours(case14_pf):
    # the slack angle is not a state, its magnitude is
    assert pinned(case14_pf, [0], "currents") == states_of([0, 1, 4])


@pytest.mark.parametrize("n_pmus", [2, 3])
def test_current_pmus_pin_exactly_the_covered_buses(case14_pf, n_pmus):
    adj = neighbours(case14_pf)
    for buses in combinations(range(14), n_pmus):
        covered = set().union(*(adj[b] for b in buses))
        assert pinned(case14_pf, buses, "currents") == states_of(covered), buses


def test_four_current_pmus_needed_on_case14(case14_pf):
    best3, _ = best_placements(case14_pf, 3)
    best4, placements = best_placements(case14_pf, 4)
    assert best3 < N_STATES
    assert best4 == N_STATES
    # buses 2, 6, 7, 9 in the IEEE numbering
    assert (1, 5, 6, 8) in placements


def test_unknown_model_raises(case14_pf):
    with pytest.raises(ValueError):
        pmu_jacobian(case14_pf, [0], "voltage")
