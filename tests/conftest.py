import numpy as np
import pandapower as pp
import pandapower.networks as pn
import pytest


@pytest.fixture
def rng():
    return np.random.default_rng(14)


@pytest.fixture(scope="session")
def case14_pf():
    net = pn.case14()
    pp.runpp(net, numba=False)
    return net
