"""Tests of the Achayuthakan et al. (2010) baseline (brief T14h, step 2).
Hand-checked examples of the paper (lossless): Fig. 1 (six nodes, acyclic; T_d printed as eq. (21)) and Fig. 4 (three
nodes with a circulating flow; T_d printed as eq. (26)); conservation, nonnegativity and the agreement of the two
solution methods on exported cases. Run: python -B -m pytest -p no:cacheprovider benchmark/baseline_achayuthakan/tests"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "timing")))
import achayuthakan as A   # noqa: E402


def lossless(pg, pd_, lines):
    """DataFrames of a lossless network: pg, pd_ per bus (bus i+1), lines [(from, to, flow)] with PF = PT = flow."""
    nb = len(pg)
    bus = pd.DataFrame({"BUS_I": np.arange(1, nb + 1), "PD": np.asarray(pd_, float)})
    gen = pd.DataFrame({"GEN_BUS": [i + 1 for i in range(nb) if pg[i] > 0], "PG": [float(p) for p in pg if p > 0]})
    branch = pd.DataFrame({"F_BUS": [l[0] for l in lines], "T_BUS": [l[1] for l in lines],
                           "PF": [float(l[2]) for l in lines], "PT": [float(l[2]) for l in lines]})
    return bus, gen, branch


# Fig. 1: P = [20 55 25 30 25 20], P_G = [20 35 0 10 0 0], P_D = [0 0 10 15 20 20]
FIG1 = lossless([20, 35, 0, 10, 0, 0], [0, 0, 10, 15, 20, 20],
                [(1, 2, 20), (2, 3, 25), (2, 4, 15), (2, 6, 15), (3, 4, 5), (3, 5, 10), (4, 5, 15), (5, 6, 5)])
TD21 = np.array([[1, 0.3636, 0.3636, 0.2424, 0.2909, 0.3455],
                 [0, 1, 1, 0.6667, 0.8, 0.95],
                 [0, 0, 1, 0.1667, 0.5, 0.125],
                 [0, 0, 0, 1, 0.6, 0.15],
                 [0, 0, 0, 0, 1, 0.25],
                 [0, 0, 0, 0, 0, 1]])
# Fig. 4: generation 200, 100, 250; demand 100, 150, 300; lines I->II 150, II->III 100, III->I 50
FIG4 = lossless([200, 100, 250], [100, 150, 300], [(1, 2, 150), (2, 3, 100), (3, 1, 50)])
TD26 = np.array([[1.0355, 0.6213, 0.1775], [0.0592, 1.0355, 0.2959], [0.2071, 0.1243, 1.0355]])


@pytest.mark.parametrize("data,expected", [(FIG1, TD21), (FIG4, TD26)])
def test_paper_examples_td(data, expected):
    _, info = A.contributions(*data, kind="net", return_info=True)
    Td = np.linalg.inv(info["A"].toarray())
    assert np.abs(Td - expected).max() <= 5e-5          # the paper prints four decimals


def test_paper_fig1_ad_and_contributions():
    _, info = A.contributions(*FIG1, kind="net", return_info=True)
    Ad = info["A"].toarray()
    assert abs(Ad[0, 1] + 20 / 55) < 1e-15 and abs(Ad[1, 3] + 15 / 30) < 1e-15 and abs(Ad[1, 5] + 15 / 20) < 1e-15   # eq. (11)
    t = A.contributions(*FIG1, kind="net")
    # eq. (14): generator I (20) to the demand at VI (20): (20/20) * [T_d]_{1,6} * 20 = 0.3455 * 20
    assert abs(t.loc["D_6", "G_1"] - 20 * TD21[0, 5]) <= 20 * 5e-5
    np.testing.assert_allclose(t.sum(axis=1).to_numpy(), [10, 15, 20, 20], atol=1e-12)


@pytest.mark.parametrize("data", [FIG1, FIG4])
def test_lossless_net_and_gross_agree(data):
    n = A.contributions(*data, kind="net").sort_index(axis=1)
    g = A.contributions(*data, kind="gross").sort_index(axis=1)
    assert list(n.index) == list(g.index) and list(n.columns) == list(g.columns)
    assert np.abs(n.to_numpy() - g.to_numpy()).max() <= 1e-12


def _real(name):
    import cases as C
    d = C.load(name, C.input_dirs())
    return d["bus"], d["gen"], d["branch"]


@pytest.mark.parametrize("name", ["case9", "case118"])
def test_conservation_nonnegativity_and_methods(name):
    bus, gen, branch = _real(name)
    for kind in ("net", "gross"):
        ti, info = A.contributions(bus, gen, branch, kind=kind, method="inverse", return_info=True)
        ts = A.contributions(bus, gen, branch, kind=kind, method="sparse")
        assert np.abs(ti.to_numpy() - ts.to_numpy()).max() <= 1e-10
        assert ts.to_numpy().min() >= -1e-12
        Td = np.linalg.inv(info["A"].toarray()); P = info["P"]; g = info["pg"] / P
        if kind == "net":
            # every bus inflow is explained by the generators: sum_i (P_G,i/P_i) [A_d^-1]_ik P_k = P_k
            np.testing.assert_allclose((g @ Td) * P, P, rtol=0, atol=1e-9)
            np.testing.assert_allclose(ts.sum(axis=1).to_numpy(), info["demand_values"], rtol=0, atol=1e-9)
        else:
            # every generation is allocated to the demands (gross: demands carry their losses)
            pg = gen.groupby("GEN_BUS")["PG"].sum()
            np.testing.assert_allclose(ts.sum(axis=0).to_numpy(), [pg[int(c[2:])] for c in ts.columns], rtol=0, atol=1e-9)
            assert (ts.sum(axis=1).to_numpy() >= info["demand_values"] - 1e-9).all()


def test_condest_1norm():
    _, info = A.contributions(*FIG4, kind="net", return_info=True)
    exact = np.linalg.cond(info["A"].toarray(), 1)
    assert abs(A.condest_1norm(info["A"]) - exact) <= 1e-9 * exact
