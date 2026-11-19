"""Regression tests of rule Q3 version 2.0 (brief T16, S4).

(a) Phase shifters that are reactive carriers: syn9_pst_carrier (branch 4-5, shift 5 degrees, pattern (+,+)) and
    syn9_pst_carrier_tap (branch 9-4, shift 5 degrees, tap 1.05, pattern (-,-), flipped to 4-9). The shifted branch is
    a carrier, so its shunt split enters the attenuation and the self-terms. Checked: the carrier and its series flows
    (Q^s = Q^f + Qsf, Q^e = Q^t - Qst with the line charging of the model), the self-term column, the row sums
    separately for P^f, P^t (APFC) and Q^f, Q^t (RPFC) within 1e-10, and the port against the Octave outputs of
    apfc.m and rpfc.m (python/tests/cases/octave/, made by benchmark/reconstruction/run_octave_cases.sh) within 1e-10.
    The old method fed the line-charging shunts agrees within 1e-10 (benchmark/reconstruction/pst_check.py; it needs
    the old method and is not run here).
(b) The circulating-flow example tap_loop (two parallel transformers, taps 1.05 and 0.95, no line charging) with the
    values of rule Q3 version 2.0: zero shunt flows, attenuations equal to the terminal-flow ratios, the 2x2 carrier
    block of N^Q, its spectral radius, det(I + N^Q), (N^Q)^2, the cycle product, and the RPFC table with every row
    attributed to G_1 (the manuscript's table: 5, 6.449, 6.111, 12.375, 11.449).
"""
import os

import numpy as np
import pandas as pd
import pytest
import scipy.io as sio

from conftest import MATRIX_DIR
from pfc_matrix import apfc, assemble_reactive, prepare_active, prepare_reactive, rpfc
from pfc_matrix.reactive import branch_shunt_flows
from test_rpfc import real_case

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10
CARRIERS = {"syn9_pst_carrier": (4, 5), "syn9_pst_carrier_tap": (4, 9)}    # the shifted branch as a carrier (after a flip)


def _case(name):
    return real_case(os.path.join(HERE, "cases", f"{name}.mat"), name)


@pytest.mark.parametrize("name", sorted(CARRIERS))
def test_phase_shifter_as_reactive_carrier(name):
    case = _case(name)
    prep = prepare_reactive(case, 1e-6)
    system = assemble_reactive(prep)
    f, t = CARRIERS[name]
    k = [i for i in range(prep.nl) if prep.f[i] + 1 == f and prep.t[i] + 1 == t]
    assert len(k) == 1, "the shifted branch is not a carrier"
    k = k[0]
    # the series flows of the carrier use the line charging of the model (after the flip the ends are swapped)
    m = sio.loadmat(os.path.join(HERE, "cases", f"{name}.mat"))
    bm = np.asarray(m["branch"], float); base = float(np.asarray(m["baseMVA"]).ravel()[0])
    j = int(np.where(bm[:, 9] != 0)[0][0])
    qsf, qst = branch_shunt_flows(case)
    vm = case.bus.set_index("BUS_I")["VM"]
    tap = abs(bm[j, 8]) or 1.0
    exp_f = bm[j, 4] / 2 * vm[int(bm[j, 0])] ** 2 / tap ** 2 * base
    exp_t = bm[j, 4] / 2 * vm[int(bm[j, 1])] ** 2 * base
    assert abs(qsf[j] - exp_f) <= TOL and abs(qst[j] - exp_t) <= TOL
    flipped = (int(bm[j, 0]), int(bm[j, 1])) != (f, t)
    sf, st = (exp_t, exp_f) if flipped else (exp_f, exp_t)
    assert abs(prep.Qsf[k] - sf) <= TOL and abs(prep.Qst[k] - st) <= TOL
    assert abs(system.Qs[k] - (prep.Qf[k] + sf)) <= TOL and abs(system.Qe[k] - (prep.Qt[k] - st)) <= TOL
    assert system.aQ[k] > 0
    # its self-term column is in the table and not zero
    table, _ = rpfc(case)
    col = f"E_{f}_{t}_s"
    assert col in table.columns and np.abs(table[col].to_numpy()).max() > 1e-6


@pytest.mark.parametrize("name", sorted(CARRIERS))
def test_phase_shifter_row_sums_and_octave(name):
    case = _case(name)
    for kind, fn in (("apfc", apfc), ("rpfc", rpfc)):
        table, _ = fn(case)
        if kind == "apfc":
            p = prepare_active(case, 1e-6); nd, nl = p.nd, p.nl
            parts = {"demand": p.Pf[nl:], "P^f": p.Pf[:nl], "P^t": p.Pt[:nl]}
        else:
            p = prepare_reactive(case, 1e-6); nd, nl = p.nd, p.nl
            hf, ht = np.asarray(p.case3_hasf, bool), np.asarray(p.case3_hast, bool)
            parts = {"demand": p.Qf[nl:], "Q^f": p.Qf[:nl], "Q^t": p.Qt[:nl], "Case-3": np.r_[-p.case3_QF[hf], p.case3_QT[ht]]}
        rs = table.to_numpy().sum(axis=1)
        start = 0
        for label, flows in parts.items():
            dev = np.abs(rs[start:start + len(flows)] - flows)
            assert (dev <= TOL).all(), (kind, label, float(dev.max()))
            start += len(flows)
        assert start == len(table)
        ref = pd.read_csv(os.path.join(HERE, "cases", "octave", f"{name}_{kind}.csv"), float_precision="round_trip",
                          dtype={"Line": str}).set_index("Line")
        assert list(table.index) == list(ref.index) and list(table.columns) == list(ref.columns)
        np.testing.assert_allclose(table.to_numpy(), ref.to_numpy(float), rtol=0, atol=TOL)


def test_tap_loop_circulating_flow_rule_q3_v2():
    case = real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "tap_loop.mat"), "tap_loop")
    qsf, qst = branch_shunt_flows(case)
    assert np.all(qsf == 0) and np.all(qst == 0)                       # no line charging, so no shunt flows
    prep = prepare_reactive(case, 1e-6)
    system = assemble_reactive(prep)
    nl = prep.nl
    assert nl == 2
    # attenuations: the terminal-flow ratios; the tap-0.95 branch is the carrier 1->2
    k95 = [i for i in range(nl) if prep.f[i] + 1 == 1 and prep.t[i] + 1 == 2][0]
    assert abs(system.aQ[k95] - 0.92516353781898142) <= 1e-12
    assert abs(abs(system.aQ[0] * system.aQ[1]) - 0.87668831860066387) <= 1e-12          # cycle product gamma
    Nc = system.N.toarray()[:nl, :nl]
    assert abs(float(np.max(np.abs(np.linalg.eigvals(Nc)))) - 0.70272805899981738) <= 1e-12
    assert abs(np.linalg.det(np.eye(nl) + Nc) - 0.50617327509434928) <= 1e-12
    np.testing.assert_allclose(Nc @ Nc, 0.4938267249056507 * np.eye(nl), rtol=0, atol=1e-12)
    table, info = rpfc(case)
    assert abs(info["spectral_radius"] - 0.70272805899981738) <= 1e-12
    # the RPFC table: five rows, everything from G_1 (the manuscript's table, to its three decimals)
    assert list(table.index) == ["D_2", "L_{f.2-1}", "L_{f.1-2}", "L_{t.2-1}", "L_{t.1-2}"]
    other = [c for c in table.columns if c != "G_1"]
    assert np.abs(table[other].to_numpy()).max() == 0
    manuscript = {"D_2": 5, "L_{f.2-1}": 6.449, "L_{t.2-1}": 6.111, "L_{f.1-2}": 12.375, "L_{t.1-2}": 11.449}
    for row, value in manuscript.items():
        assert round(float(table.loc[row, "G_1"]), 3) == value
