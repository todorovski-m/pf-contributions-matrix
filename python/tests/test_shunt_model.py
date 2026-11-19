"""Rule Q3 (version 2.0) against the branch model, on every branch of case9, syn9_carrier_tap, syn9_case3_tap, syn9_pst,
syn9_pst_tap, syn9_pst_carrier and syn9_pst_carrier_tap (brief T16, S1 and S4).

Branch model: an ideal transformer with complex ratio tau = tap e^(j shift) at the from-end, followed by the pi circuit
(series admittance ys = 1/(R + jX), line charging B/2 at each end of the circuit), that is the 2x2 stamp
    Ytt = ys + jB/2,  Yff = Ytt/|tau|^2,  Yft = -ys/conj(tau),  Ytf = -ys/tau.
From the exported bus voltages:
  - the stamp reproduces the exported terminal flows (S_f = V_f conj(I_f), S_t = V_t conj(I_t), I = Y V; the exported
    PT and QT are the delivered flows, -S_t);
  - the shunt flows of rule Q3 (branch_shunt_flows) equal the line charging at the internal node and at the to-end:
    Qsf = (B/2)|V_f|^2/|tau|^2 baseMVA, Qst = (B/2)|V_t|^2 baseMVA;
  - the series flows of the method, Q^s = Q^f + Qsf and Q^e = Q^t - Qst, equal the reactive flows of the series
    element: Im(V_i conj(I_s)) and Im(V_t conj(I_s)) baseMVA, with V_i = V_f/tau and I_s = ys (V_i - V_t).
TAP, SHIFT, R, X and B are read from the branch matrix of the exported .mat file. Branches whose terminal flows the
export set to zero (my_runpf.m zeroes flows below 1e-4) are excluded from the flow comparisons; there are none in
these cases, which the test checks.
"""
import os

import numpy as np
import pytest
import scipy.io as sio

from conftest import MATRIX_DIR
from pfc_matrix.reactive import branch_shunt_flows
from test_rpfc import real_case

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10


def _path(name):
    if name == "case9":
        return os.path.join(MATRIX_DIR, "benchmark", "inputs", "case9.mat")
    local = os.path.join(HERE, "cases", f"{name}.mat")
    if os.path.exists(local):
        return local
    d = os.environ.get("PFC_SYN_INPUTS")
    if not d or not os.path.exists(os.path.join(d, f"{name}.mat")):
        pytest.skip("PFC_SYN_INPUTS (exports of tests/export_synthetic.sh) not set")
    return os.path.join(d, f"{name}.mat")


def model(name):
    path = _path(name)
    case = real_case(path, name)
    m = sio.loadmat(path)
    bm = np.asarray(m["branch"], float)
    base = float(np.asarray(m["baseMVA"]).ravel()[0])
    tap = np.abs(bm[:, 8]); tap[tap == 0] = 1.0
    tau = tap * np.exp(1j * np.deg2rad(bm[:, 9]))
    ys = 1 / (bm[:, 2] + 1j * bm[:, 3]); B = bm[:, 4]
    bus = case.bus.set_index("BUS_I")
    V = (bus["VM"] * np.exp(1j * np.deg2rad(bus["VA"]))).to_numpy()
    pos = {b: i for i, b in enumerate(bus.index)}
    Vf = V[[pos[b] for b in case.branch["F_BUS"]]]; Vt = V[[pos[b] for b in case.branch["T_BUS"]]]
    Ytt = ys + 1j * B / 2; Yff = Ytt / np.abs(tau) ** 2; Yft = -ys / np.conj(tau); Ytf = -ys / tau
    If = Yff * Vf + Yft * Vt; It = Ytf * Vf + Ytt * Vt
    Vi = Vf / tau; Is = ys * (Vi - Vt)
    return case, base, dict(Sf=Vf * np.conj(If) * base, St=Vt * np.conj(It) * base,
                            qsf=B / 2 * np.abs(Vf) ** 2 / np.abs(tau) ** 2 * base, qst=B / 2 * np.abs(Vt) ** 2 * base,
                            qs=np.imag(Vi * np.conj(Is)) * base, qe=np.imag(Vt * np.conj(Is)) * base,
                            tapped=(tap != 1) | (bm[:, 9] != 0))


CASES = ["case9", "syn9_carrier_tap", "syn9_case3_tap", "syn9_pst", "syn9_pst_tap", "syn9_pst_carrier", "syn9_pst_carrier_tap"]


@pytest.mark.parametrize("name", CASES)
def test_shunt_rule_q3_equals_the_branch_model(name):
    case, base, mdl = model(name)
    br = case.branch
    pf, qf, pt, qt = (br[c].to_numpy(float) for c in ("PF", "QF", "PT", "QT"))
    assert ((np.abs(qf) >= 1e-4) | (np.abs(qt) >= 1e-4)).all()              # no zeroed branch in these cases
    # the stamp reproduces the exported terminal flows (PT, QT delivered)
    np.testing.assert_allclose(pf, mdl["Sf"].real, rtol=0, atol=1e-9)
    np.testing.assert_allclose(qf, mdl["Sf"].imag, rtol=0, atol=1e-9)
    np.testing.assert_allclose(pt, -mdl["St"].real, rtol=0, atol=1e-9)
    np.testing.assert_allclose(qt, -mdl["St"].imag, rtol=0, atol=1e-9)
    # rule Q3: the line charging of the model
    qsf, qst = branch_shunt_flows(case)
    np.testing.assert_allclose(qsf, mdl["qsf"], rtol=0, atol=TOL)
    np.testing.assert_allclose(qst, mdl["qst"], rtol=0, atol=TOL)
    # the series flows of the method are those of the series element of the model
    np.testing.assert_allclose(qf + qsf, mdl["qs"], rtol=0, atol=1e-9)
    np.testing.assert_allclose(qt - qst, mdl["qe"], rtol=0, atol=1e-9)
    if name != "case9":
        assert mdl["tapped"].sum() == 1


def test_phase_shift_does_not_enter_the_shunt_flows():
    # syn9_pst has tap 1 and a shift of 5 degrees on branch 8-9: its shunt flows are those of an untapped line
    case, base, mdl = model("syn9_pst")
    qsf, qst = branch_shunt_flows(case)
    k = int(np.where(mdl["tapped"])[0][0])
    m = sio.loadmat(_path("syn9_pst"))
    B = float(np.asarray(m["branch"], float)[k, 4])
    vm = case.bus.set_index("BUS_I")["VM"]
    assert abs(qsf[k] - B / 2 * vm[8] ** 2 * base) <= TOL and abs(qst[k] - B / 2 * vm[9] ** 2 * base) <= TOL
