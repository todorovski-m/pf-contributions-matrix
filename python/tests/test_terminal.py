"""Terminal classification of branches with exactly one zero terminal (brief T15b): the tests of
tests/unit/test_terminal.m for the port. The nine sign patterns on crafted cases (the exported case9 with one extra
branch 4-5, R = 0, X = 1e9, B = 0, whose terminal flows are set), the error pfc:apfc:negativeLoss, the counters,
equal tables for the dead end at the from side and at the to side, invariance to a perturbation of the zero
terminal by +/-1e-12, the equivalence with the terminal injection that a dead-end branch represents, and the port
against the Octave outputs of the dead-end cases."""
import os

import numpy as np
import pandas as pd
import pytest

from conftest import MATRIX_DIR
from pfc_matrix import Case, Options, PfcError, apfc, read_branch_params, rpfc
from pfc_matrix.inputs import BranchParams
from pfc_matrix import checks
from test_rpfc import compare_with_csv, real_case

TOL = 1e-10
DEAD = ["syn9_dead_charge_from", "syn9_dead_charge_to", "syn9_dead_reactor_from", "syn9_dead_reactor_to",
        "syn9_dead_shunt_to"]
DEAD_BUS = 10


def _case9():
    return real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "case9.mat"), "case9")


def crafted(pf, pt, qf, qt):
    """Exported case9 with an extra branch 4-5 (index 10) with the terminal flows PF, delivered PT, QF, delivered QT
    (as tests/unit/pfc_crafted_branch.m)."""
    c = _case9()
    br, par = c.branch.copy(), c.params.frame.copy()
    k = int(np.where((br.F_BUS == 4) & (br.T_BUS == 5))[0][0])
    row = br.iloc[[k]].copy()
    row[["PF", "PT", "QF", "QT"]] = [pf, pt, qf, qt]
    prow = par.iloc[[k]].copy()
    prow[["BR_R", "BR_X", "BR_B"]] = [0.0, 1e9, 0.0]
    br = pd.concat([br, row], ignore_index=True)
    par = pd.concat([par, prow], ignore_index=True)
    return Case.from_frames(c.bus, c.gen, br, "crafted", params=BranchParams(frame=par, base_mva=c.params.base_mva))


def _syn(name):
    d = os.environ.get("PFC_SYN_INPUTS")
    if not d or not os.path.exists(os.path.join(d, f"{name}.mat")):
        pytest.skip("PFC_SYN_INPUTS (exports of tests/export_synthetic.sh) not set")
    return real_case(os.path.join(d, f"{name}.mat"), name)


def raises(identifier, fn, *a):
    with pytest.raises(PfcError) as e:
        fn(*a)
    assert e.value.identifier == identifier, e.value.identifier
    return str(e.value)


def row_sum(t, label):
    return float(t.loc[label].sum())


def map_diff(tA, tB, rowpairs=(), colgroups=None, exclA=(), exclB=()):
    """tests/unit/pfc_map_diff.m: largest difference after mapping rows (pairs (labelA, labelB, sign), else the same
    label) and columns (colgroups {labelB: [labelsA]}, else the same label); every row and column that is not left
    out must be mapped."""
    colgroups = colgroups or {}
    cols = {}
    for cb in tB.columns:
        la = colgroups.get(cb, [cb])
        assert all(x in tA.columns for x in la), f"column {cb} of B has no counterpart in A"
        cols[cb] = la
    assert {x for v in cols.values() for x in v} == set(tA.columns), "a column of A is not mapped"
    pairs = {a: (b, s) for a, b, s in rowpairs}
    usedB, mx = set(exclB), 0.0
    assert tA.index.is_unique and tB.index.is_unique
    for ra in tA.index:
        if ra in exclA:
            continue
        rb, s = pairs.get(ra, (ra, 1.0))
        assert rb in tB.index, f"row {ra} of A has no counterpart {rb} in B"
        usedB.add(rb)
        for cb, la in cols.items():
            mx = max(mx, abs(s * float(tA.loc[ra, la].sum()) - float(tB.loc[rb, cb])))
    assert usedB == set(tB.index), f"rows of B not mapped: {set(tB.index) - usedB}"
    return mx


# ---- the nine patterns, APFC (reactive flows of a carrier) ----
def test_apfc_patterns():
    _, d0 = apfc(_case9())
    q = (1.0, 0.9)
    t, d = apfc(crafted(1, 0.9, *q))                                    # (+,+)
    assert "L_{f.4-5}_2" in t.index and d["n_branch"] == d0["n_branch"] + 1
    assert d["converted"]["zero_terminal_demands"] == 0 and d["converted"]["zero_terminal_sources"] == 0
    with pytest.warns(Warning):
        t, d = apfc(crafted(-0.9, -1, *q))                              # (-,-): flipped; cycle 4-5-4 (warning)
    assert "L_{f.5-4}" in t.index and d["flipped_branches"] == d0["flipped_branches"] + 1
    raises("pfc:apfc:nonPositiveAttenuation", apfc, crafted(1, -0.1, *q))     # (+,-)
    raises("pfc:apfc:nonPositiveAttenuation", apfc, crafted(-0.1, 1, *q))     # (-,+)
    t, d = apfc(crafted(1, 0, *q))                                      # (+,0)
    assert "C2_4" in t.index and "L_{f.4-5}_2" not in t.index and abs(row_sum(t, "C2_4") - 1) < TOL
    assert d["n_branch"] == d0["n_branch"] and d["converted"]["zero_terminal_demands"] == 1
    t, d = apfc(crafted(0, -1, *q))                                     # (0,-)
    assert "C2_5" in t.index and abs(row_sum(t, "C2_5") - 1) < TOL and d["converted"]["zero_terminal_demands"] == 1
    m = raises("pfc:apfc:negativeLoss", apfc, crafted(-1, 0, *q))       # (-,0)
    assert "delivered at the from-end" in m and "10 (4-5)" in m
    m = raises("pfc:apfc:negativeLoss", apfc, crafted(0, 1, *q))        # (0,+)
    assert ("Active power is delivered at the to-end of a branch with zero from-end flow (negative loss) "
            "for branch(es): 10 (4-5).") in m
    t, d = apfc(crafted(0, 0, *q))                                      # (0,0): deleted
    assert d["removed"]["branches"] == d0["removed"]["branches"] + 1


def test_negative_loss_check():
    m = raises("pfc:apfc:negativeLoss", checks.check_negative_loss, ["1 (1-2)"], ["2 (3-4)"])
    assert "to-end of a branch" in m and "from-end of a branch" in m and "1 (1-2)" in m and "2 (3-4)" in m
    checks.check_negative_loss([], [])


def test_negative_loss_error_also_with_checks_off():
    raises("pfc:apfc:negativeLoss", apfc, crafted(0, 1, 1, 0.9), Options(checks=False))


# ---- the nine patterns, RPFC (active flows of a carrier) ----
def test_rpfc_patterns():
    _, d0 = rpfc(_case9())
    p = (1.0, 0.9)
    c0 = d0["converted"]
    t, d = rpfc(crafted(*p, 1, 0.9))                                    # (+,+)
    assert "L_{f.4-5}_2" in t.index and d["n_branch"] == d0["n_branch"] + 1
    assert d["converted"]["zero_terminal_demands"] == 0 and d["converted"]["zero_terminal_sources"] == 0
    t, d = rpfc(crafted(*p, -0.9, -1))                                  # (-,-)
    assert "L_{f.5-4}" in t.index and d["flipped_branches"] == d0["flipped_branches"] + 1
    t, d = rpfc(crafted(*p, 1, -0.5))                                   # (+,-): Case 2
    assert {"C2_4", "C2_5"} <= set(t.index)
    assert d["converted"]["case2_branches_to_demands"] == c0["case2_branches_to_demands"] + 1
    assert d["converted"]["zero_terminal_demands"] == 0
    t, d = rpfc(crafted(*p, -1, 0.5))                                   # (-,+): Case 3
    assert {"E_4_5_f", "E_4_5_t"} <= set(t.columns)
    assert d["converted"]["case3_branches_to_sources"] == c0["case3_branches_to_sources"] + 1
    assert d["converted"]["zero_terminal_sources"] == 0
    t, d = rpfc(crafted(*p, 1, 0))                                      # (+,0)
    assert "C2_4" in t.index and "C2_5" not in t.index and abs(row_sum(t, "C2_4") - 1) < TOL
    assert not any(r.startswith("Tf_") for r in t.index)              # Case-2 shunt rule
    assert d["converted"]["zero_terminal_demands"] == 1
    assert d["converted"]["case2_branches_to_demands"] == c0["case2_branches_to_demands"]
    t, d = rpfc(crafted(*p, 0, -1))                                     # (0,-)
    assert "C2_5" in t.index and "C2_4" not in t.index and abs(row_sum(t, "C2_5") - 1) < TOL
    assert d["converted"]["zero_terminal_demands"] == 1
    t, d = rpfc(crafted(*p, -1, 0))                                     # (-,0)
    assert "E_4_5_f" in t.columns and "E_4_5_t" not in t.columns
    assert abs(t.loc["L_{f.4-5}_2", "E_4_5_f"] + 1) < TOL
    assert d["converted"]["zero_terminal_sources"] == 1
    assert d["converted"]["case3_branches_to_sources"] == c0["case3_branches_to_sources"]
    t, d = rpfc(crafted(*p, 0, 1))                                      # (0,+)
    assert "E_4_5_t" in t.columns and "E_4_5_f" not in t.columns
    assert abs(t.loc["L_{t.4-5}_2", "E_4_5_t"] - 1) < TOL and d["converted"]["zero_terminal_sources"] == 1
    t, d = rpfc(crafted(*p, 0, 0))                                      # (0,0): deleted
    assert d["removed"]["branches"] == d0["removed"]["branches"] + 1


# ---- counters ----
COUNTS = {"syn9_dead_charge_from": (1, 0, 0, 1), "syn9_dead_charge_to": (1, 0, 0, 1),
          "syn9_dead_reactor_from": (1, 0, 1, 0), "syn9_dead_reactor_to": (1, 0, 1, 0), "syn9_dead_shunt_to": (1, 0, 0, 0)}


@pytest.mark.parametrize("name", DEAD)
def test_counters_dead_end_cases(name):
    c = _syn(name)
    _, d1 = apfc(c); _, d2 = rpfc(c)
    got = (d1["converted"]["zero_terminal_demands"], d1["converted"]["zero_terminal_sources"],
           d2["converted"]["zero_terminal_demands"], d2["converted"]["zero_terminal_sources"])
    assert got == COUNTS[name]


@pytest.mark.parametrize("name", ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop"])
def test_counters_zero_in_the_benchmark_cases(name):
    c = real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", f"{name}.mat"), name)
    for f in (apfc, rpfc):
        _, d = f(c, Options(checks=False))
        assert d["converted"]["zero_terminal_demands"] == 0 and d["converted"]["zero_terminal_sources"] == 0


# ---- from side and to side ----
@pytest.mark.parametrize("pair", ["charge", "reactor"])
def test_from_and_to_side_equal(pair):
    cf, ct = _syn(f"syn9_dead_{pair}_from"), _syn(f"syn9_dead_{pair}_to")
    assert map_diff(apfc(cf)[0], apfc(ct)[0]) <= TOL
    if pair == "charge":
        mx = map_diff(rpfc(cf)[0], rpfc(ct)[0], [("L_{t.10-5}", "L_{f.5-10}", -1.0)], {"E_5_10_f": ["E_10_5_t"]})
    else:
        mx = map_diff(rpfc(cf)[0], rpfc(ct)[0])
    assert mx <= TOL


# ---- invariance to +/-1e-12 at the zero terminal ----
ZERO_END = {"syn9_dead_charge_from": ["PF", "QF"], "syn9_dead_charge_to": ["PT", "QT"],
            "syn9_dead_reactor_from": ["PF", "QF"], "syn9_dead_reactor_to": ["PT", "QT"], "syn9_dead_shunt_to": ["PT"]}


def _dead_row(case):
    br = case.branch
    k = np.where((br.F_BUS == DEAD_BUS) | (br.T_BUS == DEAD_BUS))[0]
    assert len(k) == 1
    return int(k[0])


@pytest.mark.parametrize("name", DEAD)
def test_invariance_to_perturbation_of_the_zero_terminal(name):
    c = _syn(name)
    a0, q0 = apfc(c)[0], rpfc(c)[0]
    k = _dead_row(c)
    for e in (1e-12, -1e-12):
        br = c.branch.copy()
        br.loc[br.index[k], ZERO_END[name]] = e
        c1 = Case.from_frames(c.bus, c.gen, br, name, params=c.params)
        assert map_diff(a0, apfc(c1)[0]) <= TOL and map_diff(q0, rpfc(c1)[0]) <= TOL


# ---- equivalence with the terminal injection ----
def replace_dead_end(case, shunt):
    """tests/unit/pfc_replace_dead_end.m on the exported frames: the dead-end bus and its branch are deleted, and the
    branch is replaced at the other bus o by a generator with PG = -P_o and QG = -Q_o (P_o, Q_o: the power the
    branch draws from bus o; PT and QT of the export are delivered), or, for the shunt case, QG = 0 and a bus shunt
    power -Q_o added to QS at bus o."""
    k = _dead_row(case)
    r = case.branch.iloc[k]
    if int(r.F_BUS) == DEAD_BUS:
        o, P_o, Q_o = int(r.T_BUS), -float(r.PT), -float(r.QT)
    else:
        o, P_o, Q_o = int(r.F_BUS), float(r.PF), float(r.QF)
    keep = np.arange(len(case.branch)) != k
    br = case.branch[keep].reset_index(drop=True)
    par = case.params.frame[keep].reset_index(drop=True)
    bus = case.bus[case.bus.BUS_I != DEAD_BUS].reset_index(drop=True)
    g = case.gen.iloc[[0]].copy()
    g["GEN_BUS"] = o; g["PG"] = -P_o; g["QG"] = 0.0 if shunt else -Q_o
    gen = pd.concat([case.gen, g], ignore_index=True)
    if shunt:
        bus.loc[bus.BUS_I == o, "QS"] += -Q_o
    rep = Case.from_frames(bus, gen, br, case.name + "_replaced", params=BranchParams(frame=par, base_mva=case.params.base_mva))
    return rep, -Q_o


EQV = {"syn9_dead_charge_from": ("L_{t.10-5}", "E_10_5_t"), "syn9_dead_charge_to": ("L_{f.5-10}", "E_5_10_f"),
       "syn9_dead_reactor_from": None, "syn9_dead_reactor_to": None, "syn9_dead_shunt_to": "shunt"}


def equivalence(name):
    """Largest differences (other rows and columns; the dead-end branch against its replacement) per kind."""
    c = _syn(name)
    shunt = EQV[name] == "shunt"
    rep, inj = replace_dead_end(c, shunt)
    out = {}
    ta, tb = apfc(c)[0], apfc(rep)[0]
    out["apfc"] = (map_diff(ta, tb, (), None, ["C2_5"], ["G_5"]),
                   map_diff(ta, tb, [("C2_5", "G_5", 1.0)], None, [r for r in ta.index if r != "C2_5"],
                            [r for r in tb.index if r != "G_5"]))
    ta, tb = rpfc(c)[0], rpfc(rep)[0]
    if shunt:
        other = map_diff(ta, tb, (), {"S_5": ["S_10", "E_10_5_s"]}, ["L_{f.10-5}", "L_{t.10-5}"], ())
        own = abs(row_sum(ta, "L_{t.10-5}") - inj)
    elif EQV[name] is None:
        other = map_diff(ta, tb, (), None, ["C2_5"], ["G_5"])
        own = map_diff(ta, tb, [("C2_5", "G_5", 1.0)], None, [r for r in ta.index if r != "C2_5"],
                       [r for r in tb.index if r != "G_5"])
    else:
        row, col = EQV[name]
        other = map_diff(ta, tb, (), {"G_5": [col]}, [row], ())
        own = abs(abs(float(ta.loc[row, col])) - inj)
    out["rpfc"] = (other, own)
    return out


@pytest.mark.parametrize("name", DEAD)
def test_equivalence_with_the_terminal_injection(name):
    for kind, (other, own) in equivalence(name).items():
        assert other <= TOL and own <= TOL, (kind, other, own)


# ---- port against the Octave outputs of the dead-end cases ----
@pytest.mark.parametrize("name", DEAD)
def test_dead_end_cases_against_octave(name):
    octave = os.environ.get("PFC_OCTAVE_OUT")
    if not octave or not os.path.exists(os.path.join(octave, f"{name}_apfc.csv")):
        pytest.skip("PFC_OCTAVE_OUT (outputs of tests/run_matrix.sh) not set")
    c = _syn(name)
    compare_with_csv(apfc(c)[0], os.path.join(octave, f"{name}_apfc.csv"))
    compare_with_csv(rpfc(c)[0], os.path.join(octave, f"{name}_rpfc.csv"))
