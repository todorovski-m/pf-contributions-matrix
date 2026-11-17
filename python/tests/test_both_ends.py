"""APFC branches whose two terminals both draw active power (brief T14c): the tests of tests/unit/test_both_ends.m
for the port. The pattern (+,-) as two terminal demands on crafted cases (test_terminal.crafted: the exported case9
with an extra branch 4-5) and on the synthetic case syn9_both_ends, the error pfc:apfc:negativeLoss for (-,+), the
counter converted.case2_branches_to_demands, the order of the terminal-demand rows, invariance to +/-1e-12 at the
terminal flows, the attenuation check at the (0,0) boundary, the equivalence with two loads, and the port against
the Octave output of syn9_both_ends."""
import os

import numpy as np
import pandas as pd
import pytest

from pfc_matrix import Case, apfc, rpfc
from pfc_matrix.inputs import BranchParams
from pfc_matrix import checks
from test_rpfc import compare_with_csv
from test_terminal import TOL, crafted, _case9, _syn, raises, row_sum, map_diff

Q = (1.0, 0.9)


def test_plus_minus_two_terminal_demands():
    _, d0 = apfc(_case9())
    t, d = apfc(crafted(1, -0.1, *Q))
    assert {"C2_4", "C2_5"} <= set(t.index) and "L_{f.4-5}_2" not in t.index
    assert abs(row_sum(t, "C2_4") - 1) < TOL and abs(row_sum(t, "C2_5") - 0.1) < TOL
    assert d["n_branch"] == d0["n_branch"]
    assert d["converted"]["case2_branches_to_demands"] == 1 and d["converted"]["zero_terminal_demands"] == 0
    assert d0["converted"]["case2_branches_to_demands"] == 0


def test_minus_plus_is_negative_loss():
    m = raises("pfc:apfc:negativeLoss", apfc, crafted(-0.1, 1, *Q))
    assert "Active power is delivered at both ends of a branch (negative loss) for branch(es): 10 (4-5)." in m


def test_negative_loss_three_sentences():
    m = raises("pfc:apfc:negativeLoss", checks.check_negative_loss, ["1 (1-2)"], ["2 (3-4)"], ["3 (5-6)"])
    assert m == ("Active power is delivered at the to-end of a branch with zero from-end flow (negative loss) for branch(es): 1 (1-2). "
                 "Active power is delivered at the from-end of a branch with zero to-end flow (negative loss) for branch(es): 2 (3-4). "
                 "Active power is delivered at both ends of a branch (negative loss) for branch(es): 3 (5-6).")
    raises("pfc:apfc:negativeLoss", checks.check_negative_loss, [], ["2 (3-4)"])


def test_order_of_the_terminal_demands():
    c = crafted(0, -0.2, *Q)                                   # branch 10: (0,-), demand at bus 5
    br, par = c.branch.copy(), c.params.frame.copy()
    row = br.iloc[[9]].copy(); row[["PF", "QF", "PT", "QT"]] = [0.3, 1.0, -0.4, 0.9]    # branch 11: (+,-) 4-5
    br = pd.concat([br, row], ignore_index=True); par = pd.concat([par, par.iloc[[9]]], ignore_index=True)
    c2 = Case.from_frames(c.bus, c.gen, br, "order", params=BranchParams(frame=par, base_mva=c.params.base_mva))
    t, d = apfc(c2)
    assert [r for r in t.index if r.startswith("C2_")] == ["C2_4", "C2_5", "C2_5_2"]
    assert abs(row_sum(t, "C2_4") - 0.3) < TOL and abs(row_sum(t, "C2_5") - 0.2) < TOL and abs(row_sum(t, "C2_5_2") - 0.4) < TOL
    assert d["converted"]["case2_branches_to_demands"] == 1 and d["converted"]["zero_terminal_demands"] == 1


def test_attenuation_check_at_the_boundary():
    m = raises("pfc:apfc:nonPositiveAttenuation", apfc, crafted(1e-6, -1e-6, *Q))
    assert "terminal flow at the tolerance; check the input" in m


def _k(case):
    br = case.branch
    return int(np.where((br.F_BUS == 2) & (br.T_BUS == 8))[0][0])


def test_syn9_both_ends():
    c = _syn("syn9_both_ends")
    k = _k(c); r = c.branch.iloc[k]
    t, d = apfc(c)
    assert d["converted"]["case2_branches_to_demands"] == 1
    assert abs(row_sum(t, "C2_2") - r.PF) < TOL and abs(row_sum(t, "C2_8") + r.PT) < TOL
    tq, _ = rpfc(c)
    assert {"E_2_8_f", "E_2_8_t"} <= set(tq.columns)


def test_invariance_to_perturbation_of_the_terminal_flows():
    c = _syn("syn9_both_ends"); k = _k(c)
    t0 = apfc(c)[0]
    for e in (1e-12, -1e-12):
        br = c.branch.copy()
        br.loc[br.index[k], ["PF", "PT"]] += e
        assert map_diff(t0, apfc(Case.from_frames(c.bus, c.gen, br, c.name, params=c.params))[0]) <= TOL


def replace_both_ends(case, k):
    """tests/unit/pfc_replace_both_ends.m on the exported frames: branch k is deleted and replaced by generators with
    negative output, PG = -PF at its from bus and PG = PT (delivered; the branch draws -PT) at its to bus."""
    r = case.branch.iloc[k]
    keep = np.arange(len(case.branch)) != k
    g = case.gen.iloc[[0, 0]].copy()
    g["GEN_BUS"] = [int(r.F_BUS), int(r.T_BUS)]; g["PG"] = [-float(r.PF), float(r.PT)]; g["QG"] = 0.0
    return Case.from_frames(case.bus, pd.concat([case.gen, g], ignore_index=True), case.branch[keep].reset_index(drop=True),
                            case.name + "_replaced",
                            params=BranchParams(frame=case.params.frame[keep].reset_index(drop=True), base_mva=case.params.base_mva))


def equivalence():
    c = _syn("syn9_both_ends"); k = _k(c)
    ta, tb = apfc(c)[0], apfc(replace_both_ends(c, k))[0]
    g2 = [x for x in tb.index if x.startswith("G_2")]; g8 = [x for x in tb.index if x.startswith("G_8")]
    other = map_diff(ta, tb, (), None, ["C2_2", "C2_8"], g2 + g8)
    own = map_diff(ta, tb, [("C2_2", g2[0], 1.0), ("C2_8", g8[0], 1.0)], None,
                   [x for x in ta.index if x not in ("C2_2", "C2_8")], [x for x in tb.index if x not in (g2[0], g8[0])])
    return other, own


def test_equivalence_with_two_loads():
    other, own = equivalence()
    assert other <= TOL and own <= TOL


def test_syn9_both_ends_against_octave():
    octave = os.environ.get("PFC_OCTAVE_OUT")
    if not octave or not os.path.exists(os.path.join(octave, "syn9_both_ends_apfc.csv")):
        pytest.skip("PFC_OCTAVE_OUT (outputs of tests/run_matrix.sh) not set")
    c = _syn("syn9_both_ends")
    compare_with_csv(apfc(c)[0], os.path.join(octave, "syn9_both_ends_apfc.csv"))
    compare_with_csv(rpfc(c)[0], os.path.join(octave, "syn9_both_ends_rpfc.csv"))
