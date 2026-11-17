"""Tests of the options and the validity checks (T13c, step 4). Each test names the Matlab test it mirrors
(tests/unit/test_options.m, test_checks.m, test_diag.m); tests marked 'T13c' are additions of the brief.
Crafted inputs are those of the Matlab tests; indices are 0-based in Python, labels and messages 1-based."""
import csv
import importlib.util
import os
import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import MATRIX_DIR
from pfc_matrix import Case, Options, PfcError, PfcWarning, apfc, read_branch_params, rpfc
from pfc_matrix import checks
from test_apfc import BASE_BRANCHES, make_case as make_active
from test_rpfc import make_case as make_reactive


def raises(identifier, fn, *a, **k):
    with pytest.raises(PfcError) as e:
        fn(*a, **k)
    assert e.value.identifier == identifier, e.value.identifier
    return str(e.value)


def no_warning(fn, *a, **k):
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        out = fn(*a, **k)
    assert not [w for w in rec if issubclass(w.category, PfcWarning)], [str(w.message) for w in rec]
    return out


# ---- test_options.m ----

def test_options_defaults():                                     # 'options: defaults'
    o = Options()
    assert (o.tol, o.tol_div, o.write, o.solved, o.name, o.checks) == (1e-6, 1e-9, True, False, "", True)


def test_options_unknown_field_is_an_error():                    # 'options: unknown field is an error'
    raises("pfc:options:unknownField", Options, tolerance=1)
    raises("pfc:options:unknownField", apfc, make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES), {"tolerance": 1})


def test_options_negative_tol_is_an_error():                     # 'options: negative tol is an error'
    raises("pfc:options:value", Options, tol=-1)
    for bad in ({"tol_div": 0}, {"tol": float("inf")}, {"tol": float("nan")}, {"tol": "1e-6"}, {"name": 3}):
        raises("pfc:options:value", Options, **bad)
    raises("pfc:options:type", apfc, make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES), 42)


def test_keyword_arguments_win_over_options():                   # T13c (step 2)
    case = make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES)
    _, info = apfc(case, Options(tol=1e-3, tol_div=1e-8), tol=1e-5)
    assert info["opts"]["tol"] == 1e-5 and info["opts"]["tol_div"] == 1e-8


def test_numeric_casedata_is_an_error():                         # 'input: numeric casedata is an error'
    raises("pfc:input:casedata", apfc, 42)
    raises("pfc:input:casedata", rpfc, 42)


def test_unsolved_input_is_an_error():                           # 'input: unsolved struct with opts.solved is an error'
    case = make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES)
    raises("pfc:input:notSolved", Case.from_frames, case.bus, case.gen, case.branch.drop(columns=["PF", "PT"]))


def test_non_consecutive_bus_numbering_is_an_error():            # T13c (Python only: apfc.m assumes 1..nb)
    bus = pd.DataFrame({"BUS_I": [1, 3], "PD": [0.0, 1.0], "QD": 0.0, "VM": 1.0, "VA": 0.0})
    gen = pd.DataFrame({"GEN_BUS": [1], "PG": [1.0], "QG": 0.0})
    br = pd.DataFrame({"F_BUS": [1], "T_BUS": [3], "PF": [1.0], "QF": 0.0, "PT": [1.0], "QT": 0.0})
    raises("pfc:input:busNumbering", Case.from_frames, bus, gen, br)


def test_write_and_solved_options_have_no_effect(tmp_path, monkeypatch):
    # 'solved struct: ... identical tables' and 'write: no file written ...': the input of the port is always solved,
    # and apfc and rpfc never write files (write_csv does), whatever write, solved and name say
    monkeypatch.chdir(tmp_path)
    case = make_reactive([1.04, 1.02, 1.0], [0, 0, 22.0], [0, 0, 0], [(1, 40.0)], [(1, 2, 0.1, 0.0, 1.0), (2, 3, 0.1, 0.0, 2.0)])
    for fn in (apfc, rpfc):
        t0, _ = fn(case)
        t1, i1 = fn(case, Options(solved=True, write=True, name="case_named"))
        pd.testing.assert_frame_equal(t0, t1)
        assert i1["opts"]["name"] == "case_named"
    assert list(tmp_path.iterdir()) == []


DIAG_FIELDS = ["n_bus", "n_branch", "n_demand", "n_source", "n_e", "nnz_N", "nnz_A", "issparse", "condest_A",
               "attenuation_range", "removed", "flipped_branches", "converted", "opts"]


def test_info_all_fields_present_and_condest_reproducible():   # 'diag: all fields present', 'diag: condest ... same'
    case = make_reactive([1.04, 1.02, 1.0], [0, 0, 22.0], [0, 0, 0], [(1, 40.0)], [(1, 2, 0.1, 0.0, 1.0), (2, 3, 0.1, 0.0, 2.0)])
    _, i1 = rpfc(case)
    _, i2 = rpfc(case)
    assert all(f in i1 for f in DIAG_FIELDS + ["acyclic", "spectral_radius", "warnings", "negative_series_flow", "case"])
    assert i1["condest_A"] == i2["condest_A"] and i1["condest_A"] >= 1
    state = np.random.get_state()
    rpfc(case)
    assert all(np.array_equal(a, b) for a, b in zip(state[1:2], np.random.get_state()[1:2]))   # global state restored


# ---- test_checks.m ----

M_CHAIN = sp.csc_matrix(([1.0, 1.0], ([1, 2], [0, 1])), shape=(3, 3))   # chain 1 -> 2 -> 3
M2 = sp.csc_matrix(([1.0, 1.0], ([0, 1], [1, 0])), shape=(2, 2))        # 2-cycle


def test_acyclic_chain_and_topological_order():                 # 'acyclic: chain is acyclic with a valid topological order'
    ac, order = checks.kahn(M_CHAIN)
    assert ac and list(order) == [0, 1, 2]


def test_spectral_radius_zero_for_acyclic():                     # 'spectral radius: exactly 0 for an acyclic digraph'
    assert checks.spectral_radius(-0.7 * M_CHAIN, True) == 0.0


def test_cycle_acyclic_active_no_warning():                      # 'cycle: acyclic active digraph gives no warning'
    ac, rho, w = no_warning(checks.check_cycle, -0.7 * M_CHAIN, M_CHAIN, "apfc")
    assert ac and rho == 0 and w == []


def test_cycle_active_rho_below_one_warns():                     # 'cycle: cyclic active digraph with rho < 1 gives warning'
    with pytest.warns(PfcWarning) as rec:
        ac, rho, w = checks.check_cycle(-0.5 * M2, M2, "apfc")
    assert not ac and abs(rho - 0.5) < 1e-15 and w == ["pfc:apfc:cyclic"]
    assert rec[0].message.identifier == "pfc:apfc:cyclic"


def test_cycle_active_rho_at_or_above_one_is_an_error():         # 'cycle: ... rho >= 1 is an error', 'rho = 1 exactly'
    raises("pfc:apfc:spectralRadius", checks.check_cycle, -1.2 * M2, M2, "apfc")
    raises("pfc:apfc:spectralRadius", checks.check_cycle, -1.0 * M2, M2, "apfc")


def test_cycle_reactive():                                       # 'cycle: cyclic reactive digraph with rho < 1 / >= 1'
    ac, rho, w = no_warning(checks.check_cycle, -0.3 * M2, M2, "rpfc")
    assert not ac and w == []
    msg = raises("pfc:rpfc:spectralRadius", checks.check_cycle, -1.5 * M2, M2, "rpfc")
    assert "rho(N^Q) < 1" in msg


def test_spectral_radius_cycle_with_long_tail():                 # 'spectral radius: cycle with a 298-node tail ... (0.4)'
    n = 300
    rows = list(range(1, n)) + [0]
    cols = list(range(0, n - 1)) + [1]
    Mt = sp.csc_matrix((np.ones(n), (rows, cols)), shape=(n, n))
    ac, _ = checks.kahn(Mt)
    assert abs(checks.spectral_radius(-0.4 * Mt, ac) - 0.4) < 1e-12


def test_spectral_radius_eigs_path():                            # 'spectral radius: eigs path for a 5000-node cyclic core'
    k = 2500
    i1, i2 = np.arange(0, 2 * k, 2), np.arange(1, 2 * k, 2)
    Me = sp.csc_matrix((np.ones(2 * k), (np.r_[i1, i2], np.r_[i2, i1])), shape=(2 * k, 2 * k))
    assert abs(checks.spectral_radius(-0.3 * Me, False) - 0.3) < 1e-10


LAB = ["1 (1-2)", "2 (2-3)", "3 (3-4)"]


def test_attenuation():                                          # the four 'attenuation: ...' tests
    msg = raises("pfc:apfc:nonPositiveAttenuation", checks.check_attenuation, [0.9, -0.1, 0.95], LAB, 1e-9)
    assert "2 (2-3)" in msg
    m0 = raises("pfc:apfc:nonPositiveAttenuation", checks.check_attenuation, [0.9, 0.0, 0.95], LAB, 1e-9)
    assert m0 == ("Active attenuation Pt/Pf <= 0 (terminal flow at the tolerance; check the input) "
                  "for branch(es): 2 (2-3).")                          # T14c
    with pytest.warns(PfcWarning) as rec:
        w = checks.check_attenuation([0.9, 1 + 1e-6, 1 + 1e-12], LAB, 1e-9)
    m = str(rec[0].message)
    assert w == ["pfc:apfc:attenuationAboveOne"] and "2 (2-3)" in m and "3 (3-4)" not in m
    assert no_warning(checks.check_attenuation, [0.9, 1.0, 1 + 1e-12], LAB, 1e-9) == []


def test_inflow():                                               # the three 'inflow: ...' tests
    msg = raises("pfc:apfc:zeroInflow", checks.check_inflow, [0, 5, 0], [0, 1], 1e-9, "apfc")
    assert "bus(es) 1," in msg
    raises("pfc:rpfc:zeroInflow", checks.check_inflow, [5, 0, 0], [0, 1], 1e-9, "rpfc")
    checks.check_inflow([5, 4, 0], [0, 1], 1e-9, "apfc")        # the sink has no outgoing branch: no error


def test_series_flow():                                          # the three 'series flow: ...' tests
    lab = ["1 (2-1)", "2 (1-2)"]
    w, idx = no_warning(checks.check_series_flow, [6.4, 12.4], [-2.1, 2.2], lab, 1e-6)
    assert list(idx) == [0] and w == []
    with pytest.warns(PfcWarning) as rec:
        w, idx = checks.check_series_flow([6.4, 12.4], [-2.1, 2.2], lab, 1e-6, True)
    m = str(rec[0].message)
    assert w == ["pfc:rpfc:negativeSeriesFlow"] and list(idx) == [0] and "1 (2-1)" in m and "2 (1-2)" not in m
    w, idx = no_warning(checks.check_series_flow, [6.4, 12.4], [2.1, 2.2], lab, 1e-6, True)
    assert w == [] and len(idx) == 0


def test_selfterm():                                             # the three 'selfterm: ...' tests (T12e)
    checks.check_selfterm([2.1, 0, -5e-10], LAB, 1e-9)
    msg = raises("pfc:rpfc:negativeSelfTerm", checks.check_selfterm, [2.1, -1e-6, 3.0], LAB, 1e-9)
    assert "2 (2-3)" in msg and "1 (1-2)" not in msg and "3 (3-4)" not in msg
    checks.check_selfterm([2.1, -1e-9, 3.0], LAB, 1e-9)          # exactly at -tol_div: no error (strict inequality)
    checks.check_inflow([1e-9, 5, 0], [0, 1], 1e-9, "rpfc")      # as the inflow check at tol_div


# ---- the checks inside apfc and rpfc (T13c) ----

def test_apfc_zero_inflow_error():
    case = make_active({1: 0, 2: 60}, [], [(1, 2, 60, 60)])     # bus 1 has an outgoing branch and no source
    msg = raises("pfc:apfc:zeroInflow", apfc, case)
    assert "bus(es) 1," in msg
    raises("pfc:apfc:zeroInflow", apfc, case, Options(checks=False))   # always on


def test_rpfc_negative_self_term_guard():
    # a to-end shunt flow of -5e-7 lies in [-tol, 0): it is not converted to a demand, and the self-term is -5e-7
    case = make_reactive([1.02, 1.0], [0, 20.0], [0, 0], [(1, 20.4)], [(1, 2, 0.1, 0.0, -5e-7)])
    msg = raises("pfc:rpfc:negativeSelfTerm", rpfc, case)
    assert "1 (1-2)" in msg
    raises("pfc:rpfc:negativeSelfTerm", rpfc, case, Options(checks=False))   # always on


def test_checks_false_switches_off_only_the_optional_checks():
    # attenuation above 1 (8e-6 / 2e-6 = 4; both ends above tol, T15b): warning with checks=True, nothing with
    # checks=False; condest_A None
    case = make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES + [(2, 3, 2e-6, 8e-6)])
    with pytest.warns(PfcWarning):
        _, on = apfc(case)
    _, off = no_warning(apfc, case, Options(checks=False))
    assert on["warnings"] == ["pfc:apfc:attenuationAboveOne"] and off["warnings"] == []
    assert isinstance(on["condest_A"], float) and off["condest_A"] is None
    # attenuation <= 0 (PT < 0 at a branch with PF > 0, pattern (+,-)): error with checks=True, none with
    # checks=False. (PT = 0 exactly is the pattern (+,0) since T15b: a terminal demand, not a carrier.)
    case0 = make_active({1: 0, 2: 60, 3: 40}, [(1, 100)], [(1, 2, 70, 70), (2, 3, 10, -1e-3), (1, 3, 30, 30)])
    raises("pfc:apfc:nonPositiveAttenuation", apfc, case0)
    no_warning(apfc, case0, Options(checks=False))
    # the cyclic warning is optional, the spectral-radius error is not
    ac, rho, w = no_warning(checks.check_cycle, -0.5 * M2, M2, "apfc", warn_cyclic=False)
    assert not ac and w == []
    raises("pfc:apfc:spectralRadius", checks.check_cycle, -1.2 * M2, M2, "apfc", warn_cyclic=False)
    # rpfc: condest_A None, the tables unchanged
    rc = make_reactive([1.04, 1.02, 1.0], [0, 0, 22.0], [0, 0, 0], [(1, 40.0)], [(1, 2, 0.1, 0.0, 1.0), (2, 3, 0.1, 0.0, 2.0)])
    t_on, i_on = rpfc(rc)
    t_off, i_off = rpfc(rc, Options(checks=False))
    pd.testing.assert_frame_equal(t_on, t_off)
    assert i_off["condest_A"] is None and isinstance(i_on["condest_A"], float)


# ---- test_diag.m (real inputs) ----

def _load_results():
    path = os.path.join(MATRIX_DIR, "external", "iterative-method", "power_flow.py")
    if not os.path.exists(path):
        pytest.skip("submodule external/iterative-method not checked out")
    spec = importlib.util.spec_from_file_location("power_flow", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_results


def _case(name):
    mat = os.path.join(MATRIX_DIR, "benchmark", "inputs", f"{name}.mat")
    if not os.path.exists(mat):
        d = os.environ.get("PFC_SYN_INPUTS")
        if not d or not os.path.exists(os.path.join(d, f"{name}.mat")):
            pytest.skip("PFC_SYN_INPUTS (exports of tests/export_synthetic.sh) not set")
        mat = os.path.join(d, f"{name}.mat")
    bus, gen, branch, _ = _load_results()(mat)
    return Case.from_frames(bus, gen, branch, name, params=read_branch_params(mat))


def _expected_carriers():
    exp, covered = {}, set()
    with open(os.path.join(MATRIX_DIR, "tests", "data", "negative_series_flow_carriers.csv")) as fh:
        for r in csv.DictReader(fh):
            covered.add(r["case"]); exp.setdefault(r["case"], [])
            if r["from_bus"]:
                exp[r["case"]].append(f"{r['from_bus']}-{r['to_bus']}")
    return exp, covered


DIAG_CASES = ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop",
              "syn9_neg_pd", "syn9_neg_pg", "syn9_neg_qd", "syn9_neg_qg", "syn9_case3_tap"]


def test_expected_carrier_data_cover_all_cases():               # 'negative series flow: expected data cover all ... cases'
    _, covered = _expected_carriers()
    assert set(DIAG_CASES) <= covered


@pytest.mark.parametrize("name", DIAG_CASES)
def test_diag_case(name):
    # 'diag: <case> apfc/rpfc acyclic with spectral radius 0' and 'diag: <case> rpfc negative_series_flow equals the
    # old method's carrier set, no warning'
    case = _case(name)
    for fn in (apfc, rpfc):
        _, info = no_warning(fn, case)
        if name in ("case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118"):
            assert info["acyclic"] and info["spectral_radius"] == 0
        if fn is rpfc:
            exp, _ = _expected_carriers()
            got = [b.split(" ", 1)[1].strip("()") for b in info["negative_series_flow"]["branches"]]
            assert info["negative_series_flow"]["count"] == len(exp[name]) and set(got) == set(exp[name])
            assert "pfc:rpfc:negativeSeriesFlow" not in info["warnings"]


def test_tap_loop_cyclic_reactive_digraph():                     # 'diag: tap_loop rpfc rho = 0.17937561102916488, ...'
    _, info = no_warning(rpfc, _case("tap_loop"))
    assert abs(info["spectral_radius"] - 0.17937561102916488) <= 1e-12
    assert not info["acyclic"] and info["negative_series_flow"]["branches"] == ["1 (2-1)"]
