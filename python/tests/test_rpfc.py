"""Edge cases of the RPFC port (T13b, acceptance 4.6) on small crafted inputs and on real inputs.

Crafted inputs: all voltage angles are 0 and R = 0, so the series flow of a branch from bus f to bus t with
reactance X and tap 1 is Qs = Vf (Vf - Vt)/X * baseMVA at the from-end and Qe = Vt (Vf - Vt)/X * baseMVA at the
to-end (the formula of rpfc.m). The terminal flows are set to QF = Qs - Qsf and QT = Qe + Qst for the wanted shunt
flows Qsf, Qst. Real inputs: benchmark/inputs (load_results imported by path) compared with reference_results/.
Optional (skipped if not set): PFC_SYN_INPUTS = folder with the exports of tests/export_synthetic.sh, and
PFC_OCTAVE_OUT = folder with the Octave outputs of tests/run_matrix.sh.
"""
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest

import pfc_matrix.reactive as reactive
from conftest import MATRIX_DIR
from pfc_matrix import (BranchParams, Case, assemble_reactive, build_reactive_table, prepare_reactive,
                        read_branch_params, rpfc, solve_reactive)

BASE = 100.0


def make_case(vm, qd, qs, gens, branches, name="t"):
    """vm, qd, qs: per bus (bus i+1); gens: [(bus, QG)]; branches: [(f, t, X, Qsf, Qst)] with R = 0, tap 1."""
    nb = len(vm)
    bus = pd.DataFrame({"BUS_I": np.arange(1, nb + 1), "PD": 0.0, "QD": np.asarray(qd, float), "VM": np.asarray(vm, float),
                        "VA": 0.0, "QS": np.asarray(qs, float)})
    gen = pd.DataFrame({"GEN_BUS": [g[0] for g in gens], "PG": 0.0, "QG": [float(g[1]) for g in gens]})
    rows = []
    for f, t, x, qsf, qst in branches:
        vf, vt = vm[f - 1], vm[t - 1]
        qs_ser, qe_ser = vf * (vf - vt) / x * BASE, vt * (vf - vt) / x * BASE
        rows.append({"F_BUS": f, "T_BUS": t, "PF": 0.0, "QF": qs_ser - qsf, "PT": 0.0, "QT": qe_ser + qst, "TAP": 1.0})
    branch = pd.DataFrame(rows)
    params = BranchParams(frame=pd.DataFrame({"F_BUS": branch.F_BUS, "T_BUS": branch.T_BUS, "BR_R": 0.0,
                                              "BR_X": [b[2] for b in branches], "BR_B": 0.0}), base_mva=BASE)
    return Case.from_frames(bus, gen, branch, name, params=params)


def expected_row_sums(case, tol=1e-6):
    prep = prepare_reactive(case, tol)
    nl = prep.nl
    return np.r_[prep.Qf[nl:], prep.Qf[:nl], prep.Qt[:nl], -prep.case3_QF, prep.case3_QT]


def check_row_sums(case, table):
    np.testing.assert_allclose(table.to_numpy().sum(axis=1), expected_row_sums(case), rtol=0, atol=1e-10)


def test_case2_terminal_demands_and_no_transformer_shunt_rows():
    # branch 1-2 consumes at both ends (QF > 0, QT < 0) and has an inductive to-end shunt (Qst = -25)
    case = make_case([1.02, 1.0], [0, 0], [0, 0], [(1, 20.4), (2, 5)], [(1, 2, 0.1, 0.0, -25.0)])
    table, info = rpfc(case)
    assert list(table.index) == ["C2_1", "C2_2"]
    assert not any(r.startswith(("Tf_", "Tt_")) for r in table.index)
    assert list(table.columns) == ["G_1", "G_2"]
    np.testing.assert_allclose(table.loc["C2_1"].to_numpy(), [20.4, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(table.loc["C2_2"].to_numpy(), [0, 5], rtol=0, atol=1e-12)
    assert info["converted"]["case2_branches_to_demands"] == 1
    assert info["converted"]["inductive_to_end_shunt_to_demand"] == 0


def test_several_case3_branches_interleaved_and_appended_rows():
    # branches 1-2 and 2-3 inject at both ends (QF < 0, QT > 0): Qsf = 30 > Qs = 20.4, Qst = 10
    case = make_case([1.02, 1.0, 0.98], [9.6, 30.0 + 9.4, 30.0], [0, 0, 0], [],
                     [(1, 2, 0.1, 30.0, 10.0), (2, 3, 0.1, 30.0, 10.0)])
    table, info = rpfc(case)
    assert list(table.columns) == ["E_1_2_f", "E_1_2_t", "E_2_3_f", "E_2_3_t"]
    assert list(table.index[-4:]) == ["L_{f.1-2}", "L_{f.2-3}", "L_{t.1-2}", "L_{t.2-3}"]
    br = case.branch
    np.testing.assert_allclose(table.loc["L_{f.1-2}"].to_numpy(), [br.QF[0], 0, 0, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(table.loc["L_{f.2-3}"].to_numpy(), [0, 0, br.QF[1], 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(table.loc["L_{t.1-2}"].to_numpy(), [0, br.QT[0], 0, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(table.loc["L_{t.2-3}"].to_numpy(), [0, 0, 0, br.QT[1]], rtol=0, atol=1e-12)
    assert info["converted"]["case3_branches_to_sources"] == 2
    check_row_sums(case, table)


def test_case3_branch_with_inductive_shunt_keeps_its_row():
    # Case-3 branch (QF = 20.4 - 30 < 0, QT = 20 - 5 > 0) with an inductive to-end shunt: row Tt_2
    case = make_case([1.02, 1.0], [9.6, 15.0 + 5.0], [0, 0], [], [(1, 2, 0.1, 30.0, -5.0)])
    table, info = rpfc(case)
    assert "Tt_2" in table.index
    assert info["converted"]["inductive_to_end_shunt_to_demand"] == 1
    assert info["converted"]["case3_branches_to_sources"] == 1
    np.testing.assert_allclose(table.loc["Tt_2"].sum(), 5.0, rtol=0, atol=1e-12)
    check_row_sums(case, table)


def test_negative_qd_negative_qg_and_bus_shunts():
    # bus 2: QD = -10 (source D_2), absorbing bus shunt QS = -3 (row Bs_2); bus 3: injecting shunt QS = 5 (S_3),
    # generator with QG = -4 (row G_3)
    case = make_case([1.02, 1.0, 0.98], [0, -10.0, 8.0], [0, -3.0, 5.0], [(1, 50.0), (3, -4.0)],
                     [(1, 2, 0.1, 0.0, 0.0), (2, 3, 0.1, 0.0, 0.0)])
    table, info = rpfc(case)
    assert list(table.columns[:3]) == ["G_1", "S_3", "D_2"]
    assert list(table.columns[3:]) == ["E_1_2_s", "E_2_3_s"]
    assert list(table.index[:3]) == ["D_3", "G_3", "Bs_2"]
    np.testing.assert_allclose(table.loc[["D_3", "G_3", "Bs_2"]].sum(axis=1).to_numpy(), [8.0, 4.0, 3.0], rtol=0, atol=1e-12)
    conv = info["converted"]
    assert (conv["negative_demand_to_source"], conv["negative_generation_to_demand"], conv["inductive_bus_shunt_to_demand"]) == (1, 1, 1)
    check_row_sums(case, table)


def test_inductive_branch_end_shunts_on_a_carrier():
    case = make_case([1.02, 1.0], [0, 18.0], [0, 0], [(1, 23.4)], [(1, 2, 0.1, -3.0, -2.0)])
    table, info = rpfc(case)
    assert list(table.index) == ["D_2", "Tf_1", "Tt_2", "L_{f.1-2}", "L_{t.1-2}"]
    np.testing.assert_allclose(table.loc[["Tf_1", "Tt_2"]].sum(axis=1).to_numpy(), [3.0, 2.0], rtol=0, atol=1e-12)
    assert (info["converted"]["inductive_from_end_shunt_to_demand"], info["converted"]["inductive_to_end_shunt_to_demand"]) == (1, 1)
    check_row_sums(case, table)


def test_self_term_threshold():
    # carrier 1-2 has the self-term 5e-10 <= tol_div: added at its own to-end, not propagated to branch 2-3;
    # carrier 2-3 has the self-term 2 > tol_div: propagated to D_3
    case = make_case([1.04, 1.02, 1.0], [0, 0, 22.0], [0, 0, 0], [(1, 40.0)],
                     [(1, 2, 0.1, 0.0, 5e-10), (2, 3, 0.1, 0.0, 2.0)])
    table, info = rpfc(case)
    assert info["self_toend"][0] == pytest.approx(5e-10, abs=1e-13)
    assert info["self_toend"][1] == pytest.approx(2.0, abs=1e-12)
    assert table.loc["L_{t.1-2}", "E_1_2_s"] == pytest.approx(5e-10, abs=1e-13)
    assert table.loc["L_{f.2-3}", "E_1_2_s"] == 0.0 and table.loc["D_3", "E_1_2_s"] == 0.0
    assert table.loc["L_{t.2-3}", "E_2_3_s"] == pytest.approx(2.0, abs=1e-12)
    assert table.loc["D_3", "E_2_3_s"] > 0


def test_parallel_carriers_get_suffixes():
    case = make_case([1.02, 1.0], [0, 40.0], [0, 0], [(1, 40.8)], [(1, 2, 0.1, 0.0, 0.0), (1, 2, 0.1, 0.0, 0.0)])
    table, _ = rpfc(case)
    assert list(table.columns) == ["G_1", "E_1_2_s", "E_1_2_s_2"]
    assert list(table.index) == ["D_2", "L_{f.1-2}", "L_{f.1-2}_2", "L_{t.1-2}", "L_{t.1-2}_2"]


def test_one_lu_factorisation_for_both_solves(monkeypatch):
    calls = []
    real = reactive.spla.splu
    monkeypatch.setattr(reactive.spla, "splu", lambda *a, **k: calls.append(1) or real(*a, **k))
    case = make_case([1.04, 1.02, 1.0], [0, 0, 22.0], [0, 0, 0], [(1, 40.0)], [(1, 2, 0.1, 0.0, 1.0), (2, 3, 0.1, 0.0, 2.0)])
    prep = prepare_reactive(case)
    system = assemble_reactive(prep)
    assert system.B.shape[1] > 0 and system.B_self.shape[1] == 2
    X, X_self = solve_reactive(system)
    assert len(calls) == 1
    pd.testing.assert_frame_equal(build_reactive_table(prep, system, X, X_self), rpfc(case)[0])


def test_rpfc_needs_branch_parameters():
    case = make_case([1.02, 1.0], [0, 20.0], [0, 0], [(1, 20.4)], [(1, 2, 0.1, 0.0, 0.0)])
    with pytest.raises(ValueError):
        rpfc(Case.from_frames(case.bus, case.gen, case.branch, "t"))


def test_misaligned_parameters_are_rejected():
    case = make_case([1.02, 1.0], [0, 20.0], [0, 0], [(1, 20.4)], [(1, 2, 0.1, 0.0, 0.0)])
    bad = BranchParams(frame=case.params.frame.assign(T_BUS=[1]), base_mva=BASE)
    with pytest.raises(ValueError):
        Case.from_frames(case.bus, case.gen, case.branch, "t", params=bad)


# ---- real inputs ----

def _load_results():
    path = os.path.join(MATRIX_DIR, "external", "iterative-method", "power_flow.py")
    if not os.path.exists(path):
        pytest.skip("submodule external/iterative-method not checked out")
    spec = importlib.util.spec_from_file_location("power_flow", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_results


def real_case(mat, name):
    bus, gen, branch, _ = _load_results()(mat)
    return Case.from_frames(bus, gen, branch, name, params=read_branch_params(mat))


def compare_with_csv(table, path):
    ref = pd.read_csv(path, float_precision="round_trip", dtype={"Line": str})
    assert list(table.index) == list(ref["Line"])
    assert list(table.columns) == list(ref.columns[1:])
    np.testing.assert_allclose(table.to_numpy(), ref.iloc[:, 1:].to_numpy(float), rtol=0, atol=1e-10)


BENCH = ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop"]


@pytest.mark.parametrize("name", BENCH)
def test_benchmark_inputs_row_sums_and_reference_results(name):
    case = real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", f"{name}.mat"), name)
    table, info = rpfc(case)
    check_row_sums(case, table)
    assert info["n_e"] == prepare_reactive(case).ne
    ref = os.path.join(MATRIX_DIR, "reference_results", f"{name}_rpfc.csv")
    if os.path.exists(ref):
        compare_with_csv(table, ref)


def test_repeated_labels_case24_case14_case118():
    t24, _ = rpfc(real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "case24_ieee_rts.mat"), "case24_ieee_rts"))
    # six absorbing generators at bus 22 are demand rows; five injecting generators at bus 15 are columns
    assert [r for r in t24.index if r.startswith("G_22")] == ["G_22", "G_22_2", "G_22_3", "G_22_4", "G_22_5", "G_22_6"]
    assert [c for c in t24.columns if c.startswith("G_15")] == ["G_15", "G_15_2", "G_15_3", "G_15_4", "G_15_5"]
    assert {"C2_3", "C2_24", "Tt_11", "Tt_11_2"} <= set(t24.index)
    t14, _ = rpfc(real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "case14.mat"), "case14"))
    assert {"Tf_4", "Tf_4_2"} <= set(t14.index)
    t118, _ = rpfc(real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "case118.mat"), "case118"))
    assert {"C2_49", "C2_49_2", "C2_49_3", "C2_49_4"} <= set(t118.index)
    for t in (t24, t14, t118):
        assert t.index.is_unique and t.columns.is_unique


def test_tap_loop_negative_attenuation_and_cyclic_digraph():
    case = real_case(os.path.join(MATRIX_DIR, "benchmark", "inputs", "tap_loop.mat"), "tap_loop")
    table, info = rpfc(case)
    assert info["negative_series_flow"] == {"count": 1, "branches": ["1 (2-1)"]}
    system = assemble_reactive(prepare_reactive(case))
    assert (system.aQ[:prepare_reactive(case).nl] < 0).any()
    assert np.isfinite(table.to_numpy()).all()
    check_row_sums(case, table)


def test_extra_case_carrier_with_tap_and_charging():
    # syn9_carrier_tap (T13b): branch 8 (8-9) of case9 with tap 1.05 and B = 0.306 stays a carrier (flipped to 9-8)
    case = real_case(os.path.join(os.path.dirname(__file__), "cases", "syn9_carrier_tap.mat"), "syn9_carrier_tap")
    table, info = rpfc(case)
    assert "E_9_8_s" in table.columns and "L_{t.9-8}" in table.index
    check_row_sums(case, table)
    octave = os.environ.get("PFC_OCTAVE_OUT")
    ref = os.path.join(octave, "syn9_carrier_tap_rpfc.csv") if octave else ""
    if ref and os.path.exists(ref):
        compare_with_csv(table, ref)


def _syn(name):
    d = os.environ.get("PFC_SYN_INPUTS")
    if not d or not os.path.exists(os.path.join(d, f"{name}.mat")):
        pytest.skip("PFC_SYN_INPUTS (exports of tests/export_synthetic.sh) not set")
    return real_case(os.path.join(d, f"{name}.mat"), name)


def test_syn9_case3_tap_case3_branch_with_inductive_shunt():
    case = _syn("syn9_case3_tap")
    table, info = rpfc(case)
    assert info["converted"]["case3_branches_to_sources"] >= 1
    assert info["converted"]["inductive_from_end_shunt_to_demand"] + info["converted"]["inductive_to_end_shunt_to_demand"] >= 1
    assert any(r.startswith(("Tf_", "Tt_")) for r in table.index)
    check_row_sums(case, table)
    octave = os.environ.get("PFC_OCTAVE_OUT")
    if octave and os.path.exists(os.path.join(octave, "syn9_case3_tap_rpfc.csv")):
        compare_with_csv(table, os.path.join(octave, "syn9_case3_tap_rpfc.csv"))


@pytest.mark.parametrize("name", ["syn9_neg_pd", "syn9_neg_pg", "syn9_neg_qd", "syn9_neg_qg"])
def test_synthetic_cases_against_octave(name):
    case = _syn(name)
    table, _ = rpfc(case)
    check_row_sums(case, table)
    octave = os.environ.get("PFC_OCTAVE_OUT")
    if not octave:
        pytest.skip("PFC_OCTAVE_OUT not set")
    compare_with_csv(table, os.path.join(octave, f"{name}_rpfc.csv"))
