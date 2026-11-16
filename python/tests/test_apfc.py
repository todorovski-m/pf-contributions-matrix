"""Edge cases of the APFC port (T13a, acceptance 6), on small hand-made networks, plus checks of the labels,
the phases and write_csv. The input frames follow load_results: PT is the flow delivered at the to-end."""
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from conftest import MATRIX_DIR
from pfc_matrix import Case, apfc, assemble_active, build_active_table, prepare_active, solve_active, write_csv
from pfc_matrix.labels import unique_labels


def make_case(buses, gens, branches, name="t"):
    """buses: {bus: PD}; gens: [(bus, PG)]; branches: [(F, T, PF, PT)] with PT delivered at the to-end."""
    bus = pd.DataFrame({"BUS_I": list(buses), "PD": [float(v) for v in buses.values()],
                        "QD": 0.0, "VM": 1.0, "VA": 0.0})
    gen = pd.DataFrame({"GEN_BUS": [g[0] for g in gens], "PG": [float(g[1]) for g in gens], "QG": 0.0})
    branch = pd.DataFrame({"F_BUS": [b[0] for b in branches], "T_BUS": [b[1] for b in branches],
                           "PF": [float(b[2]) for b in branches], "QF": 0.0,
                           "PT": [float(b[3]) for b in branches], "QT": 0.0})
    return Case.from_frames(bus, gen, branch, name)


BASE_BRANCHES = [(1, 2, 70, 70), (2, 3, 10, 10), (1, 3, 30, 30)]


def test_base_network():
    table, info = apfc(make_case({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES))
    assert list(table.columns) == ["G_1"]
    assert list(table.index) == ["D_2", "D_3", "L_{f.1-2}", "L_{f.2-3}", "L_{f.1-3}",
                                 "L_{t.1-2}", "L_{t.2-3}", "L_{t.1-3}"]
    np.testing.assert_allclose(table["G_1"].to_numpy(), [60, 40, 70, 10, 30, 70, 10, 30], rtol=0, atol=1e-12)
    assert info["n_e"] == 5 and info["n_bus"] == 3 and info["n_branch"] == 3 and info["n_demand"] == 2
    assert info["issparse"] == {"N": True, "A": True, "B": True}
    assert info["removed"] == {"branches": 0, "generators": 0} and info["flipped_branches"] == 0


def test_negative_pd_is_source_column():
    # bus 3 has PD = -20: a source (column D_3); bus 3 sends 30 + 20 = 50 to bus 2
    table, info = apfc(make_case({1: 0, 2: 120, 3: -20}, [(1, 100)], [(1, 2, 70, 70), (1, 3, 30, 30), (3, 2, 50, 50)]))
    assert list(table.columns) == ["G_1", "D_3"]
    assert "D_3" not in table.index
    np.testing.assert_allclose(table.loc["D_2"].to_numpy(), [100, 20], rtol=0, atol=1e-12)
    np.testing.assert_allclose(table.loc["L_{f.3-2}"].to_numpy(), [30, 20], rtol=0, atol=1e-12)
    assert info["converted"] == {"negative_demand_to_source": 1, "negative_generation_to_demand": 0}


def test_negative_pg_is_demand_row():
    # the generator at bus 3 has PG = -10: a demand (row G_3), placed after the D rows
    table, info = apfc(make_case({1: 0, 2: 60, 3: 30}, [(1, 100), (3, -10)], [(1, 2, 60, 60), (1, 3, 40, 40)]))
    assert list(table.columns) == ["G_1"]
    assert list(table.index[:3]) == ["D_2", "D_3", "G_3"]
    assert table.loc["G_3", "G_1"] == pytest.approx(10, abs=1e-12)
    assert info["converted"] == {"negative_demand_to_source": 0, "negative_generation_to_demand": 1}
    assert info["n_demand"] == 3


def test_parallel_branches_get_suffixes():
    table, _ = apfc(make_case({1: 0, 2: 90}, [(1, 90)], [(1, 2, 30, 30), (1, 2, 30, 30), (1, 2, 30, 30)]))
    assert list(table.index) == ["D_2", "L_{f.1-2}", "L_{f.1-2}_2", "L_{f.1-2}_3",
                                 "L_{t.1-2}", "L_{t.1-2}_2", "L_{t.1-2}_3"]
    np.testing.assert_allclose(table["G_1"].to_numpy(), [90, 30, 30, 30, 30, 30, 30], rtol=0, atol=1e-12)


def test_branch_with_both_flows_below_tol_is_deleted():
    branches = BASE_BRANCHES + [(2, 3, 5e-7, -5e-7)]
    table, info = apfc(make_case({1: 0, 2: 60, 3: 40}, [(1, 100)], branches))
    assert info["removed"]["branches"] == 1 and info["n_branch"] == 3
    assert list(table.index).count("L_{f.2-3}") == 1 and "L_{f.2-3}_2" not in table.index
    # only one end below tol: kept
    branches = BASE_BRANCHES + [(2, 3, 5e-7, 2e-6)]
    _, info = apfc(make_case({1: 0, 2: 60, 3: 40}, [(1, 100)], branches))
    assert info["removed"]["branches"] == 0 and info["n_branch"] == 4


def test_generator_with_small_pg_is_deleted():
    table, info = apfc(make_case({1: 0, 2: 60, 3: 40}, [(1, 100), (2, 5e-7), (3, -5e-7)], BASE_BRANCHES))
    assert info["removed"]["generators"] == 2
    assert list(table.columns) == ["G_1"] and "G_3" not in table.index
    assert info["converted"]["negative_generation_to_demand"] == 0


def test_branch_with_negative_pf_is_flipped():
    # branch stored as 2 -> 1 with power flowing 1 -> 2: 71 leaves bus 1, 70 arrives at bus 2
    table, info = apfc(make_case({1: 0, 2: 70}, [(1, 71)], [(2, 1, -70, -71)]))
    assert info["flipped_branches"] == 1
    assert list(table.index) == ["D_2", "L_{f.1-2}", "L_{t.1-2}"]
    np.testing.assert_allclose(table["G_1"].to_numpy(), [70, 71, 70], rtol=0, atol=1e-12)


def test_unique_labels_skips_labels_in_use():
    assert unique_labels(["G_2", "G_2", "G_2_2", "G_2"]) == ["G_2", "G_2_3", "G_2_2", "G_2_4"]
    assert unique_labels(["a", "b"]) == ["a", "b"]


def test_bus_numbering_must_be_consecutive():
    bus = pd.DataFrame({"BUS_I": [1, 3], "PD": [0.0, 1.0], "QD": 0.0, "VM": 1.0, "VA": 0.0})
    gen = pd.DataFrame({"GEN_BUS": [1], "PG": [1.0], "QG": 0.0})
    br = pd.DataFrame({"F_BUS": [1], "T_BUS": [3], "PF": [1.0], "QF": 0.0, "PT": [1.0], "QT": 0.0})
    with pytest.raises(ValueError):
        Case.from_frames(bus, gen, br)


def test_phases_compose_to_apfc_and_are_sparse():
    case = make_case({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES)
    prep = prepare_active(case, tol=1e-6)
    system = assemble_active(prep)
    assert all(sp.issparse(m) for m in (system.M, system.N, system.A, system.B))
    X = solve_active(system)
    pd.testing.assert_frame_equal(build_active_table(prep, X, system.a), apfc(case)[0])


def test_write_csv_layout(tmp_path):
    table, _ = apfc(make_case({1: 0, 2: 60, 3: 40}, [(1, 100)], BASE_BRANCHES))
    path = tmp_path / "x_apfc.csv"
    write_csv(table, path)
    lines = path.read_text().split("\n")
    assert lines[0] == "Line,G_1"
    assert lines[1].startswith("D_2,") and lines[-1] == ""
    back = pd.read_csv(path, float_precision="round_trip")
    np.testing.assert_array_equal(back.iloc[:, 1:].to_numpy(), table.to_numpy())


def _load_results():
    path = os.path.join(MATRIX_DIR, "external", "iterative-method", "power_flow.py")
    if not os.path.exists(path):
        pytest.skip("submodule external/iterative-method not checked out")
    spec = importlib.util.spec_from_file_location("power_flow", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.load_results


@pytest.mark.parametrize("name", ["case1", "case9", "case14", "case24_ieee_rts", "case30", "case33", "case118", "tap_loop"])
def test_benchmark_inputs_row_sums_equal_flows(name):
    bus, gen, branch, _ = _load_results()(os.path.join(MATRIX_DIR, "benchmark", "inputs", f"{name}.mat"))
    case = Case.from_frames(bus, gen, branch, name)
    table, info = apfc(case)
    prep = prepare_active(case)
    flows = np.r_[prep.Pf[prep.nl:], prep.Pf[:prep.nl], prep.Pt[:prep.nl]]
    np.testing.assert_allclose(table.to_numpy().sum(axis=1), flows, rtol=0, atol=1e-10)
    assert info["n_e"] == prep.ne
