"""Tests of the timing machinery (not of the methods): statistics, the time-out and iteration-cap paths, the repetition
rules, the environment record, the quiet flag, the QSF replacement and the guard on a deliberately wrong input.
Run: python -B -m pytest -p no:cacheprovider benchmark/timing/tests   (from matrix/)"""
import os
import sys
import time

import numpy as np
import pytest

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import timing as T   # noqa: E402
import cases as C    # noqa: E402

CARRIER = os.path.join(C.MATRIX, "python", "tests", "cases")


def test_stats_on_known_samples():
    s = T.stats([1, 2, 3, 4, 100])
    assert (s["n"], s["min"], s["median"], s["mean"]) == (5, 1.0, 3.0, 22.0)
    assert s["std"] == pytest.approx(np.std([1, 2, 3, 4, 100], ddof=1))
    assert s["iqr"] == pytest.approx(4.0 - 2.0)
    assert T.stats([7])["std"] == 0.0 and T.stats([])["n"] == 0


def test_timeout_path_with_a_sleeping_function():
    with pytest.raises(T.RunTimeout):
        T.call_with_timeout(lambda: time.sleep(2), 0.1)
    r = T.measure_cell(lambda: (time.sleep(2), None), budget_s=1, min_rep=2, max_rep=3, warmup=0, timeout_s=0.1)
    assert r["status"] == "did not finish (time-out)" and r["samples"] == []
    assert T.call_with_timeout(lambda: 5, 1.0) == 5            # the timer is cleared after a normal return
    time.sleep(1.1)


def test_iteration_cap_path():
    r = T.measure_cell(lambda: ("x", None), budget_s=0, min_rep=2, max_rep=3, warmup=1, timeout_s=1, finished=lambda res: False)
    assert r["status"] == "did not finish (iteration cap)"


def test_repetition_rules():
    fast = lambda: (None, {"a": 1})
    assert len(T.measure_cell(fast, budget_s=0, min_rep=5, max_rep=200, warmup=2, timeout_s=1)["samples"]) == 5
    r = T.measure_cell(fast, budget_s=1e6, min_rep=5, max_rep=7, warmup=0, timeout_s=1)
    assert len(r["samples"]) == 7 and r["phases"][0] == {"a": 1} and r["status"] == "ok"
    slow = lambda: (time.sleep(0.02), None)
    r = T.measure_cell(slow, budget_s=0.1, min_rep=2, max_rep=200, warmup=0, timeout_s=1)
    assert 2 <= len(r["samples"]) <= 7                         # stops once the budget is used


def test_environment_record_fields():
    rec = T.environment_record(smoke=True)
    for k in ("label", "smoke", "date_utc", "cpu_model", "cores_logical", "memory_bytes", "os", "python", "numpy", "scipy",
              "pandas", "networkx", "threads", "blas", "git_matrix", "git_submodule", "power", "caffeinate_env",
              "caffeinate_process", "ru_maxrss_unit"):
        assert k in rec, k
    assert rec["label"] == T.SMOKE_LABEL and all(v == "1" for v in rec["threads"].values())
    assert T.environment_record(smoke=False)["label"] == "series"


def test_quiet_flag_logic():
    assert T.quiet_flag(1.9, 8) and not T.quiet_flag(2.0, 8) and not T.quiet_flag(3.0, 8)
    assert set(T.load_record()) == {"load1", "load5", "load15"}


def test_qsf_replacement():
    data = C.load("syn9_carrier_tap", [CARRIER])
    branch, rows, old, new = C.correct_qsf(data)
    assert list(rows) == [7] and old[0] != new[0]
    other = [i for i in range(len(branch)) if i != 7]
    assert np.array_equal(branch["QSF"].to_numpy()[other], data["branch"]["QSF"].to_numpy()[other])
    assert branch["QSF"].iloc[7] == new[0] and data["branch"]["QSF"].iloc[7] == old[0]     # the input is not changed
    d9 = C.load("case9", C.input_dirs())
    b9, rows9, _, _ = C.correct_qsf(d9)
    assert len(rows9) == 0 and b9.equals(d9["branch"])


def test_guard_passes_and_fails_on_a_wrong_input():
    data = C.load("syn9_carrier_tap", [CARRIER])
    branch_old, _, _, _ = C.correct_qsf(data)
    good = T.guard(data, "rpfc", 10000, 60, branch_old)
    assert good["ok"] and good["old_status"] == "ok" and good["max_cell"] <= 1e-10
    bad = T.guard(data, "rpfc", 10000, 60, data["branch"])      # the old method fed the uncorrected QSF
    assert not bad["ok"] and bad["max_cell"] > 1e-3


def test_guard_records_an_error_of_the_old_method(monkeypatch):
    data = C.load("case9", C.input_dirs())
    branch_old, _, _, _ = C.correct_qsf(data)

    def raising(*a, **k):
        return 1 / 0
    monkeypatch.setattr(T, "old_with_iterations", raising)
    g = T.guard(data, "apfc", 10000, 60, branch_old)
    assert g["ok"] and g["iterations"] is None
    assert g["old_status"].startswith("did not finish (error: ZeroDivisionError: division by zero (at test_timing.py:")


def test_port_runner_phases():
    data = C.load("case9", C.input_dirs())
    table, ph = T.port_runner(data, "rpfc", True)()
    assert set(ph) == {"input", "prepare", "assemble", "solve", "table", "condest"} and all(v >= 0 for v in ph.values())
    _, ph = T.port_runner(data, "apfc", False)()
    assert "condest" not in ph


def test_memory_limit_and_port_errors_skip_cells(tmp_path):
    # a memory fraction of 1e-12 skips every cell of case9; the non-consecutive numbering of a case is recorded
    T.main(["run", "--out", str(tmp_path / "r"), "--cases", "case9", "--budget", "0", "--min-rep", "1", "--max-rep", "1",
            "--warmup", "0", "--timeout", "30", "--no-memory", "--mem-fraction", "1e-12", "--smoke"])
    import csv
    st = [r["status"] for r in csv.DictReader(open(tmp_path / "r" / "cells.csv"))]
    assert len(st) == 2 and all(x.startswith("not timed (estimated dense tables") for x in st)


def test_dense_bytes_formula():
    class P: ne, ng, nd, nl, ns, n_case3 = 10, 2, 3, 4, 5, 1
    assert T.dense_bytes(P, "apfc") == 8 * (10 * 2 + (3 + 8) * 2)
    assert T.dense_bytes(P, "rpfc") == 8 * (2 * 10 * 9 + (3 + 8 + 2) * 9)


def test_slow_cell_rule():
    # a single run above the threshold: no further warm-up, exactly 3 recorded repetitions, flagged slow; min_rep and
    # the budget do not apply. Below the threshold the usual rules apply and the cell is not slow.
    calls = []

    def slow():
        calls.append(1); time.sleep(0.03)
        return None, None
    r = T.measure_cell(slow, budget_s=1e6, min_rep=5, max_rep=200, warmup=2, timeout_s=5, slow_s=0.01)
    assert r["slow"] and r["status"] == "ok" and len(r["samples"]) == 3 and r["warmup"] == 1 and len(calls) == 4
    fast = lambda: (None, None)
    r = T.measure_cell(fast, budget_s=0, min_rep=5, max_rep=200, warmup=2, timeout_s=1, slow_s=0.01)
    assert not r["slow"] and len(r["samples"]) == 5
    # a cell that becomes slow only in a recorded run also stops at 3
    n = []

    def grow():
        n.append(1); time.sleep(0.03 if len(n) > 1 else 0)
        return None, None
    r = T.measure_cell(grow, budget_s=1e6, min_rep=5, max_rep=200, warmup=0, timeout_s=5, slow_s=0.01)
    assert r["slow"] and len(r["samples"]) == 3


def test_combination_list():
    assert C.parse_combos("case9:apfc,case9:rpfc,case3120sp:rpfc") == [("case9", ("apfc", "rpfc")), ("case3120sp", ("rpfc",))]
    assert all(c in C.SMALL for c, k in C.COMBOS[:len(C.SMALL)]) and all(k == ("apfc", "rpfc") for _, k in C.COMBOS[:len(C.SMALL)])
    assert ("case3120sp", ("rpfc",)) in C.COMBOS and ("case_ACTIVSg10k", ("apfc",)) in C.COMBOS
    assert "case9241pegase" not in C.LARGE
    assert C.parse_combos(C.combos_text(C.COMBOS)) == [(c, tuple(k)) for c, k in C.COMBOS]


@pytest.mark.skipif(not __import__("shutil").which("octave"), reason="octave not available")
def test_renumbering_of_the_export(tmp_path):
    # prepare_case.m on a copy of case9 with every bus number shifted by +100: the busmap is a bijection onto 1..nb
    # (old numbers 101..109), and the prepared case has the same power flow after mapping back (neutrality_check.m)
    import subprocess, csv
    script = ("define_constants; m = loadcase('case9'); m.bus(:, BUS_I) = m.bus(:, BUS_I) + 100; "
              "m.gen(:, GEN_BUS) = m.gen(:, GEN_BUS) + 100; m.branch(:, [F_BUS T_BUS]) = m.branch(:, [F_BUS T_BUS]) + 100; "
              "savecase('case9_shift.m', m); prepare_case('case9_shift', 'prepared.mat', 'busmap.csv'); "
              "neutrality_check('case9_shift', 'prepared.mat', 'busmap.csv', 'neutrality.csv');")
    (tmp_path / "neutrality.csv").write_text("case,max_abs_vm,max_abs_va_deg,max_abs_pf,max_abs_qf,max_abs_pt,max_abs_qt,"
                                            "iterations_original,iterations_prepared,success_original,success_prepared,busmap_bijection,maps_back\n")
    subprocess.run(["octave", "--no-gui", "--eval", f"addpath('{os.path.dirname(HERE)}'); {script}"], cwd=tmp_path,
                   check=True, capture_output=True)
    bm = list(csv.DictReader(open(tmp_path / "busmap.csv")))
    assert sorted(int(r["new_bus"]) for r in bm) == list(range(1, 10))
    assert sorted(int(r["old_bus"]) for r in bm) == list(range(101, 110))
    r = list(csv.DictReader(open(tmp_path / "neutrality.csv")))[-1]
    assert float(r["max_abs_vm"]) <= 1e-9 and float(r["max_abs_va_deg"]) <= 1e-9
    assert max(float(r[k]) for k in ("max_abs_pf", "max_abs_qf", "max_abs_pt", "max_abs_qt")) <= 1e-6
    assert r["iterations_original"] == r["iterations_prepared"] and r["busmap_bijection"] == "1" and r["maps_back"] == "1"
