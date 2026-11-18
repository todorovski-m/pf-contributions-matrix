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


@pytest.mark.skipif(not os.environ.get("PFC_SYN_INPUTS"), reason="needs the synthetic exports (PFC_SYN_INPUTS)")
def test_guard_flows_with_omitted_zero_terminals():
    # a reactive Case-3 branch with a zero terminal has one row only (T15b; syn9_dead_charge_to): the guard's flows follow the table rows
    data = C.load("syn9_dead_charge_to", [os.environ["PFC_SYN_INPUTS"]])
    table, pub, flows, _ = T.port_flows(data, "rpfc")
    assert len(flows) == len(table)
    assert np.abs(table.to_numpy().sum(axis=1) - flows).max() <= 1e-10


def _mock_old(monkeypatch, done, max_cell, unmapped, rowsum):
    """The old method and the comparison are mocked: the given finish flag, largest cell difference, unmapped port
    rows and old row-sum deviation."""
    import correctness
    monkeypatch.setattr(T, "old_with_iterations", lambda *a, **k: (None, None, done, 10000 if not done else 7))
    monkeypatch.setattr(correctness, "compare", lambda *a, **k: {
        "unmapped_port_rows": unmapped, "max_cell": max_cell, "port_rowsum_vs_old_flow": rowsum,
        "excluded_rows": 0, "unmapped_old_columns": [], "unmapped_port_columns": []})


def test_guard_stalled_old_method_checks_the_port_only(monkeypatch):
    # T14e: an unfinished old table (unmapped rows, row-sum deviation 267) is not compared; the port is checked
    data = C.load("case9", C.input_dirs())
    branch_old, _, _, _ = C.correct_qsf(data)
    _mock_old(monkeypatch, done=False, max_cell=300.0, unmapped=["Tf_1", "Tt_2"], rowsum=267.0)
    g = T.guard(data, "rpfc", 10000, 60, branch_old)
    assert g["ok"] and g["old_status"] == "did not finish (iteration cap)" and g["iterations"] == 10000
    assert g["guard_checks"] == "port only (old did not finish)"
    assert "old_rowsum" not in g and "unmapped_matrix_rows" not in g
    assert g["port_equal_public"] and g["port_rowsum"] <= 1e-10


def test_guard_finished_old_method_with_a_cell_difference_fails(monkeypatch):
    data = C.load("case9", C.input_dirs())
    branch_old, _, _, _ = C.correct_qsf(data)
    _mock_old(monkeypatch, done=True, max_cell=2e-10, unmapped=[], rowsum=1e-12)  # only the cell differs
    g = T.guard(data, "rpfc", 10000, 60, branch_old)
    assert not g["ok"] and g["old_status"] == "ok" and g["guard_checks"] == "port and old (old finished)"
    assert g["max_cell"] == 2e-10


class _FakeClock:
    def __init__(self):
        self.t = 0.0
        self.slept = []

    def now(self):
        return self.t

    def sleep(self, s):
        self.slept.append(s); self.t += s


def test_quiet_wait_waits_then_continues(monkeypatch):
    clk = _FakeClock(); loads = iter([3.0, 2.5, 1.0])          # 8 cores: threshold 2.0
    monkeypatch.setattr(T, "CLOCK_FN", clk.now); monkeypatch.setattr(T, "SLEEP_FN", clk.sleep)
    monkeypatch.setattr(T, "LOAD_FN", lambda: next(loads))
    q = T.wait_quiet(8, max_wait_s=300, poll_s=10)
    assert q["quiet"] and q["wait_s"] == 20 and clk.slept == [10, 10]
    assert q["load1_start"] == 3.0 and q["load1_end"] == 1.0 and q["top_processes"] == []


def test_quiet_wait_expires_and_flags(monkeypatch):
    clk = _FakeClock()
    monkeypatch.setattr(T, "CLOCK_FN", clk.now); monkeypatch.setattr(T, "SLEEP_FN", clk.sleep)
    monkeypatch.setattr(T, "LOAD_FN", lambda: 5.0)
    monkeypatch.setattr(T, "top_processes", lambda n=3: ["1 99.0 a", "2 50.0 b", "3 10.0 c"])
    q = T.wait_quiet(8, max_wait_s=30, poll_s=10)
    assert not q["quiet"] and q["wait_s"] == 30 and sum(clk.slept) == 30
    assert q["top_processes"] == ["1 99.0 a", "2 50.0 b", "3 10.0 c"]


def test_quiet_wait_is_outside_the_samples(tmp_path, monkeypatch):
    # a real wait of about 0.6 s before each cell (load high for three polls of 0.2 s), then a tiny series of case9
    # APFC: every sample is far below the wait, the wait is recorded per cell, and every cell is quiet
    calls = {"n": 0}

    def load():
        calls["n"] += 1
        return 99.0 if calls["n"] % 4 != 0 else 0.0
    monkeypatch.setattr(T, "LOAD_FN", load)
    out = tmp_path / "res"
    T.main(["run", "--out", str(out), "--combos", "case9:apfc", "--budget", "0.05", "--min-rep", "2", "--max-rep", "3",
            "--warmup", "1", "--no-memory", "--smoke", "--quiet-wait", "5", "--quiet-poll", "0.2"])
    import pandas as pd
    cells = pd.read_csv(out / "cells.csv")
    timed = cells[cells["n"] > 0]
    assert len(timed) == 3 and (timed["wait_s"] >= 0.6).all() and timed["quiet"].all()
    smp = pd.read_csv(out / "samples.csv")
    assert (smp[smp["phase"].isna()]["ns"] < 0.3e9).all()


def _crafted_folder(path):
    """A results folder with known medians: port off 1, 4, 16 ms and old 100, 800, 6400 ms at n_e 500, 1000, 2000
    (local slopes 2 and 3), plus a small case at n_e 10; the old cell at n_e 2000 is slow and not quiet."""
    import csv as _csv, json as _json
    os.makedirs(path, exist_ok=True)
    cases = {"s": 10, "a": 500, "b": 1000, "c": 2000}
    port = {"s": 1e6, "a": 1e6, "b": 4e6, "c": 16e6}; old = {"s": 2e6, "a": 100e6, "b": 800e6, "c": 6400e6}
    env = {"label": "series", "smoke": False, "args": {"timeout": 1800.0},
           "cases": {c: {"n_e_apfc": n, "dense_bytes_apfc": 8 * n, "quiet": True} for c, n in cases.items()}}
    _json.dump(env, open(os.path.join(path, "environment.json"), "w"))
    with open(os.path.join(path, "samples.csv"), "w", newline="") as fh:
        w = _csv.writer(fh); w.writerow(["label", "case", "method", "kind", "checks", "repetition", "ns", "phase"])
        for c in cases:
            for k in range(3):
                w.writerow(["series", c, "port", "apfc", "off", k + 1, int(port[c]), ""])
                w.writerow(["series", c, "old", "apfc", "", k + 1, int(old[c]), ""])
    with open(os.path.join(path, "cells.csv"), "w", newline="") as fh:
        w = _csv.writer(fh); w.writerow(["label", "case", "checks", "kind", "method", "n", "slow", "status", "warmup", "quiet"])
        for c in cases:
            w.writerow(["series", c, "off", "apfc", "port", 3, False, "ok", 2, True])
            w.writerow(["series", c, "", "apfc", "old", 3, c == "c", "ok", 2, c != "c"])
    with open(os.path.join(path, "guard.csv"), "w", newline="") as fh:
        w = _csv.writer(fh); w.writerow(["label", "case", "kind", "columns"])
        for c in cases:
            w.writerow(["series", c, "apfc", 7])


def test_summary_ratio_slow_cells_and_slopes(tmp_path):
    import summarize as S
    f = tmp_path / "res"; _crafted_folder(str(f))
    rows, slopes, locs = S.summarize(str(f))
    r = {(x["case"], x["method"]): x for x in rows}
    assert r[("c", "old")]["status"] == "ok (slow)" and abs(r[("c", "old")]["ratio_old_over_port_off"] - 400.0) < 1e-12
    assert r[("c", "old")]["quiet"] is False and r[("c", "port")]["dense_table_bytes"] == 16000 and r[("c", "port")]["table_columns"] == 7
    port = [l for l in locs if l["method"] == "port"]
    assert [(l["case_from"], l["case_to"]) for l in port] == [("a", "b"), ("b", "c")]
    assert all(abs(l["slope"] - 2.0) < 1e-12 for l in port)
    oldl = [l for l in locs if l["method"] == "old"]
    assert all(abs(l["slope"] - 3.0) < 1e-12 for l in oldl) and oldl[-1]["not_quiet"] == "c"
    sp = {(x["method"], x["checks"], x["kind"]): x for x in slopes}
    assert sp[("port", "off", "apfc")]["points_ne_ge_500"] == 3 and abs(sp[("port", "off", "apfc")]["slope_ne_ge_500"] - 2.0) < 1e-12
    assert sp[("old", "", "apfc")]["points_ne_ge_500"] == 3 and sp[("old", "", "apfc")]["not_quiet_ne_ge_500"] == "c"
    # the global slope keeps its earlier rule (status exactly ok): the slow old cell is not in it
    assert sp[("old", "", "apfc")]["cases"] == "s a b"
    import plot_timing as P
    for out in P.plot(str(f)):
        assert os.path.getsize(out) > 0


def test_old_stall_predictor():
    # tap_loop RPFC: one 2-cycle, both buses never processed (T14d); case9 and tap_loop APFC can finish
    import old_stall_predictor as P
    dirs = C.input_dirs()
    tl = C.load("tap_loop", dirs)
    r = P.predict(tl, "rpfc", C.correct_qsf(tl)[0])
    assert r["label"] == "predicted stall" and r["predicted"] == 2 and r["cyclic_buses"] == 2 and r["orphan_buses"] == 0
    assert P.predict(tl, "apfc")["label"] == "can finish"
    d9 = C.load("case9", dirs)
    assert P.predict(d9, "rpfc")["predicted"] == 0 and P.predict(d9, "apfc")["predicted"] == 0
    # an orphan inflow: a branch with zero from-end flow and positive to-end flow is in no outflow set
    br = d9["branch"].copy(); i = br.index[0]
    br.loc[i, "QF"] = 0.0; br.loc[i, "QT"] = abs(br.loc[i, "QT"]) + 1.0
    r = P.predict(d9, "rpfc", br)
    assert r["orphan_buses"] >= 1 and r["label"] == "predicted stall"


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
