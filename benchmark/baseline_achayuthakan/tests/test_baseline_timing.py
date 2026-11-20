"""Machinery tests of baseline_timing.py (brief T14h, step 4). No time is measured: timing.measure_cell and the quiet
wait are replaced by stubs that return fixed samples, so the tests check the guard, the cells, the CSV schema, the
labels and the summary only. Run: python -B -m pytest -p no:cacheprovider benchmark/baseline_achayuthakan/tests"""
import csv
import os
import sys

import pytest

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import baseline_timing as BT   # noqa: E402


def _stub(monkeypatch, calls):
    def measure(fn, *a, **k):
        res, _ = fn()                            # the timed function runs once, to check that it works
        calls.append(res.frame.shape if hasattr(res, "frame") else res.shape)   # read-only returns BranchParams
        return {"status": "ok", "samples": [1000, 2000, 3000], "phases": [None] * 3, "warmup": 0, "slow": False}
    monkeypatch.setattr(BT.T, "measure_cell", measure)
    monkeypatch.setattr(BT.T, "wait_quiet", lambda *a, **k: {"quiet": True, "load1_end": 0.5, "wait_s": 0.0,
                                                              "load1_start": 0.5, "threshold": 2.0, "top_processes": []})


@pytest.mark.parametrize("name", ["case9", "case118", "case_ACTIVSg2000"])
def test_guard_on_real_data(name):
    # real exports, no stub; case_ACTIVSg2000 has terminal-demand rows with a repeated label (C2_1088), which the
    # guard of the first version could not reindex
    g = BT.guard(BT.CMP.load(name))
    assert g["ok"] and g["net_vs_apfc"] <= 1e-9 and g["rows_missing"] == 0 and g["inverse_vs_sparse"] <= 1e-10


def test_guard_only_prints_all_networks(capsys):
    assert BT.main(["run", "--guard-only", "--cases", "case1,case9"]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "| " + " | ".join(BT.GUARD_COLUMNS) + " |" and len(out) == 4
    assert out[2].startswith("| case1 |") and out[3].startswith("| case9 |") and all(l.endswith("| True |") for l in out[2:])


def test_run_schema_labels_and_summary(tmp_path, monkeypatch):
    calls = []
    _stub(monkeypatch, calls)
    out = tmp_path / "res"
    BT.main(["run", "--out", str(out), "--cases", "case1,case9", "--smoke"])
    rows = list(csv.DictReader(open(out / "baseline_timing.csv")))
    assert list(rows[0].keys()) == BT.CSV_COLUMNS
    assert [(r["network"], r["method"]) for r in rows] == [(c, m) for c in ("case1", "case9") for m in BT.METHODS]
    assert all(r["label"] == "smoke, not for publication" and r["status"] == "ok" and r["median_ns"] == "2000.0" for r in rows)
    assert len(calls) == 2 * len(BT.METHODS) == 14 and all(s[0] > 0 for s in calls)
    gd = list(csv.DictReader(open(out / "guard.csv")))
    assert [g["ok"] for g in gd] == ["True", "True"]
    s = (out / "summary.md").read_text()
    assert s.startswith("**smoke, not for publication**") and "| case9 | net-sparse | ok |" in s
    assert "| case9 | port-demand | ok | 3 | 2000.0 | 1000.0 | 1000.0 | False | True | 1 | 1 |" in s


def test_guard_failure_skips_the_network(tmp_path, monkeypatch):
    calls = []
    _stub(monkeypatch, calls)
    monkeypatch.setattr(BT, "guard", lambda data: {"net_vs_apfc": 1.0, "inverse_vs_sparse": 0.0, "ok": False})
    out = tmp_path / "res"
    BT.main(["run", "--out", str(out), "--cases", "case9", "--smoke"])
    rows = list(csv.DictReader(open(out / "baseline_timing.csv")))
    assert all(r["status"] == "not timed (guard failed)" for r in rows) and calls == []


def test_estimate():
    s = BT.estimate(BT.CASES, 20, 300)
    assert "| cells (networks x 7 cells) | 49 |" in s and "| lower bound: cells x budget of 20 s (s) | 980 |" in s


@pytest.mark.parametrize("name", ["case9", "case118"])
def test_port_cells_on_real_data(name, monkeypatch):
    """The three port-side cells on the real exports (no stub of the computation): built with checks off and ending at
    the shapes of reports/T14h_scope_note.md. Checks off is shown directly: the condition estimate and the attenuation
    range check raise if called, and the cycle check must be called with warn_cyclic=False."""
    import pfc_matrix.checks as K
    from pfc_matrix import Case, Options, apfc, read_branch_params
    data = BT.CMP.load(name)
    reference_full = apfc(Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"],
                                           params=read_branch_params(data["path"])), Options(checks=False))[0]
    reference_demand = BT.CMP.apfc_bus_table(data)[0]          # the guard's table (checks on); values do not depend on checks

    def forbidden(*a, **k):
        raise AssertionError("a check that checks=False switches off was called")
    cycle = K.check_cycle
    seen = []

    def cycle_spy(*a, **k):
        seen.append(k.get("warn_cyclic", True))
        return cycle(*a, **k)
    monkeypatch.setattr(K, "condest", forbidden)
    monkeypatch.setattr(K, "check_attenuation", forbidden)
    monkeypatch.setattr(K, "check_cycle", cycle_spy)

    full, ph_full = BT.runner(data, "port-full")()
    demand, ph_dem = BT.runner(data, "port-demand")()
    params, ph_read = BT.runner(data, "read-only")()
    assert seen == [False, False]
    assert "condest" not in ph_full and "condest" not in ph_dem and ph_read is None
    assert set(ph_full) == {"input", "prepare", "assemble", "solve", "table"}

    # port-full: APFC's full table (demand rows and both end rows of every branch), as apfc with checks off
    assert list(full.index) == list(reference_full.index) and list(full.columns) == list(reference_full.columns)
    assert (full.to_numpy() == reference_full.to_numpy()).all()
    assert any(r.startswith("L_{f.") for r in full.index) and any(r.startswith("L_{t.") for r in full.index)
    # port-demand: the demand rows D_<bus>, generator columns G_<bus> summed per bus, as compare.apfc_bus_table
    assert list(demand.index) == list(reference_demand.index) and list(demand.columns) == list(reference_demand.columns)
    assert abs(demand.to_numpy() - reference_demand.to_numpy()).max() <= 1e-12
    assert all(r.startswith("D_") for r in demand.index) and all(c.startswith("G_") for c in demand.columns)
    # the baseline cells keep their rows; in these networks they are exactly port-demand's rows and columns
    net = BT.runner(data, "net-sparse")()[0]
    assert set(net.index) == set(demand.index) and set(net.columns) == set(demand.columns)
    # read-only: the branch parameters, one row per branch of the export
    assert params.frame.shape == (len(data["branch"]), 5) and params.base_mva > 0
