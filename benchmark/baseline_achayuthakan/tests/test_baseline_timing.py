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
        table, _ = fn()                          # the timed function runs once, to check that it works
        calls.append(table.shape)
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
    assert len(calls) == 10 and all(s[0] > 0 for s in calls)
    gd = list(csv.DictReader(open(out / "guard.csv")))
    assert [g["ok"] for g in gd] == ["True", "True"]
    s = (out / "summary.md").read_text()
    assert s.startswith("**smoke, not for publication**") and "| case9 | net-sparse | ok |" in s


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
    assert "| cells (networks x 5 methods) | 35 |" in s and "| lower bound: cells x budget of 20 s (s) | 700 |" in s
