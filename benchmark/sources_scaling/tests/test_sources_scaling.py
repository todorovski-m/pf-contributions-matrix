"""Machinery tests of benchmark/sources_scaling/sources_scaling.py: subset selection, k values, the correctness of the
k-column computation, extrapolation labels and the CSV schema. No time is checked or recorded.
Run: python -B -m pytest -p no:cacheprovider benchmark/sources_scaling/tests"""
import csv
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import sources_scaling as S   # noqa: E402


def test_k_values():
    assert S.k_values(1) == [1]
    assert S.k_values(8) == [1, 2, 4, 8]
    assert S.k_values(19) == [1, 2, 4, 8, 16, 19]
    assert S.k_values(19, kmax=4) == [1, 2, 4]


def test_subsets_deterministic_and_nested():
    p1, p2 = S.subsets(71), S.subsets(71)
    assert np.array_equal(p1, p2) and sorted(p1.tolist()) == list(range(71))
    for k in (1, 2, 4, 8):
        assert set(p1[:k]) <= set(p1[:2 * k])
    assert not np.array_equal(S.subsets(71, seed=1), p1)


def test_k_columns_equal_stack_and_public_table():
    res = S.check(["case9"], [1, 2, 4])
    assert res and all(r["ok"] for r in res)
    assert {(r["kind"], r["k"]) for r in res} == {("apfc", 1), ("apfc", 2), ("rpfc", 1), ("rpfc", 2), ("rpfc", 4)}


def test_run_schema_and_extrapolation_labels(tmp_path, monkeypatch):
    # with the repeated mode measured only up to k = 1, the repeated rows for k = 2 and 4 are extrapolated
    monkeypatch.setattr(S, "KMAX_REPEATED", 1)
    out = tmp_path / "res"
    S.main(["run", "--out", str(out), "--cases", "case9", "--kmax", "4", "--budget", "0.01", "--min-rep", "2",
            "--max-rep", "2", "--warmup", "0", "--quiet-wait", "0", "--smoke"])
    rows = list(csv.DictReader(open(out / "sources.csv")))
    assert list(rows[0].keys()) == S.CSV_COLUMNS
    assert all(r["label"] == "smoke, not for publication" for r in rows)
    rep = {(r["kind"], int(r["k"])): r for r in rows if r["mode"] == "repeated"}
    one = {(r["kind"], int(r["k"])) for r in rows if r["mode"] == "one run"}
    assert rep[("rpfc", 1)]["extrapolated"] == "False"
    for k in (2, 4):
        r = rep[("rpfc", k)]
        assert r["extrapolated"] == "True"
        assert float(r["median_ns"]) == k * float(rep[("rpfc", 1)]["median_ns"])
    # case9: 3 APFC sources (k = 1, 2, 3) and 8 RPFC sources (k = 1, 2, 4 with --kmax 4)
    assert one == {("apfc", 1), ("apfc", 2), ("apfc", 3), ("rpfc", 1), ("rpfc", 2), ("rpfc", 4)}
    for f in ("summary.md", "environment.json", "sources_time.pdf", "sources_ratio.pdf"):
        assert (out / f).exists()
