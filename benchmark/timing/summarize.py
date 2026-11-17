"""Summarise a results folder of timing.py: per-cell statistics, the phase fractions of the port, the ratio old/port of
the medians, and, for the cases ordered by n_e, the slope of a log-log fit of the median time against n_e for each
method and kind (with the cases used). Writes summary.md and summary.csv into the folder. When the environment record
says smoke, every output carries the label "smoke, not for publication".

Usage: python -B benchmark/timing/summarize.py <results folder>
"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from timing import stats, SMOKE_LABEL, PORT_PHASES   # noqa: E402

g = lambda x: "" if x is None else f"{x:.17g}"


def load(folder):
    env = json.load(open(os.path.join(folder, "environment.json")))
    total, phase = defaultdict(list), defaultdict(lambda: defaultdict(list))
    with open(os.path.join(folder, "samples.csv")) as fh:
        for r in csv.DictReader(fh):
            key = (r["case"], r["method"], r["kind"], r["checks"])
            if r["phase"]:
                phase[key][r["phase"]].append(int(r["ns"]))
            else:
                total[key].append(int(r["ns"]))
    cells = list(csv.DictReader(open(os.path.join(folder, "cells.csv"))))
    return env, total, phase, cells


def loglog_slope(points):
    """Least-squares slope of log10(median) against log10(n_e); None for fewer than two distinct n_e."""
    if len({p[0] for p in points}) < 2:
        return None
    x = np.log10([p[0] for p in points]); y = np.log10([p[1] for p in points])
    return float(np.polyfit(x, y, 1)[0])


def summarize(folder):
    env, total, phase, cells = load(folder)
    smoke = bool(env.get("smoke"))
    label = SMOKE_LABEL if smoke else "series"
    status = {(c["case"], c["method"], c["kind"], c["checks"]): c["status"] for c in cells}
    n_e = {(c, k): env["cases"][c].get(f"n_e_{k}") for c in env["cases"] for k in ("apfc", "rpfc")}
    rows = []
    for key in sorted(set(status) | set(total), key=lambda k: (k[0], k[2], k[1], k[3])):
        case, method, kind, checks = key
        if method == "all":
            continue
        s = stats(total.get(key, []))
        row = {"label": label, "case": case, "method": method, "kind": kind, "checks": checks,
               "status": status.get(key, ""), "n_e": n_e.get((case, kind)), "quiet": env["cases"].get(case, {}).get("quiet"),
               **{f"{k}_ns": s[k] for k in ("min", "median", "mean", "std", "iqr")}, "n": s["n"]}
        if method == "port" and s["n"]:
            for p in PORT_PHASES:
                v = phase[key].get(p)
                row[f"fraction_{p}"] = (float(np.median(v)) / s["median"]) if v else None
        rows.append(row)
    med = {(r["case"], r["method"], r["kind"], r["checks"]): r["median_ns"] for r in rows if r["status"] == "ok"}
    for r in rows:
        if r["method"] == "old" and r["status"] == "ok":
            p = med.get((r["case"], "port", r["kind"], "off"))
            r["ratio_old_over_port_off"] = (r["median_ns"] / p) if p else None
    slopes = []
    for method, checks in (("port", "off"), ("port", "on"), ("old", "")):
        for kind in ("apfc", "rpfc"):
            pts = sorted((r["n_e"], r["median_ns"], r["case"]) for r in rows
                         if r["method"] == method and r["checks"] == checks and r["kind"] == kind and r["status"] == "ok" and r["n_e"])
            slopes.append({"label": label, "method": method, "checks": checks, "kind": kind,
                           "slope": loglog_slope(pts), "cases": " ".join(p[2] for p in pts)})
    cols = ["label", "case", "method", "kind", "checks", "status", "n_e", "quiet", "n", "min_ns", "median_ns", "mean_ns", "std_ns",
            "iqr_ns", "ratio_old_over_port_off"] + [f"fraction_{p}" for p in PORT_PHASES]
    with open(os.path.join(folder, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    with open(os.path.join(folder, "summary_slopes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "method", "checks", "kind", "slope", "cases"]); w.writeheader()
        for s in slopes:
            w.writerow(s)
    md = [f"**{label}**" if smoke else "", "",
          f"Environment: {env.get('cpu_model')}, {env.get('cores_logical')} cores, {env.get('os')}; git {env.get('git_matrix')}; "
          f"date {env.get('date_local')}.", "",
          f"| case | method | kind | checks | status | n_e | quiet | n | min (ns) | median (ns) | mean (ns) | std (ns) | IQR (ns) | old/port median (checks off) | {' | '.join('fraction ' + p for p in PORT_PHASES)} |",
          "|" + "---|" * (14 + len(PORT_PHASES))]
    for r in rows:
        md.append("| " + " | ".join(str(x) if not isinstance(x, float) else g(x) for x in
                                     [r["case"], r["method"], r["kind"], r["checks"] or "-", r["status"], r["n_e"], r["quiet"], r["n"],
                                      r["min_ns"], r["median_ns"], r["mean_ns"], r["std_ns"], r["iqr_ns"],
                                      r.get("ratio_old_over_port_off")] + [r.get(f"fraction_{p}") for p in PORT_PHASES]).replace("None", "") + " |")
    md += ["", "Slope of log10(median time) against log10(n_e), cases ordered by n_e:", "",
           "| method | checks | kind | slope | cases used |", "|---|---|---|---|---|"]
    for s in slopes:
        md.append(f"| {s['method']} | {s['checks'] or '-'} | {s['kind']} | {g(s['slope'])} | {s['cases']} |")
    if smoke:
        md.insert(1, f"All numbers in this table are from a {SMOKE_LABEL} run.")
    open(os.path.join(folder, "summary.md"), "w").write("\n".join(md) + "\n")
    return rows, slopes


if __name__ == "__main__":
    summarize(sys.argv[1])
    print(f"wrote summary.md, summary.csv and summary_slopes.csv in {sys.argv[1]}")
