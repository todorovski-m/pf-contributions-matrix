"""Summarise a results folder of timing.py: per-cell statistics, the phase fractions of the port, the ratio old/port of
the medians (for every cell where both finished, slow cells included), the dense table size and the number of table
columns, and the slopes of log10(median time) against log10(n_e) per method, checks and kind:
  - local slopes between consecutive cases (ordered by n_e) with n_e >= 500, each with its two cases (T14e);
  - one least-squares slope over the cases with n_e >= 500, with the number of points (T14e);
  - the earlier global slope over all cases, labelled "all cases, dominated by fixed overhead".
A median whose cell (or, for folders without per-cell flags, whose case) is not quiet is marked.
T14f: the same slopes against the dense table size (rows x columns from the guard; for folders without the column
count, the estimated dense bytes / 8), in summary_slopes.csv and summary_local_slopes_cells.csv.
Writes summary.md, summary.csv, summary_slopes.csv, summary_local_slopes.csv and summary_local_slopes_cells.csv. When the environment
record says smoke, every output carries the label "smoke, not for publication".

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
    gp = os.path.join(folder, "guard.csv")
    guard = list(csv.DictReader(open(gp))) if os.path.exists(gp) else []
    return env, total, phase, cells, guard


def loglog_slope(points):
    """Least-squares slope of log10(median) against log10(n_e); None for fewer than two distinct n_e."""
    if len({p[0] for p in points}) < 2:
        return None
    x = np.log10([p[0] for p in points]); y = np.log10([p[1] for p in points])
    return float(np.polyfit(x, y, 1)[0])


NE_MIN = 500


def base_status(st):
    """The status without the flags: 'ok (slow)' and 'ok (predicted stall, but finished)' -> 'ok'."""
    for flag in (" (slow)", " (predicted stall, but finished)"):
        if st.endswith(flag):
            st = st[:-len(flag)]
    return st


def local_slopes(points):
    """Slopes between consecutive points (n_e, median, case, quiet) ordered by n_e, for n_e >= NE_MIN; equal n_e skipped."""
    pts = sorted(p for p in points if p[0] >= NE_MIN)
    out = []
    for a, b in zip(pts, pts[1:]):
        if b[0] == a[0]:
            continue
        out.append({"case_from": a[2], "case_to": b[2], "n_e_from": a[0], "n_e_to": b[0],
                    "slope": float(np.log10(b[1] / a[1]) / np.log10(b[0] / a[0])),
                    "not_quiet": ", ".join(p[2] for p in (a, b) if p[3] is False)})
    return out


def local_slopes_x(points):
    """Slopes between consecutive points (x, median, case, quiet) ordered by x (already filtered); equal x skipped."""
    pts = sorted(points)
    return [{"case_from": a[2], "case_to": b[2], "x_from": a[0], "x_to": b[0],
             "slope": float(np.log10(b[1] / a[1]) / np.log10(b[0] / a[0])),
             "not_quiet": ", ".join(p[2] for p in (a, b) if p[3] is False)} for a, b in zip(pts, pts[1:]) if b[0] != a[0]]


def summarize(folder):
    env, total, phase, cells, guard = load(folder)
    smoke = bool(env.get("smoke"))
    label = SMOKE_LABEL if smoke else "series"
    # a slow cell (slow-cell rule: a single run above the threshold, 3 repetitions) is flagged in its status
    status = {(c["case"], c["method"], c["kind"], c["checks"]): c["status"] + (" (slow)" if c.get("slow") == "True" else "")
              for c in cells}
    n_e = {(c, k): env["cases"][c].get(f"n_e_{k}") for c in env["cases"] for k in ("apfc", "rpfc")}
    dense = {(c, k): env["cases"][c].get(f"dense_bytes_{k}") for c in env["cases"] for k in ("apfc", "rpfc")}
    ncols = {(r["case"], r["kind"]): (int(r["columns"]) if r.get("columns") else None) for r in guard}
    nrows = {(r["case"], r["kind"]): (int(r["rows"]) if r.get("rows") else None) for r in guard}

    def table_cells(case, kind):
        """The size of the dense result table: rows x columns from the guard (T14f), else, for folders without the
        column count (the first series), the estimated dense bytes of environment.json divided by 8 (an estimate of
        all dense arrays of one run, not only the table)."""
        r, c = nrows.get((case, kind)), ncols.get((case, kind))
        if r and c:
            return r * c, "table rows x columns"
        b = dense.get((case, kind))
        return (b / 8, "dense bytes / 8") if b else (None, "")
    # quiet per cell when the folder has the T14e flags, else the flag of the case
    cq = {(c["case"], c["method"], c["kind"], c["checks"]): (c["quiet"] == "True") for c in cells if c.get("quiet") not in (None, "")}
    rows = []
    for key in sorted(set(status) | set(total), key=lambda k: (k[0], k[2], k[1], k[3])):
        case, method, kind, checks = key
        if method == "all":
            continue
        s = stats(total.get(key, []))
        row = {"label": label, "case": case, "method": method, "kind": kind, "checks": checks,
               "status": status.get(key, ""), "n_e": n_e.get((case, kind)),
               "quiet": cq.get(key, env["cases"].get(case, {}).get("quiet")),
               "dense_table_bytes": dense.get((case, kind)), "table_columns": ncols.get((case, kind)),
               "table_cells": table_cells(case, kind)[0], "table_cells_source": table_cells(case, kind)[1],
               **{f"{k}_ns": s[k] for k in ("min", "median", "mean", "std", "iqr")}, "n": s["n"]}
        if method == "port" and s["n"]:
            for p in PORT_PHASES:
                v = phase[key].get(p)
                row[f"fraction_{p}"] = (float(np.median(v)) / s["median"]) if v else None
        rows.append(row)
    med = {(r["case"], r["method"], r["kind"], r["checks"]): r["median_ns"] for r in rows if base_status(r["status"]) == "ok"}
    for r in rows:
        if r["method"] == "old" and base_status(r["status"]) == "ok":
            p = med.get((r["case"], "port", r["kind"], "off"))
            r["ratio_old_over_port_off"] = (r["median_ns"] / p) if p else None
    slopes, locs, locs_cells = [], [], []
    for method, checks in (("port", "off"), ("port", "on"), ("old", "")):
        for kind in ("apfc", "rpfc"):
            # the global slope over all cases: the cells with status exactly ok, as before T14e (reproducible)
            pts = sorted((r["n_e"], r["median_ns"], r["case"]) for r in rows
                         if r["method"] == method and r["checks"] == checks and r["kind"] == kind and r["status"] == "ok" and r["n_e"])
            # the large end: every finished cell (slow ones included) with n_e >= NE_MIN
            big = sorted((r["n_e"], r["median_ns"], r["case"], r["quiet"]) for r in rows
                         if r["method"] == method and r["checks"] == checks and r["kind"] == kind
                         and base_status(r["status"]) == "ok" and r["n_e"] and r["n_e"] >= NE_MIN)
            slopes.append({"label": label, "method": method, "checks": checks, "kind": kind,
                           "slope": loglog_slope(pts), "cases": " ".join(p[2] for p in pts),
                           "slope_ne_ge_500": loglog_slope(big), "points_ne_ge_500": len(big),
                           "cases_ne_ge_500": " ".join(p[2] for p in big),
                           "not_quiet_ne_ge_500": " ".join(p[2] for p in big if p[3] is False)})
            for l in local_slopes(big):
                locs.append({"label": label, "method": method, "checks": checks, "kind": kind, **l})
            # against the dense table size (T14f): the same cells (n_e >= NE_MIN, finished), ordered by table cells
            bigc = sorted((r["table_cells"], r["median_ns"], r["case"], r["quiet"]) for r in rows
                          if r["method"] == method and r["checks"] == checks and r["kind"] == kind
                          and base_status(r["status"]) == "ok" and r["n_e"] and r["n_e"] >= NE_MIN and r["table_cells"])
            slopes[-1].update({"slope_cells_ne_ge_500": loglog_slope(bigc), "points_cells_ne_ge_500": len(bigc),
                               "cells_source": " ".join(sorted({r["table_cells_source"] for r in rows if r["table_cells_source"]}))})
            for l in local_slopes_x(bigc):
                locs_cells.append({"label": label, "method": method, "checks": checks, "kind": kind, **l})
    cols = ["label", "case", "method", "kind", "checks", "status", "n_e", "quiet", "n", "min_ns", "median_ns", "mean_ns", "std_ns",
            "iqr_ns", "ratio_old_over_port_off"] + [f"fraction_{p}" for p in PORT_PHASES] + ["dense_table_bytes", "table_columns",
                                                                                              "table_cells", "table_cells_source"]
    with open(os.path.join(folder, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    with open(os.path.join(folder, "summary_slopes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "method", "checks", "kind", "slope", "cases", "slope_ne_ge_500",
                                           "points_ne_ge_500", "cases_ne_ge_500", "not_quiet_ne_ge_500",
                                           "slope_cells_ne_ge_500", "points_cells_ne_ge_500", "cells_source"]); w.writeheader()
        for s in slopes:
            w.writerow(s)
    with open(os.path.join(folder, "summary_local_slopes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "method", "checks", "kind", "case_from", "case_to", "n_e_from", "n_e_to",
                                           "slope", "not_quiet"]); w.writeheader()
        for l in locs:
            w.writerow(l)
    with open(os.path.join(folder, "summary_local_slopes_cells.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label", "method", "checks", "kind", "case_from", "case_to", "x_from", "x_to",
                                           "slope", "not_quiet"]); w.writeheader()
        for l in locs_cells:
            w.writerow(l)
    md = [f"**{label}**" if smoke else "", "",
          f"Environment: {env.get('cpu_model')}, {env.get('cores_logical')} cores, {env.get('os')}; git {env.get('git_matrix')}; "
          f"date {env.get('date_local')}.", "",
          f"| case | method | kind | checks | status | n_e | quiet | n | min (ns) | median (ns) | mean (ns) | std (ns) | IQR (ns) | old/port median (checks off) | {' | '.join('fraction ' + p for p in PORT_PHASES)} | dense table bytes | table columns |",
          "|" + "---|" * (16 + len(PORT_PHASES))]
    for r in rows:
        md.append("| " + " | ".join(str(x) if not isinstance(x, float) else g(x) for x in
                                     [r["case"], r["method"], r["kind"], r["checks"] or "-", r["status"], r["n_e"], r["quiet"], r["n"],
                                      r["min_ns"], r["median_ns"], r["mean_ns"], r["std_ns"], r["iqr_ns"],
                                      r.get("ratio_old_over_port_off")] + [r.get(f"fraction_{p}") for p in PORT_PHASES]
                                     + [r.get("dense_table_bytes"), r.get("table_columns")]).replace("None", "") + " |")
    md += ["", f"Slope of log10(median time) against log10(n_e): one least-squares fit over the cases with n_e >= {NE_MIN} "
           "(finished cells, slow ones included), and, separately, the fit over all cases (cells with status ok), which "
           "is dominated by the fixed overhead of the small cases:", "",
           f"| method | checks | kind | slope, n_e >= {NE_MIN} | points | cases (n_e >= {NE_MIN}) | not quiet among them | "
           "slope, all cases, dominated by fixed overhead | cases (all) |", "|---|---|---|---|---|---|---|---|---|"]
    for s in slopes:
        md.append(f"| {s['method']} | {s['checks'] or '-'} | {s['kind']} | {g(s['slope_ne_ge_500'])} | {s['points_ne_ge_500']} | "
                  f"{s['cases_ne_ge_500']} | {s['not_quiet_ne_ge_500']} | {g(s['slope'])} | {s['cases']} |")
    md += ["", f"Local slopes between consecutive cases (ordered by n_e, n_e >= {NE_MIN}):", "",
           "| method | checks | kind | from case (n_e) | to case (n_e) | local slope | not quiet |", "|---|---|---|---|---|---|---|"]
    for l in locs:
        md.append(f"| {l['method']} | {l['checks'] or '-'} | {l['kind']} | {l['case_from']} ({l['n_e_from']}) | "
                  f"{l['case_to']} ({l['n_e_to']}) | {g(l['slope'])} | {l['not_quiet']} |")
    md += ["", f"Slope of log10(median time) against log10(dense table cells), the same cells (n_e >= {NE_MIN}, finished); "
           "table cells: " + (slopes[0]["cells_source"] if slopes else "") + ":", "",
           "| method | checks | kind | slope against table cells | points |", "|---|---|---|---|---|"]
    for s in slopes:
        md.append(f"| {s['method']} | {s['checks'] or '-'} | {s['kind']} | {g(s['slope_cells_ne_ge_500'])} | {s['points_cells_ne_ge_500']} |")
    md += ["", "Local slopes against the table cells between consecutive cases (ordered by table cells, n_e >= "
           f"{NE_MIN}):", "", "| method | checks | kind | from case (cells) | to case (cells) | local slope | not quiet |",
           "|---|---|---|---|---|---|---|"]
    for l in locs_cells:
        md.append(f"| {l['method']} | {l['checks'] or '-'} | {l['kind']} | {l['case_from']} ({g(float(l['x_from']))}) | "
                  f"{l['case_to']} ({g(float(l['x_to']))}) | {g(l['slope'])} | {l['not_quiet']} |")
    if smoke:
        md.insert(1, f"All numbers in this table are from a {SMOKE_LABEL} run.")
    open(os.path.join(folder, "summary.md"), "w").write("\n".join(md) + "\n")
    return rows, slopes, locs


if __name__ == "__main__":
    summarize(sys.argv[1])
    print(f"wrote summary.md, summary.csv, summary_slopes.csv, summary_local_slopes.csv and summary_local_slopes_cells.csv in {sys.argv[1]}")
