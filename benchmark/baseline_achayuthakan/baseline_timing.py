"""Runtime of the Achayuthakan et al. (2010) baseline against the port (brief T14h, step 4).

For each network of the comparison (case1, case9, case14, case30, case118, case_ACTIVSg500, case_ACTIVSg2000), the
computation from the power-flow quantities (the DataFrames of the harness's loader) is timed in seven cells
(reports/T14h_scope_note.md):
- "port-full": timing.port_runner(data, "apfc", False), exactly the port cell of the final series: checks off (no
  condition estimate, no attenuation range check, no cycle warning), the read of the branch parameters from the .mat
  file and Case.from_frames inside, ending at APFC's full table (demand rows and both end rows of every branch);
- "port-demand": the same call, followed by the row selection and per-bus summation of compare.apfc_bus_table
  (demand_table below), ending at the demand rows D_<bus> with the generator columns G_<bus> summed per bus;
- "read-only": pfc_matrix.read_branch_params alone (the file read contained in both port cells);
- "net-inverse", "net-sparse", "gross-inverse", "gross-sparse": the baseline (achayuthakan.contributions), both loss
  variants, with the explicit inverse (variant (a), as published) and with one sparse LU (variant (b)), unchanged;
  their tables keep the terminal-demand rows C2_<bus> (6 extra rows in case_ACTIVSg2000) and read no file.
The machinery of benchmark/timing/ is imported and not changed: measure_cell (warm-up, repetitions within a budget,
median, time-out, slow-cell rule), wait_quiet before each cell (outside the timed region), stats and the environment
record. A guard runs first for each network, on the same tables as compare.py (bus-demand rows D_<bus>, APFC's
generator columns summed per bus): the net variant must equal APFC within 1e-9 MW with no demand row missing, and the
two solution methods must agree within 1e-10; a network that fails it is not timed.

Usage (from the matrix root):
    python -B benchmark/baseline_achayuthakan/baseline_timing.py estimate [--budget 20]
    python -B benchmark/baseline_achayuthakan/baseline_timing.py run --out DIR [--cases ...] [--budget 20] [--smoke]
    python -B benchmark/baseline_achayuthakan/baseline_timing.py run --guard-only [--cases ...]   (guards only, no timing)
Outputs: baseline_timing.csv (label, network, method, status, n, median_ns, min_ns, iqr_ns, slow, quiet, load1),
guard.csv, environment.json, summary.md (with the ratio of each median to the port-full and port-demand medians).
"""
import argparse
import csv
import json
import os
import re
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
import timing as T          # noqa: E402  (thread settings before numpy)
import numpy as np          # noqa: E402
import pandas as pd         # noqa: E402
import compare as CMP       # noqa: E402
import achayuthakan as A    # noqa: E402

CASES = list(CMP.CASES)
METHODS = ["port-full", "port-demand", "read-only", "net-inverse", "net-sparse", "gross-inverse", "gross-sparse"]
PORT_REFS = ["port-full", "port-demand"]
CSV_COLUMNS = ["label", "network", "method", "status", "n", "median_ns", "min_ns", "iqr_ns", "slow", "quiet", "load1"]


def demand_table(t):
    """The row selection and per-bus summation of compare.apfc_bus_table (compare.py lines 55-64), applied to APFC's
    full table t: the demand rows D_<bus> (repeated rows of one bus merged) with the generator columns summed per bus.
    The test checks that it equals compare.apfc_bus_table on real data."""
    dem = t.loc[[r for r in t.index if re.fullmatch(r"D_\d+(_\d+)?", r)]]
    base = [re.match(r"(D_\d+)", r).group(1) for r in dem.index]
    dem = dem.groupby(base, sort=False).sum()
    gcols = {}
    for c in dem.columns:
        m = re.fullmatch(r"(G_\d+)(_\d+)?", c)
        if m:
            gcols.setdefault(m.group(1), []).append(c)
    return pd.DataFrame({k: dem[v].sum(axis=1) for k, v in gcols.items()}, index=dem.index)


def runner(data, method):
    """The timed function of one cell, returning (result, phases)."""
    if method == "port-full":
        return T.port_runner(data, "apfc", False)
    if method == "port-demand":
        full = T.port_runner(data, "apfc", False)

        def run():
            table, ph = full()
            return demand_table(table), ph
        return run
    if method == "read-only":
        from pfc_matrix import read_branch_params
        path = data["path"]
        return lambda: (read_branch_params(path), None)
    kind, how = method.split("-")
    return lambda: (A.contributions(data["bus"], data["gen"], data["branch"], kind=kind, method=how), None)


GUARD_COLUMNS = ["network", "net_vs_apfc", "rows_missing", "inverse_vs_sparse", "ok"]


def guard(data):
    """The tables of compare.py: APFC's demand rows D_<bus> with the generator columns summed per bus
    (CMP.apfc_bus_table), the baseline's bus-demand rows D_<bus> only (its terminal-demand rows C2_<bus> have no
    counterpart there and may repeat a label, as C2_1088 in case_ACTIVSg2000), and the difference by CMP.diff."""
    ap = CMP.apfc_bus_table(data)[0]
    bus = data["bus"].set_index("BUS_I")
    demand = np.array([float(bus.loc[int(r[2:]), "PD"]) for r in ap.index])
    net = A.contributions(data["bus"], data["gen"], data["branch"], kind="net", method="sparse")
    dd = CMP.diff(net.loc[[r for r in net.index if r.startswith("D_")]], ap, demand)
    d_meth = max(float(np.abs(A.contributions(data["bus"], data["gen"], data["branch"], kind=k, method="inverse").to_numpy()
                              - A.contributions(data["bus"], data["gen"], data["branch"], kind=k, method="sparse").to_numpy()).max())
                 for k in ("net", "gross"))
    missing = len(dd["rows_missing"])
    return {"net_vs_apfc": dd["max_abs"], "rows_missing": missing, "inverse_vs_sparse": d_meth,
            "ok": dd["max_abs"] <= 1e-9 and missing == 0 and d_meth <= 1e-10}


def guard_only(cases):
    """The guards alone (no timing): a Markdown table, one row per network."""
    L = ["| " + " | ".join(GUARD_COLUMNS) + " |", "|" + "---|" * len(GUARD_COLUMNS)]
    for c in cases:
        gd = guard(CMP.load(c)); gd["network"] = c
        L.append("| " + " | ".join(f"{gd[k]:.17g}" if isinstance(gd[k], float) else str(gd[k]) for k in GUARD_COLUMNS) + " |")
    return "\n".join(L)


def estimate(cases, budget, quiet_wait):
    cells = len(cases) * len(METHODS)
    return "\n".join(["| quantity | value |", "|---|---|", f"| networks | {len(cases)} |", f"| cells (networks x {len(METHODS)} cells) | {cells} |",
                      f"| lower bound: cells x budget of {budget:g} s (s) | {cells * budget:g} |",
                      f"| added at most by the quiet wait, cells x {quiet_wait:g} s (s) | {cells * quiet_wait:g} |", "",
                      "Not estimated (depends on run times): warm-up runs, cells whose single run exceeds the budget, the guards "
                      "and the loading of the cases."])


def run(a):
    os.makedirs(a.out, exist_ok=True)
    env = T.environment_record(a.smoke, a)
    ncpu = os.cpu_count()
    rows, guards = [], []
    for c in a.cases.split(","):
        data = CMP.load(c)
        gd = guard(data); gd["network"] = c; guards.append(gd)
        for m in METHODS:
            if not gd["ok"]:
                rows.append({"label": env["label"], "network": c, "method": m, "status": "not timed (guard failed)", "n": 0,
                             "median_ns": None, "min_ns": None, "iqr_ns": None, "slow": None, "quiet": None, "load1": None})
                continue
            q = T.wait_quiet(ncpu, a.quiet_wait, a.quiet_poll)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = T.measure_cell(runner(data, m), a.budget, a.min_rep, a.max_rep, a.warmup, a.timeout, slow_s=a.slow)
            s = T.stats(r["samples"])
            rows.append({"label": env["label"], "network": c, "method": m, "status": r["status"], "n": s["n"],
                         "median_ns": s["median"], "min_ns": s["min"], "iqr_ns": s["iqr"], "slow": r["slow"],
                         "quiet": q["quiet"], "load1": q["load1_end"]})
    with open(os.path.join(a.out, "baseline_timing.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(a.out, "guard.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=GUARD_COLUMNS, lineterminator="\n"); w.writeheader()
        for r in guards:
            w.writerow(r)
    with open(os.path.join(a.out, "environment.json"), "w") as fh:
        json.dump(env, fh, indent=1, default=str)
    summarize(a.out)
    return rows


def summarize(folder):
    rows = list(csv.DictReader(open(os.path.join(folder, "baseline_timing.csv"))))
    env = json.load(open(os.path.join(folder, "environment.json")))
    ref = {m: {r["network"]: float(r["median_ns"]) for r in rows if r["method"] == m and r["median_ns"]} for m in PORT_REFS}
    L = [f"**{T.SMOKE_LABEL}**" if env.get("smoke") else "", "",
         f"Environment: {env.get('cpu_model')}, {env.get('cores_logical')} cores; git {env.get('git_matrix')}; date {env.get('date_local')}.", "",
         "| network | method | status | n | median (ns) | min (ns) | IQR (ns) | slow | quiet | median / port-full median | median / port-demand median |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        ratio = [f"{float(r['median_ns']) / ref[m][r['network']]:.17g}" if r["median_ns"] and r["network"] in ref[m] else ""
                 for m in PORT_REFS]
        L.append(f"| {r['network']} | {r['method']} | {r['status']} | {r['n']} | {r['median_ns']} | {r['min_ns']} | {r['iqr_ns']} | "
                 f"{r['slow']} | {r['quiet']} | {ratio[0]} | {ratio[1]} |")
    open(os.path.join(folder, "summary.md"), "w").write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("estimate"); e.add_argument("--cases", default=",".join(CASES)); e.add_argument("--budget", type=float, default=20.0)
    e.add_argument("--quiet-wait", type=float, default=300.0)
    r = sub.add_parser("run"); r.add_argument("--out"); r.add_argument("--cases", default=",".join(CASES))
    r.add_argument("--budget", type=float, default=20.0); r.add_argument("--min-rep", type=int, default=5)
    r.add_argument("--max-rep", type=int, default=200); r.add_argument("--warmup", type=int, default=2)
    r.add_argument("--timeout", type=float, default=600.0); r.add_argument("--slow", type=float, default=T.SLOW_S)
    r.add_argument("--quiet-wait", type=float, default=300.0); r.add_argument("--quiet-poll", type=float, default=10.0)
    r.add_argument("--smoke", action="store_true")
    r.add_argument("--guard-only", action="store_true", help="run the guards only and print them; no timing")
    a = ap.parse_args(argv)
    if a.cmd == "estimate":
        print(estimate(a.cases.split(","), a.budget, a.quiet_wait)); return 0
    if a.guard_only:
        print(guard_only(a.cases.split(","))); return 0
    if not a.out:
        ap.error("run needs --out (unless --guard-only)")
    run(a); return 0


if __name__ == "__main__":
    sys.exit(main())
