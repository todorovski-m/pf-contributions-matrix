"""Runtime of the Achayuthakan et al. (2010) baseline against the port (brief T14h, step 4).

For each network of the comparison (case1, case9, case14, case30, case118, case_ACTIVSg500, case_ACTIVSg2000), the
complete computation from the power-flow quantities (the DataFrames of the harness's loader) to the demand-level
contribution table is timed for
- "port": the port with checks off (Case.from_frames and read_branch_params, apfc, then the demand rows with the
  generator columns summed per bus, as compare.apfc_bus_table);
- "net-inverse", "net-sparse", "gross-inverse", "gross-sparse": the baseline (achayuthakan.contributions), both loss
  variants, with the explicit inverse (variant (a), as published) and with one sparse LU (variant (b)).
The machinery of benchmark/timing/ is imported and not changed: measure_cell (warm-up, repetitions within a budget,
median, time-out, slow-cell rule), wait_quiet before each cell (outside the timed region), stats and the environment
record. A guard runs first for each network: the net variant must equal APFC at demand level within 1e-9 MW and the
two solution methods within 1e-10; a network that fails it is not timed.

Usage (from the matrix root):
    python -B benchmark/baseline_achayuthakan/baseline_timing.py estimate [--budget 20]
    python -B benchmark/baseline_achayuthakan/baseline_timing.py run --out DIR [--cases ...] [--budget 20] [--smoke]
Outputs: baseline_timing.csv (label, network, method, status, n, median_ns, min_ns, iqr_ns, slow, quiet, load1),
guard.csv, environment.json, summary.md (with the ratio of each baseline method to the port).
"""
import argparse
import csv
import json
import os
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
import timing as T          # noqa: E402  (thread settings before numpy)
import numpy as np          # noqa: E402
import compare as CMP       # noqa: E402
import achayuthakan as A    # noqa: E402

CASES = list(CMP.CASES)
METHODS = ["port", "net-inverse", "net-sparse", "gross-inverse", "gross-sparse"]
CSV_COLUMNS = ["label", "network", "method", "status", "n", "median_ns", "min_ns", "iqr_ns", "slow", "quiet", "load1"]


def runner(data, method):
    """The timed function of one cell: the complete computation, returning (table, None)."""
    if method == "port":
        return lambda: (CMP.apfc_bus_table(data)[0], None)
    kind, how = method.split("-")
    return lambda: (A.contributions(data["bus"], data["gen"], data["branch"], kind=kind, method=how), None)


def guard(data):
    ap = CMP.apfc_bus_table(data)[0]
    net = A.contributions(data["bus"], data["gen"], data["branch"], kind="net", method="sparse")
    d_apfc = float(np.abs(net.reindex(index=ap.index, columns=ap.columns, fill_value=0.0).to_numpy() - ap.to_numpy()).max())
    d_meth = max(float(np.abs(A.contributions(data["bus"], data["gen"], data["branch"], kind=k, method="inverse").to_numpy()
                              - A.contributions(data["bus"], data["gen"], data["branch"], kind=k, method="sparse").to_numpy()).max())
                 for k in ("net", "gross"))
    return {"net_vs_apfc": d_apfc, "inverse_vs_sparse": d_meth, "ok": d_apfc <= 1e-9 and d_meth <= 1e-10}


def estimate(cases, budget, quiet_wait):
    cells = len(cases) * len(METHODS)
    return "\n".join(["| quantity | value |", "|---|---|", f"| networks | {len(cases)} |", f"| cells (networks x 5 methods) | {cells} |",
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
        w = csv.DictWriter(fh, fieldnames=["network", "net_vs_apfc", "inverse_vs_sparse", "ok"], lineterminator="\n"); w.writeheader()
        for r in guards:
            w.writerow(r)
    with open(os.path.join(a.out, "environment.json"), "w") as fh:
        json.dump(env, fh, indent=1, default=str)
    summarize(a.out)
    return rows


def summarize(folder):
    rows = list(csv.DictReader(open(os.path.join(folder, "baseline_timing.csv"))))
    env = json.load(open(os.path.join(folder, "environment.json")))
    port = {r["network"]: float(r["median_ns"]) for r in rows if r["method"] == "port" and r["median_ns"]}
    L = [f"**{T.SMOKE_LABEL}**" if env.get("smoke") else "", "",
         f"Environment: {env.get('cpu_model')}, {env.get('cores_logical')} cores; git {env.get('git_matrix')}; date {env.get('date_local')}.", "",
         "| network | method | status | n | median (ns) | min (ns) | IQR (ns) | slow | quiet | median / port median |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        ratio = f"{float(r['median_ns']) / port[r['network']]:.17g}" if r["median_ns"] and r["network"] in port else ""
        L.append(f"| {r['network']} | {r['method']} | {r['status']} | {r['n']} | {r['median_ns']} | {r['min_ns']} | {r['iqr_ns']} | "
                 f"{r['slow']} | {r['quiet']} | {ratio} |")
    open(os.path.join(folder, "summary.md"), "w").write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("estimate"); e.add_argument("--cases", default=",".join(CASES)); e.add_argument("--budget", type=float, default=20.0)
    e.add_argument("--quiet-wait", type=float, default=300.0)
    r = sub.add_parser("run"); r.add_argument("--out", required=True); r.add_argument("--cases", default=",".join(CASES))
    r.add_argument("--budget", type=float, default=20.0); r.add_argument("--min-rep", type=int, default=5)
    r.add_argument("--max-rep", type=int, default=200); r.add_argument("--warmup", type=int, default=2)
    r.add_argument("--timeout", type=float, default=600.0); r.add_argument("--slow", type=float, default=T.SLOW_S)
    r.add_argument("--quiet-wait", type=float, default=300.0); r.add_argument("--quiet-poll", type=float, default=10.0)
    r.add_argument("--smoke", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "estimate":
        print(estimate(a.cases.split(","), a.budget, a.quiet_wait)); return 0
    run(a); return 0


if __name__ == "__main__":
    sys.exit(main())
