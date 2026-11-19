"""Scaling of the matrix method (Python port) with the number of sources (brief T14g, part A; T17 item 4).

Question: how does the time of the port change with the number k of sources, and what does one assembly and one
factorization with k right-hand sides save against k separate single-source computations?

Sources. A source is a column of the right-hand side B, which is also a column of the output table: in APFC the
generators with positive output and the buses with negative demand; in RPFC the generators, the injecting bus shunts,
the buses with negative demand and the terminal injections of the Case-3 branches. The self-term columns E_i_j_s of RPFC
are not columns of B (they come from the separate self-term solve) and are not counted as sources.

Subsets. For a network and kind with ns sources, the subsets are the first k entries of one fixed permutation of the
ns columns (numpy default_rng(SEED)), so they are deterministic and nested; k = 1, 2, 4, 8, ... below ns, and ns.

Modes (both through the port's public phase functions; the method code is not changed):
- "one run": one complete computation for k sources: input (Case.from_frames, read_branch_params), prepare, assemble,
  solve for the k selected columns of B (solve_active/solve_reactive on the system with B restricted to those columns
  and, in RPFC, without the self-term solve), and the table block of the k columns;
- "repeated": k such computations with one column each. It is measured up to k = 16 (KMAX_REPEATED); above, it is
  reported as k times the measured single-column time and labelled extrapolated.
The table block of the selected columns is formed with the formula of build_active_table / build_reactive_table:
demand rows, from-end rows X and to-end rows a X (and, in RPFC, the Case-3 rows with the own injection of a selected
Case-3 column), because the public builders always label every source. `check` verifies that the block equals the same
columns of the public table (apfc / rpfc) and that the k-column block equals the stack of the k single-column blocks.

Subcommands (from the matrix root):
    python -B benchmark/sources_scaling/sources_scaling.py check [--cases case9,case118] [--k 1,2,4]
    python -B benchmark/sources_scaling/sources_scaling.py estimate [--cases ...] [--budget 20]
    python -B benchmark/sources_scaling/sources_scaling.py run --out DIR [--cases ...] [--budget 20] [--smoke] [--kmax K]
    python -B benchmark/sources_scaling/sources_scaling.py plot DIR
"""
import argparse
import csv
import dataclasses
import json
import os
import sys
import time
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
import timing as T          # noqa: E402  (sets the thread counts before numpy is imported)
import cases as C           # noqa: E402
import numpy as np          # noqa: E402
from pfc_matrix import (Case, Options, read_branch_params, apfc, rpfc, prepare_active, assemble_active,  # noqa: E402
                        solve_active, prepare_reactive, assemble_reactive, solve_reactive)

SEED = 20261119
KMAX_REPEATED = 16
CASES = ["case118", "case_ACTIVSg500", "case_ACTIVSg2000"]
CSV_COLUMNS = ["label", "network", "kind", "k", "mode", "n", "median_ns", "min_ns", "iqr_ns", "quiet", "load1", "extrapolated"]


def k_values(ns, kmax=None):
    """1, 2, 4, ... below ns, and ns (each at most kmax if given)."""
    ks, k = [], 1
    while k < ns:
        ks.append(k); k *= 2
    ks.append(ns)
    return [k for k in ks if kmax is None or k <= kmax]


def subsets(ns, seed=SEED):
    """The fixed permutation of the ns source columns; subset k is its first k entries (nested)."""
    return np.random.default_rng(seed).permutation(ns)


def _block(kind, prep, system, X, idx):
    """Table rows of the source columns idx, as build_active_table / build_reactive_table compute them."""
    nl = prep.nl
    a = system.a if kind == "apfc" else system.aQ
    rows = [X[nl:], X[:nl], (a[:, None] * X)[:nl]]
    if kind == "rpfc" and prep.n_case3:
        hf, ht = np.asarray(prep.case3_hasf, bool), np.asarray(prep.case3_hast, bool)
        col_end = prep.n_pre + np.cumsum(hf.astype(int) + ht.astype(int))
        colf, colt = col_end - 1 - ht.astype(int), col_end - 1
        pos = {int(c): j for j, c in enumerate(idx)}
        c3f = np.zeros((int(hf.sum()), len(idx))); c3t = np.zeros((int(ht.sum()), len(idx)))
        for r, i in enumerate(np.where(hf)[0]):
            if int(colf[i]) in pos:
                c3f[r, pos[int(colf[i])]] = -prep.case3_QF[i]
        for r, i in enumerate(np.where(ht)[0]):
            if int(colt[i]) in pos:
                c3t[r, pos[int(colt[i])]] = prep.case3_QT[i]
        rows += [c3f, c3t]
    return np.vstack(rows)


def one_run(data, kind, idx):
    """One complete computation for the source columns idx (array of column indices of B). Returns (block, phases)."""
    t0 = time.perf_counter_ns()
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    t1 = time.perf_counter_ns()
    if kind == "apfc":
        prep = prepare_active(case, 1e-6); t2 = time.perf_counter_ns()
        system = assemble_active(prep, run_checks=False); t3 = time.perf_counter_ns()
        X = solve_active(dataclasses.replace(system, B=system.B[:, idx]))
    else:
        prep = prepare_reactive(case, 1e-6); t2 = time.perf_counter_ns()
        system = assemble_reactive(prep); t3 = time.perf_counter_ns()
        sub = dataclasses.replace(system, B=system.B[:, idx], B_self=system.B_self[:, :0])
        X, _ = solve_reactive(sub)
    t4 = time.perf_counter_ns()
    block = _block(kind, prep, system, X, idx)
    t5 = time.perf_counter_ns()
    return block, {"input": t1 - t0, "prepare": t2 - t1, "assemble": t3 - t2, "solve": t4 - t3, "table": t5 - t4}


def repeated(data, kind, idx):
    """k complete computations with one column each; returns the blocks side by side (and no phases)."""
    return np.hstack([one_run(data, kind, idx[j:j + 1])[0] for j in range(len(idx))]), None


def n_sources(data, kind):
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if kind == "apfc":
            return assemble_active(prepare_active(case, 1e-6), run_checks=False).B.shape[1]
        return assemble_reactive(prepare_reactive(case, 1e-6)).B.shape[1]


# ---------------------------------------------------------------- check

def check(cases, ks, tol=1e-12):
    """For each case, kind and k: the k-column block against the stack of the k single-column blocks, and against the
    same columns of the public table (apfc / rpfc, which compute every column). Returns the rows of the result."""
    out = []
    for c in cases:
        data = C.load(c, C.input_dirs())
        case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
        for kind in ("apfc", "rpfc"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                full, _ = (apfc if kind == "apfc" else rpfc)(case, Options(checks=False))
                ns = n_sources(data, kind)
                perm = subsets(ns)
                for k in ks:
                    if k > ns:
                        continue
                    idx = perm[:k]
                    blk, _ = one_run(data, kind, idx)
                    stack, _ = repeated(data, kind, idx)
                    pub = full.to_numpy()[:, idx]
                    d_stack = float(np.abs(blk - stack).max())
                    d_pub = float(np.abs(blk - pub).max()) if blk.shape == pub.shape else float("inf")
                    out.append({"case": c, "kind": kind, "ns": ns, "k": k, "rows": blk.shape[0],
                                "max_abs_k_vs_stack": d_stack, "max_abs_k_vs_public_table": d_pub,
                                "ok": d_stack <= tol and d_pub <= tol})
    return out


# ---------------------------------------------------------------- estimate and run

def plan(cases, kmax=None):
    """The cells of a run: (network, kind, k, mode) with mode 'one run' for every k and 'repeated' for k <= 16."""
    cells, ns_of = [], {}
    for c in cases:
        data = C.load(c, C.input_dirs())
        for kind in ("apfc", "rpfc"):
            ns = n_sources(data, kind); ns_of[(c, kind)] = ns
            for k in k_values(ns, kmax):
                cells.append((c, kind, k, "one run"))
                if k <= KMAX_REPEATED:
                    cells.append((c, kind, k, "repeated"))
    return cells, ns_of


def estimate(cases, budget, quiet_wait, kmax=None):
    cells, ns_of = plan(cases, kmax)
    lines = ["| quantity | value |", "|---|---|"]
    for (c, kind), ns in ns_of.items():
        lines.append(f"| sources of {c} {kind} | {ns} |")
    lines += [f"| measured cells | {len(cells)} |",
              f"| lower bound: cells x budget of {budget:g} s (s) | {len(cells) * budget:g} |",
              f"| added at most by the quiet wait, cells x {quiet_wait:g} s (s) | {len(cells) * quiet_wait:g} |"]
    lines.append("\nNot estimated (depends on run times): the warm-up runs, cells whose single run exceeds the budget "
                 "(at least 5 repetitions each, or 3 under the slow-cell rule), and the loading of the cases.")
    return "\n".join(lines)


def run(a):
    os.makedirs(a.out, exist_ok=True)
    env = T.environment_record(a.smoke, a)
    env["seed"] = SEED; env["kmax_repeated"] = KMAX_REPEATED
    ncpu = os.cpu_count()
    rows = []
    for c in a.cases.split(","):
        data = C.load(c, C.input_dirs())
        for kind in ("apfc", "rpfc"):
            ns = n_sources(data, kind)
            perm = subsets(ns)
            single = None
            for k in k_values(ns, a.kmax):
                idx = perm[:k]
                modes = [("one run", lambda idx=idx: one_run(data, kind, idx))]
                if k <= KMAX_REPEATED:
                    modes.append(("repeated", lambda idx=idx: repeated(data, kind, idx)))
                for mode, fn in modes:
                    q = T.wait_quiet(ncpu, a.quiet_wait, a.quiet_poll)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        r = T.measure_cell(fn, a.budget, a.min_rep, a.max_rep, a.warmup, a.timeout, slow_s=a.slow)
                    s = T.stats(r["samples"])
                    row = {"label": env["label"], "network": c, "kind": kind, "k": k, "mode": mode, "n": s["n"],
                           "median_ns": s["median"], "min_ns": s["min"], "iqr_ns": s["iqr"], "quiet": q["quiet"],
                           "load1": q["load1_end"], "extrapolated": False}
                    rows.append(row)
                    if mode == "repeated" and k == 1:
                        single = row
                if k > KMAX_REPEATED and single is not None:
                    rows.append({"label": env["label"], "network": c, "kind": kind, "k": k, "mode": "repeated",
                                 "n": single["n"], "median_ns": k * single["median_ns"], "min_ns": k * single["min_ns"],
                                 "iqr_ns": k * single["iqr_ns"], "quiet": single["quiet"], "load1": single["load1"],
                                 "extrapolated": True})
    with open(os.path.join(a.out, "sources.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(a.out, "environment.json"), "w") as fh:
        json.dump(env, fh, indent=1, default=str)
    summarize(a.out)
    return rows


def summarize(folder):
    rows = list(csv.DictReader(open(os.path.join(folder, "sources.csv"))))
    env = json.load(open(os.path.join(folder, "environment.json")))
    smoke = bool(env.get("smoke"))
    md = [f"**{T.SMOKE_LABEL}**" if smoke else "", "",
          f"Environment: {env.get('cpu_model')}, {env.get('cores_logical')} cores; git {env.get('git_matrix')}; date {env.get('date_local')}.", "",
          "| network | kind | k | mode | n | median (ns) | min (ns) | IQR (ns) | quiet | extrapolated | repeated / one run (medians) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    one = {(r["network"], r["kind"], r["k"]): float(r["median_ns"]) for r in rows if r["mode"] == "one run" and r["median_ns"]}
    for r in rows:
        ratio = ""
        if r["mode"] == "repeated" and r["median_ns"] and (r["network"], r["kind"], r["k"]) in one:
            ratio = f"{float(r['median_ns']) / one[(r['network'], r['kind'], r['k'])]:.17g}"
        md.append(f"| {r['network']} | {r['kind']} | {r['k']} | {r['mode']} | {r['n']} | {r['median_ns']} | {r['min_ns']} | "
                  f"{r['iqr_ns']} | {r['quiet']} | {r['extrapolated']} | {ratio} |")
    open(os.path.join(folder, "summary.md"), "w").write("\n".join(md) + "\n")
    plot(folder)


def plot(folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader(open(os.path.join(folder, "sources.csv"))))
    env = json.load(open(os.path.join(folder, "environment.json")))
    suffix = f" ({T.SMOKE_LABEL})" if env.get("smoke") else ""
    keys = sorted({(r["network"], r["kind"]) for r in rows})
    fig, ax = plt.subplots(figsize=(7, 4.8)); fig2, ax2 = plt.subplots(figsize=(7, 4.8))
    for net, kind in keys:
        for mode, ls in (("one run", "-"), ("repeated", "--")):
            pts = sorted((int(r["k"]), float(r["median_ns"]) / 1e9, r["extrapolated"] == "True") for r in rows
                         if r["network"] == net and r["kind"] == kind and r["mode"] == mode and r["median_ns"])
            if not pts:
                continue
            line, = ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, marker="o", ms=3, label=f"{net} {kind.upper()}, {mode}")
            ext = [p for p in pts if p[2]]
            if ext:
                ax.plot([p[0] for p in ext], [p[1] for p in ext], "x", color=line.get_color(), ms=7,
                        label=f"{net} {kind.upper()}, repeated, extrapolated (k x single-column time)")
        one = {int(r["k"]): float(r["median_ns"]) for r in rows if r["network"] == net and r["kind"] == kind and r["mode"] == "one run" and r["median_ns"]}
        rep = sorted((int(r["k"]), float(r["median_ns"]), r["extrapolated"] == "True") for r in rows
                     if r["network"] == net and r["kind"] == kind and r["mode"] == "repeated" and r["median_ns"] and int(r["k"]) in one)
        if rep:
            line, = ax2.plot([p[0] for p in rep], [p[1] / one[p[0]] for p in rep], marker="o", ms=3, label=f"{net} {kind.upper()}")
            ext = [p for p in rep if p[2]]
            if ext:
                ax2.plot([p[0] for p in ext], [p[1] / one[p[0]] for p in ext], "x", color=line.get_color(), ms=7)
    for a_, yl, title, name, f in ((ax, "median time (s)", "time against the number of sources", "sources_time", fig),
                                   (ax2, "repeated / one run (medians)", "ratio; x: extrapolated", "sources_ratio", fig2)):
        a_.set_xscale("log", base=2); a_.set_yscale("log"); a_.set_xlabel("number of sources k")
        a_.set_ylabel(yl); a_.set_title(title + suffix, fontsize=9); a_.grid(True, which="both", lw=0.3)
        if a_.get_legend_handles_labels()[0]:
            a_.legend(fontsize=6)
        f.tight_layout()
        for ext in ("pdf", "png"):
            f.savefig(os.path.join(folder, f"{name}.{ext}"), dpi=200)
        plt.close(f)


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check"); c.add_argument("--cases", default="case9,case118"); c.add_argument("--k", default="1,2,4")
    e = sub.add_parser("estimate"); e.add_argument("--cases", default=",".join(CASES)); e.add_argument("--budget", type=float, default=20.0)
    e.add_argument("--quiet-wait", type=float, default=300.0); e.add_argument("--kmax", type=int)
    r = sub.add_parser("run"); r.add_argument("--out", required=True); r.add_argument("--cases", default=",".join(CASES))
    r.add_argument("--budget", type=float, default=20.0); r.add_argument("--min-rep", type=int, default=5)
    r.add_argument("--max-rep", type=int, default=200); r.add_argument("--warmup", type=int, default=2)
    r.add_argument("--timeout", type=float, default=600.0); r.add_argument("--slow", type=float, default=T.SLOW_S)
    r.add_argument("--quiet-wait", type=float, default=300.0); r.add_argument("--quiet-poll", type=float, default=10.0)
    r.add_argument("--kmax", type=int, help="largest k (for a smoke run)"); r.add_argument("--smoke", action="store_true")
    p = sub.add_parser("plot"); p.add_argument("folder")
    a = ap.parse_args(argv)
    if a.cmd == "check":
        res = check(a.cases.split(","), [int(x) for x in a.k.split(",")])
        print("| case | kind | sources | k | rows | max abs, k columns against the stack of k single columns | max abs, against the public table | ok |")
        print("|---|---|---|---|---|---|---|---|")
        for x in res:
            print(f"| {x['case']} | {x['kind']} | {x['ns']} | {x['k']} | {x['rows']} | {x['max_abs_k_vs_stack']:.17g} | {x['max_abs_k_vs_public_table']:.17g} | {x['ok']} |")
        bad = sum(not x["ok"] for x in res)
        print(f"\nfailing: {bad}")
        return 1 if bad else 0
    if a.cmd == "estimate":
        print(estimate(a.cases.split(","), a.budget, a.quiet_wait, a.kmax)); return 0
    if a.cmd == "run":
        run(a); return 0
    plot(a.folder); return 0


if __name__ == "__main__":
    sys.exit(main())
