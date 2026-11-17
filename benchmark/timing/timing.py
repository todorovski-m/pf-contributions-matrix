"""Timing harness: the matrix method (Python port, pfc_matrix) against the old iterative method (Python), both from
the same in-memory DataFrames of load_results. See README.md for the timed regions and the statistics.

Usage (from matrix/, with the global Python):
    python -B benchmark/timing/timing.py run --out DIR [--cases c1,c2,...] [--budget 30] [--min-rep 5]
        [--max-rep 200] [--warmup 2] [--timeout 600] [--itermax 10000] [--smoke] [--no-memory]
        [--syn-inputs DIR] [--large-inputs DIR]
    python -B benchmark/timing/timing.py memory --case C --method port|old --kind apfc|rpfc --checks 0|1 ...
The thread counts of the numerical libraries are set to 1 before numpy is imported.
"""
import os

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS")
for _v in THREAD_VARS:
    os.environ[_v] = "1"

import argparse          # noqa: E402
import csv               # noqa: E402
import datetime          # noqa: E402
import gc                # noqa: E402
import io                # noqa: E402
import json              # noqa: E402
import logging           # noqa: E402
import platform          # noqa: E402
import re                # noqa: E402
import resource          # noqa: E402
import signal            # noqa: E402
import subprocess        # noqa: E402
import sys               # noqa: E402
import time              # noqa: E402
import warnings          # noqa: E402
from contextlib import redirect_stdout  # noqa: E402

import numpy as np       # noqa: E402

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases as C        # noqa: E402

SMOKE_LABEL = "smoke, not for publication"
TOL, TOL_DIV, GUARD_TOL = 1e-6, 1e-9, 1e-10
PORT_PHASES = ("input", "prepare", "assemble", "solve", "table", "condest")


# ---------------------------------------------------------------- statistics and measurement

def stats(samples):
    """n, min, median, mean, standard deviation (ddof = 1; 0 for one sample) and interquartile range (75th minus
    25th percentile, linear interpolation) of a list of nanosecond samples."""
    x = np.asarray(samples, dtype=float)
    if x.size == 0:
        return {"n": 0, "min": None, "median": None, "mean": None, "std": None, "iqr": None}
    q75, q25 = np.percentile(x, [75, 25])
    return {"n": int(x.size), "min": float(x.min()), "median": float(np.median(x)), "mean": float(x.mean()),
            "std": float(x.std(ddof=1)) if x.size > 1 else 0.0, "iqr": float(q75 - q25)}


class RunTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise RunTimeout()


def call_with_timeout(fn, timeout_s):
    """fn() with a time-out of timeout_s seconds (SIGALRM; interrupts Python code and sleeps). Raises RunTimeout."""
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, timeout_s)
    try:
        return fn()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


SLOW_S, SLOW_REP = 30.0, 3     # the slow-cell rule


def measure_cell(fn, budget_s=30.0, min_rep=5, max_rep=200, warmup=2, timeout_s=600.0, finished=None, slow_s=SLOW_S):
    """Measure one cell. fn() returns (result, phases), phases a dict of nanoseconds or None. gc.collect() before the
    cell, the garbage collector disabled during each single run, time.perf_counter_ns around fn. ``warmup`` runs are
    not recorded. Then repetitions until their total time reaches budget_s with at least min_rep, at most max_rep.
    A single run longer than timeout_s ends the cell with status 'did not finish (time-out)'; if finished(result) is
    False, the status is 'did not finish (iteration cap)'.
    Slow-cell rule: if a single run (warm-up or recorded) takes longer than slow_s seconds, the cell is 'slow': no
    further warm-up run, and exactly SLOW_REP = 3 repetitions are recorded (fewer only on a time-out); min_rep and
    the budget do not apply. Returns dict(status, samples, phases, warmup, slow)."""
    gc.collect()
    samples, phases = [], []

    def one():
        gc.disable()
        try:
            t0 = time.perf_counter_ns()
            res, ph = call_with_timeout(fn, timeout_s)
            t1 = time.perf_counter_ns()
        finally:
            gc.enable()
        return t1 - t0, res, ph

    slow, done_warm = False, 0
    out = lambda st: {"status": st, "samples": samples, "phases": phases, "warmup": done_warm, "slow": slow}
    try:
        for _ in range(warmup):
            ns, res, _ = one(); done_warm += 1
            if finished is not None and not finished(res):
                return out("did not finish (iteration cap)")
            if ns > slow_s * 1e9:
                slow = True
                break
        total = 0
        while True:
            if slow:
                if len(samples) >= SLOW_REP:
                    break
            elif not (len(samples) < max_rep and (len(samples) < min_rep or total < budget_s * 1e9)):
                break
            ns, res, ph = one()
            if finished is not None and not finished(res):
                return out("did not finish (iteration cap)")
            samples.append(ns); phases.append(ph); total += ns
            if ns > slow_s * 1e9 and not slow:
                slow = True
    except RunTimeout:
        return out("did not finish (time-out)")
    return out("ok")


# ---------------------------------------------------------------- the two methods

def port_runner(data, kind, checks):
    """The timed function of the port: Case.from_frames and read_branch_params (included, conservative against the
    port), then the public phase functions (a) prepare, (b) assemble, (c) solve, (d) table, and, with checks, the
    condition estimate (as apfc/rpfc with Options(checks=True)). No CSV writing. Returns (table, phases in ns)."""
    from pfc_matrix import (Case, read_branch_params, prepare_active, assemble_active, solve_active,
                            build_active_table, prepare_reactive, assemble_reactive, solve_reactive,
                            build_reactive_table)
    from pfc_matrix.checks import condest
    bus, gen, branch, path, name = data["bus"], data["gen"], data["branch"], data["path"], data["case"]

    def run():
        pc = time.perf_counter_ns
        t0 = pc()
        case = Case.from_frames(bus, gen, branch, name, params=read_branch_params(path))
        t1 = pc()
        if kind == "apfc":
            prep = prepare_active(case, TOL); t2 = pc()
            system = assemble_active(prep, TOL_DIV, checks); t3 = pc()
            X = solve_active(system); t4 = pc()
            table = build_active_table(prep, X, system.a); t5 = pc()
        else:
            prep = prepare_reactive(case, TOL); t2 = pc()
            system = assemble_reactive(prep, TOL, TOL_DIV); t3 = pc()
            X, X_self = solve_reactive(system); t4 = pc()
            table = build_reactive_table(prep, system, X, X_self); t5 = pc()
        t6 = t5
        if checks:
            condest(system.A); t6 = pc()
        ph = {"input": t1 - t0, "prepare": t2 - t1, "assemble": t3 - t2, "solve": t4 - t3, "table": t5 - t4}
        if checks:
            ph["condest"] = t6 - t5
        return table, ph
    return run


def old_runner(data, kind, itermax, branch=None):
    """The timed function of the old method: makeAPFC or makeRPFC on the DataFrames (bus_pos None, graph False). It
    copies the branch table itself, so repeated runs see the same input. Returns ((flows, table), None)."""
    _, pc, _ = C.import_old()
    bus, gen = data["bus"], data["gen"]
    br = data["branch"] if branch is None else branch
    f = pc.makeAPFC if kind == "apfc" else pc.makeRPFC

    def run():
        return f(bus, gen, br, None, itermax, False), None
    return run


def old_finished(res):
    (_, table) = res
    return bool(table["DONE"].all())


class _Collect(logging.Handler):
    def __init__(self):
        super().__init__(logging.INFO); self.lines = []

    def emit(self, r):
        self.lines.append(r.getMessage())


def old_with_iterations(data, kind, itermax, branch, timeout_s):
    """One untimed run of the old method with an INFO collector (no output), for the iteration count and the guard.
    Returns (flows, table, done, iterations) or raises RunTimeout."""
    _, pc, _ = C.import_old()
    root = logging.getLogger(); h = _Collect(); root.addHandler(h); level = root.level; root.setLevel(logging.INFO)
    try:
        (flows, table), _ = call_with_timeout(old_runner(data, kind, itermax, branch), timeout_s)
    finally:
        root.removeHandler(h); root.setLevel(level)
    name = kind.upper()
    m = [re.search(rf"{name}: nr = \d+, iters = (\d+), DONE = (\w+)", l) for l in h.lines]
    m = [x for x in m if x][-1]
    return flows, table, m.group(2) == "True", int(m.group(1))


# ---------------------------------------------------------------- memory estimate of the dense result tables

def dense_bytes(prep, kind):
    """Estimate of the memory of the dense arrays of one port run, 8 bytes per value: apfc the solution X (ne x ng)
    and the table (nd + 2 nl rows, ng columns); rpfc the from-end and to-end blocks (ne x (ns + nl), twice) and the
    table (nd + 2 nl + 2 n_case3 rows, ns + nl columns). The old method's dense array is of the same order."""
    if kind == "apfc":
        return 8 * (prep.ne * prep.ng + (prep.nd + 2 * prep.nl) * prep.ng)
    cols = prep.ns + prep.nl
    return 8 * (2 * prep.ne * cols + (prep.nd + 2 * prep.nl + 2 * prep.n_case3) * cols)


def estimate_case(data, kind):
    from pfc_matrix import Case, read_branch_params, prepare_active, prepare_reactive
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    prep = (prepare_active if kind == "apfc" else prepare_reactive)(case, TOL)
    return dense_bytes(prep, kind)


def physical_memory():
    try:
        return int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout.strip())
    except Exception:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")


# ---------------------------------------------------------------- guard

def port_flows(data, kind):
    """The table of the phase composition, the table of the public function (must be equal), and the flows of the rows."""
    from pfc_matrix import apfc, rpfc, Options, Case, read_branch_params, prepare_active, prepare_reactive
    table, _ = port_runner(data, kind, False)()
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pub, info = (apfc if kind == "apfc" else rpfc)(case, Options(checks=False))
    if kind == "apfc":
        p = prepare_active(case, TOL); nl = p.nl
        flows = np.r_[p.Pf[nl:], p.Pf[:nl], p.Pt[:nl]]
    else:
        p = prepare_reactive(case, TOL); nl = p.nl
        flows = np.r_[p.Qf[nl:], p.Qf[:nl], p.Qt[:nl], -p.case3_QF, p.case3_QT]
    return table, pub, flows, info


def compare_old_in_memory(table_old, fq, branch_raw, port_table, kind, done):
    """The comparison of tests/compare_old.py (its helper functions, unchanged), on in-memory tables: rows mapped, the
    largest cell difference (only if the old method finished) and the largest row-sum deviation of the old method."""
    _, _, co = C.import_old()
    M = port_table.reset_index()
    mk = co.matrix_keys(M["Line"], kind); mpos = {k: i for i, k in enumerate(mk)}
    mcols = list(M.columns[1:])
    src_end = list(table_old.columns).index("F_BUS")
    cm, unm_o, unm_m = co.colmap(list(table_old.columns[:src_end]), mcols, kind, fq)
    branch_from = [int(x) for x in branch_raw.F_BUS]
    okeys = [co.old_row_key(i, r, kind, branch_from) for i, r in table_old.iterrows()]
    matched, mx, rsmax = set(), 0.0, 0.0
    Mv = M.iloc[:, 1:].to_numpy(float); col = {c: k for k, c in enumerate(mcols)}
    for j, ((k, n), idx) in enumerate(zip(co.occ(okeys), table_old.index)):
        if k in (None, "skip") or (k, n) not in mpos:
            continue
        i = mpos[(k, n)]; matched.add(i)
        orow = table_old.iloc[j]
        flow = float(orow["Total"]) - float(orow["Error"])
        rsmax = max(rsmax, abs(float(Mv[i].sum()) - flow))
        if done:
            for oc, mc in cm.items():
                mx = max(mx, abs(float(orow[oc]) - float(sum(Mv[i, col[c]] for c in mc))))
    unmapped = [M["Line"][i] for i in range(len(M)) if i not in matched]
    return {"rows_mapped": len(matched), "rows": len(M), "unmapped_matrix_rows": unmapped, "max_cell": mx,
            "old_rowsum": rsmax, "unmapped_old_columns": unm_o, "unmapped_matrix_columns": unm_m}


def guard(data, kind, itermax, timeout_s, branch_old):
    """Before timing: the port (phase composition equal to the public function; row sums equal the flows within 1e-10)
    and the old method in memory (with the corrected QSF; cells within 1e-10 if it finishes; row sums within 1e-10).
    Returns a dict with 'ok', 'old_status', 'iterations' and the numbers."""
    out = {"case": data["case"], "kind": kind}
    table, pub, flows, info = port_flows(data, kind)
    out["n_e"] = info["n_e"]; out["nnz_N"] = info["nnz_N"]
    out["port_equal_public"] = bool(table.equals(pub))
    out["port_rowsum"] = float(np.abs(table.to_numpy().sum(axis=1) - flows).max()) if len(flows) else 0.0
    try:
        import correctness           # the row mapping with the terminal classification (T15b, T14c)
        flows_old, t_old, done, iters = old_with_iterations(data, kind, itermax, branch_old, timeout_s)
        out["old_status"] = "ok" if done else "did not finish (iteration cap)"
        out["iterations"] = iters
        c = correctness.compare(data, kind, t_old, flows_old, table)
        cmpd = {"rows": len(table), "unmapped_matrix_rows": c["unmapped_port_rows"], "max_cell": c["max_cell"],
                "old_rowsum": c["port_rowsum_vs_old_flow"], "excluded_rows": c["excluded_rows"],
                "unmapped_old_columns": c["unmapped_old_columns"], "unmapped_matrix_columns": c["unmapped_port_columns"]}
        out.update(cmpd)
    except RunTimeout:
        out["old_status"] = "did not finish (time-out)"; out["iterations"] = None
        cmpd = None
    ok = out["port_equal_public"] and out["port_rowsum"] <= GUARD_TOL
    if cmpd is not None:
        ok = ok and not cmpd["unmapped_matrix_rows"] and cmpd["old_rowsum"] <= GUARD_TOL
        if out["old_status"] == "ok":
            ok = ok and cmpd["max_cell"] <= GUARD_TOL
    out["ok"] = bool(ok)
    return out


# ---------------------------------------------------------------- environment

def _cmd(args):
    try:
        return subprocess.run(args, capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception as e:      # the command is missing on this platform
        return f"unavailable ({type(e).__name__})"


def blas_config():
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            cfg = np.show_config(mode="dicts")
        return cfg if cfg is not None else buf.getvalue()
    except TypeError:
        with redirect_stdout(buf):
            np.show_config()
        return buf.getvalue()


def load_record():
    l1, l5, l15 = os.getloadavg()
    return {"load1": l1, "load5": l5, "load15": l15}


def quiet_flag(load1_before, ncpu):
    """quiet: the 1-minute load average before the case is below 0.25 times the number of cores."""
    return bool(load1_before < 0.25 * ncpu)


def environment_record(smoke=False, args=None):
    import scipy, pandas, networkx
    ncpu = os.cpu_count()
    rec = {
        "label": SMOKE_LABEL if smoke else "series",
        "smoke": bool(smoke),
        "date_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "date_local": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "cpu_model": _cmd(["sysctl", "-n", "machdep.cpu.brand_string"]) or platform.processor(),
        "cores_logical": ncpu,
        "cores_physical": _cmd(["sysctl", "-n", "hw.physicalcpu"]),
        "memory_bytes": _cmd(["sysctl", "-n", "hw.memsize"]),
        "os": f"{platform.system()} {platform.release()}; " + _cmd(["sw_vers"]).replace("\n", "; "),
        "python": sys.version, "python_executable": sys.executable,
        "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pandas.__version__, "networkx": networkx.__version__,
        "threads": {v: os.environ.get(v) for v in THREAD_VARS},
        "blas": blas_config(),
        "git_matrix": _cmd(["git", "-C", C.MATRIX, "rev-parse", "HEAD"]),
        "git_matrix_dirty": bool(_cmd(["git", "-C", C.MATRIX, "status", "--porcelain"])),
        "git_submodule": _cmd(["git", "-C", C.SUB, "rev-parse", "HEAD"]),
        "power": _cmd(["pmset", "-g", "batt"]),
        "caffeinate_env": os.environ.get("PFC_CAFFEINATE", "0") == "1",
        "caffeinate_process": bool(_cmd(["pgrep", "-x", "caffeinate"])),
        "ru_maxrss_unit": "bytes" if sys.platform == "darwin" else "kibibytes",
        "args": vars(args) if args is not None else None,
    }
    return rec


# ---------------------------------------------------------------- memory (fresh subprocess per method and case)

def memory_child(a):
    """One run of one method in this (fresh) process; prints JSON with ru_maxrss before and after the run."""
    data = C.load(a.case, C.input_dirs(a.syn_inputs, a.large_inputs))
    if a.method == "port":
        fn = port_runner(data, a.kind, bool(a.checks))
    else:
        branch, _, _, _ = C.correct_qsf(data)
        fn = old_runner(data, a.kind, a.itermax, branch)
    gc.collect()
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fn()
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({"ru_maxrss_before": before, "ru_maxrss_after": after}))


def memory_cell(a, case, method, kind, checks):
    cmd = [sys.executable, "-B", os.path.abspath(__file__), "memory", "--case", case, "--method", method, "--kind", kind,
           "--checks", str(int(checks)), "--itermax", str(a.itermax)]
    if a.syn_inputs:
        cmd += ["--syn-inputs", a.syn_inputs]
    if a.large_inputs:
        cmd += ["--large-inputs", a.large_inputs]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout + 60)
        if p.returncode != 0:
            return {"status": f"error (exit {p.returncode})", "ru_maxrss_before": None, "ru_maxrss_after": None}
        d = json.loads(p.stdout.strip().splitlines()[-1]); d["status"] = "ok"
        return d
    except subprocess.TimeoutExpired:
        return {"status": "did not finish (time-out)", "ru_maxrss_before": None, "ru_maxrss_after": None}


# ---------------------------------------------------------------- the series

def run_series(a):
    os.makedirs(a.out, exist_ok=True)
    env = environment_record(a.smoke, a)
    dirs = C.input_dirs(a.syn_inputs, a.large_inputs)
    combos = C.parse_combos(a.combos) if a.combos else [(c, ("apfc", "rpfc")) for c in (a.cases.split(",") if a.cases else C.SMALL)]
    ncpu = os.cpu_count()
    samples_f = open(os.path.join(a.out, "samples.csv"), "w", newline="")
    sw = csv.writer(samples_f)
    sw.writerow(["label", "case", "method", "kind", "checks", "repetition", "ns", "phase"])
    cells, guards, iters, mem, per_case = [], [], [], [], {}
    for c, kinds in combos:
        before = load_record()
        rec = {"load_before": before, "quiet": quiet_flag(before["load1"], ncpu)}
        try:
            data = C.load(c, dirs)
            branch_old, rows, q_old, q_new = C.correct_qsf(data)
            rec["qsf_replaced_rows"] = [int(r) for r in rows]
            rec["qsf_replaced"] = [{"row": int(r), "from": int(data["branch"].F_BUS.iloc[r]), "to": int(data["branch"].T_BUS.iloc[r]),
                                    "exported": a_, "rpfc_m": b_} for r, a_, b_ in zip(rows, q_old, q_new)]
        except Exception as e:        # the case cannot be loaded or the port rejects the input
            rec["error"] = f"{type(e).__name__}: {getattr(e, 'identifier', '')} {str(e)[:300]}"
            rec["load_after"] = load_record(); per_case[c] = rec
            cells.append({"case": c, "method": "all", "kind": "", "checks": "", "status": f"not timed ({rec['error']})"})
            continue
        limit = a.mem_fraction * physical_memory()
        for kind in kinds:
            try:
                est = estimate_case(data, kind)
            except Exception as e:
                cells.append({"case": c, "method": "all", "kind": kind, "checks": "",
                              "status": f"not timed (port raises {getattr(e, 'identifier', type(e).__name__)})"})
                continue
            rec[f"dense_bytes_{kind}"] = est
            if est > limit:
                cells.append({"case": c, "method": "all", "kind": kind, "checks": "",
                              "status": f"not timed (estimated dense tables {est} bytes > limit {int(limit)} bytes)"})
                continue
            try:
                g = guard(data, kind, a.itermax, a.timeout, branch_old)
            except Exception as e:
                cells.append({"case": c, "method": "all", "kind": kind, "checks": "",
                              "status": f"not timed (port raises {getattr(e, 'identifier', type(e).__name__)}: {str(e)[:200]})"})
                continue
            guards.append(g)
            iters.append({"case": c, "kind": kind, "old_status": g["old_status"], "iterations": g["iterations"], "itermax": a.itermax})
            rec[f"n_e_{kind}"] = g["n_e"]
            if not g["ok"]:
                cells.append({"case": c, "method": "all", "kind": kind, "checks": "", "status": "not timed (guard failed)"})
                continue
            series = [("port", False), ("port", True)]
            if g["old_status"] == "ok":
                series.append(("old", None))
            else:
                cells.append({"case": c, "method": "old", "kind": kind, "checks": "", "status": g["old_status"] + " (guard run)", "n": 0})
            for method, checks in series:
                fn = port_runner(data, kind, checks) if method == "port" else old_runner(data, kind, a.itermax, branch_old)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r = measure_cell(fn, a.budget, a.min_rep, a.max_rep, a.warmup, a.timeout,
                                     finished=old_finished if method == "old" else None, slow_s=a.slow)
                ch = "" if checks is None else ("on" if checks else "off")
                for k, (ns, ph) in enumerate(zip(r["samples"], r["phases"]), 1):
                    sw.writerow([env["label"], c, method, kind, ch, k, ns, ""])
                    for p_, v in (ph or {}).items():
                        sw.writerow([env["label"], c, method, kind, ch, k, v, p_])
                samples_f.flush()
                cells.append({"case": c, "method": method, "kind": kind, "checks": ch, "status": r["status"],
                              "n": len(r["samples"]), "warmup": r["warmup"], "slow": r["slow"]})
                if a.memory and r["status"] == "ok":
                    m = memory_cell(a, c, method, kind, bool(checks))
                    mem.append({"case": c, "method": method, "kind": kind, "checks": ch, **m, "unit": env["ru_maxrss_unit"]})
        rec["load_after"] = load_record()
        per_case[c] = rec
    samples_f.close()
    env["cases"] = per_case
    env["qsf_replacement_cases"] = [c for c, r in per_case.items() if r.get("qsf_replaced_rows")]
    with open(os.path.join(a.out, "environment.json"), "w") as fh:
        json.dump(env, fh, indent=1, default=str)
    for name, rows in (("cells", cells), ("guard", guards), ("iterations", iters), ("memory", mem)):
        keys = sorted({k for r in rows for k in r}) if rows else ["case"]
        with open(os.path.join(a.out, f"{name}.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["label"] + keys)
            w.writeheader()
            for r in rows:
                w.writerow({"label": env["label"], **{k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()}})
    return env, cells, guards


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "memory"):
        p = sub.add_parser(name)
        p.add_argument("--syn-inputs"); p.add_argument("--large-inputs")
        p.add_argument("--itermax", type=int, default=10000)
        p.add_argument("--timeout", type=float, default=600.0)
    r = sub.choices["run"]
    r.add_argument("--out", required=True); r.add_argument("--cases")
    r.add_argument("--combos", help="case:kind list, for example case9:apfc,case9:rpfc,case3120sp:rpfc (wins over --cases)")
    r.add_argument("--slow", type=float, default=SLOW_S, help="slow-cell threshold in seconds for a single run")
    r.add_argument("--budget", type=float, default=30.0); r.add_argument("--min-rep", type=int, default=5)
    r.add_argument("--max-rep", type=int, default=200); r.add_argument("--warmup", type=int, default=2)
    r.add_argument("--smoke", action="store_true"); r.add_argument("--no-memory", dest="memory", action="store_false")
    r.add_argument("--mem-fraction", type=float, default=0.5,
                   help="skip a cell whose estimated dense tables exceed this fraction of the physical memory")
    m = sub.choices["memory"]
    m.add_argument("--case", required=True); m.add_argument("--method", choices=["port", "old"], required=True)
    m.add_argument("--kind", choices=["apfc", "rpfc"], required=True); m.add_argument("--checks", type=int, default=0)
    a = ap.parse_args(argv)
    if a.cmd == "memory":
        memory_child(a)
    else:
        run_series(a)


if __name__ == "__main__":
    main()
