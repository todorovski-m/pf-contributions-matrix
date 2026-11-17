"""Scan of the cases for the timing study: sizes, taps and line charging, the branches deleted by my_runpf.m, the
reactive Cases 2 and 3, and, where the port accepts the input, n_e, nnz(N), acyclicity, spectral radius and
condest_A (port with checks=True, phases (a) and (b) and condest only: no dense result table is built), and an
estimate of the memory of the dense result tables of one run.

Usage (from matrix/): python -B benchmark/timing/scan.py --cases c1,c2 [--raw-counts CSV] [--syn-inputs DIR]
    [--out-md FILE] [--out-csv FILE]
"""
import argparse, csv, os, sys, warnings
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np                      # noqa: E402
import cases as C                       # noqa: E402
import timing                           # noqa: E402,F401  (sets the thread counts before numpy is used)

TOL = 1e-6
g = lambda x: "" if x is None else (f"{x:.17g}" if isinstance(x, float) else str(x))


def reactive_classes(branch):
    """Case 2 and Case 3 counts and the carrier mask (in the rows of the exported table), by the rules of rpfc.m:
    deletion of |QF| < tol and |QT| < tol, flip of both negative, Case 3 = QF < 0 and QT > 0, Case 2 = QF > 0 and QT < 0."""
    qf, qt = branch["QF"].to_numpy(float), branch["QT"].to_numpy(float)
    keep = ~((np.abs(qf) < TOL) & (np.abs(qt) < TOL))
    flip = (qf < 0) & (qt < 0)
    qf2, qt2 = np.where(flip, -qt, qf), np.where(flip, -qf, qt)
    c3 = keep & (qf2 < 0) & (qt2 > 0)
    c2 = keep & (qf2 > 0) & (qt2 < 0)
    return int(c2.sum()), int(c3.sum()), keep & ~c2 & ~c3


def scan_case(name, dirs, raw):
    from pfc_matrix import Case, Options, read_branch_params, prepare_active, assemble_active, prepare_reactive, assemble_reactive
    from pfc_matrix.checks import condest
    try:
        data = C.load(name, dirs)
    except Exception as e:      # load_results of the old pipeline cannot read the export
        import scipy.io as sio
        m = sio.loadmat(C.find_input(name, dirs))
        return {"case": name, "buses": m["bus"].shape[0], "branches": m["branch"].shape[0], "generators": m["gen"].shape[0],
                "port": f"not loaded: load_results raises {type(e).__name__}: {str(e)[:120]} "
                        f"(bus, gen, branch columns in the export: {m['bus'].shape[1]}, {m['gen'].shape[1]}, {m['branch'].shape[1]})"}
    bus, gen, branch = data["bus"], data["gen"], data["branch"]
    params = read_branch_params(data["path"])
    tap = branch["TAP"].to_numpy(float); B = params.frame["BR_B"].to_numpy(float)
    both = (tap != 1) & (B != 0)
    n2, n3, carrier = reactive_classes(branch)
    row = {"case": name, "buses": len(bus), "branches": len(branch), "generators": len(gen),
           "buses_numbered_1_to_nb": bool(np.array_equal(bus["BUS_I"].to_numpy(int), np.arange(1, len(bus) + 1))),
           "tap_ne_1": int((tap != 1).sum()), "br_b_ne_0": int((B != 0).sum()), "tap_and_br_b": int(both.sum()),
           "tap_and_br_b_carriers": int((both & carrier).sum()),
           "deleted_by_my_runpf": (raw[name] - len(branch)) if name in raw else None,
           "case2_branches": n2, "case3_branches": n3}
    try:
        case = Case.from_frames(bus, gen, branch, name, params=params)
    except Exception as e:
        row["port"] = f"rejected: {getattr(e, 'identifier', type(e).__name__)}"
        return row
    row["port"] = "accepted"
    pa = pr = None
    try:      # buses with an outgoing branch and a total active inflow below tol_div (the APFC zero-inflow guard)
        p0 = prepare_active(case, TOL)
        inflow = np.zeros(p0.nb + 1)
        np.add.at(inflow, p0.t, p0.Pt); np.add.at(inflow, p0.g, p0.Pg)
        has_out = np.zeros(p0.nb + 1, bool); has_out[p0.f] = True
        row["zero_inflow_buses_apfc"] = int((has_out & (inflow < 1e-9)).sum())
    except Exception as e:
        row["zero_inflow_buses_apfc"] = f"not computed ({getattr(e, 'identifier', type(e).__name__)})"
    for kind, prep_f, asm in (("apfc", prepare_active, lambda p: assemble_active(p, 1e-9, True)),
                              ("rpfc", prepare_reactive, lambda p: assemble_reactive(p, TOL, 1e-9))):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = prep_f(case, TOL)
                if kind == "apfc": pa = p
                else: pr = p
                sy = asm(p)
        except Exception as e:
            row[f"error_{kind}"] = f"{getattr(e, 'identifier', type(e).__name__)}: {str(e)[:160]}"
            continue
        row.update({f"n_e_{kind}": p.ne, f"nnz_N_{kind}": int(np.count_nonzero(sy.N.data)), f"acyclic_{kind}": sy.acyclic,
                    f"spectral_radius_{kind}": sy.spectral_radius, f"condest_A_{kind}": condest(sy.A)})
        if kind == "rpfc":
            row["negative_series_flow_carriers"] = len(sy.negative_series_flow)
    if pa is None or pr is None:
        return row
    from timing import dense_bytes          # the estimate used by the harness to skip cells
    row["dense_table_bytes_apfc"] = dense_bytes(pa, "apfc")
    row["dense_table_bytes_rpfc"] = dense_bytes(pr, "rpfc")
    return row


def usability(row, limit):
    """usable for both kinds / for apfc only / for rpfc only / not usable, with the reasons: a port error, or estimated
    dense tables above the memory limit of the harness."""
    if row.get("port") != "accepted":
        return f"not usable ({row.get('port')})"
    ok, why = {}, {}
    for k in ("apfc", "rpfc"):
        err = row.get(f"error_{k}")
        size = row.get(f"dense_table_bytes_{k}")
        if err:
            ok[k], why[k] = False, err.split(":")[0] + ":" + err.split(":")[1] + ":" + err.split(":")[2].split(" ")[0] if err.count(":") >= 2 else err
        elif size is not None and size > limit:
            ok[k], why[k] = False, f"dense tables {size} bytes > {int(limit)} bytes"
        else:
            ok[k] = True
    if ok["apfc"] and ok["rpfc"]:
        return "usable for both kinds"
    if ok["apfc"]:
        return f"usable for apfc only (rpfc: {why['rpfc']})"
    if ok["rpfc"]:
        return f"usable for rpfc only (apfc: {why['apfc']})"
    return f"not usable (apfc: {why['apfc']}; rpfc: {why['rpfc']})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True); ap.add_argument("--raw-counts"); ap.add_argument("--syn-inputs")
    ap.add_argument("--large-inputs"); ap.add_argument("--out-md"); ap.add_argument("--out-csv")
    ap.add_argument("--mem-fraction", type=float, default=0.5)
    a = ap.parse_args()
    raw = {}
    if a.raw_counts:
        for r in csv.DictReader(open(a.raw_counts)):
            raw[r["case"]] = int(r["branches_raw"])
    dirs = C.input_dirs(a.syn_inputs, a.large_inputs)
    rows = [scan_case(c, dirs, raw) for c in a.cases.split(",")]
    limit = a.mem_fraction * timing.physical_memory()
    for r in rows:
        r["status"] = usability(r, limit)
    cols = ["case", "status", "buses", "branches", "generators", "buses_numbered_1_to_nb", "tap_ne_1", "br_b_ne_0", "tap_and_br_b",
            "tap_and_br_b_carriers", "deleted_by_my_runpf", "case2_branches", "case3_branches", "port", "zero_inflow_buses_apfc",
            "negative_series_flow_carriers", "n_e_apfc", "nnz_N_apfc", "acyclic_apfc", "spectral_radius_apfc", "condest_A_apfc",
            "n_e_rpfc", "nnz_N_rpfc", "acyclic_rpfc", "spectral_radius_rpfc", "condest_A_rpfc",
            "dense_table_bytes_apfc", "dense_table_bytes_rpfc", "error_apfc", "error_rpfc"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(g(r.get(k)) for k in cols) + " |")
    text = "\n".join(lines)
    print(text)
    if a.out_md:
        open(a.out_md, "w").write(text + "\n")
    if a.out_csv:
        with open(a.out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); [w.writerow(r) for r in rows]


if __name__ == "__main__":
    main()
