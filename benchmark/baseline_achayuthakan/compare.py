"""Comparison of the Achayuthakan et al. (2010) baseline with APFC at demand level (brief T14h, step 3). No timing.

For each network: the contribution of every generator bus to the active power of every bus demand (rows D_<bus>),
- by APFC: the port (python/pfc_matrix.apfc, checks off), rows D_<bus>; its generator columns G_<bus>, G_<bus>_2, ...
  are summed per bus, because the baseline is at bus level (each bus has at most one demand row, checked);
- by the baseline, net-downstream and gross-upstream (achayuthakan.py), sparse solution; the explicit inverse is
  compared with it.
Reported: the largest absolute difference (MW), the largest difference relative to the demand of its row, the mean
absolute difference over all cells, and the generator and demand of the largest difference; the residuals (row sum
against the demand) of each method; the evidence for the hypotheses on the differences; the matrix order (n_b, n_e),
nnz, density and the 1-norm condition estimate of each system matrix; and the peak resident memory of each run, in a
fresh process per run (resource.getrusage, as the timing harness's memory measurement).
Inputs: the exported power flows of the final timing series (benchmark/inputs, benchmark/inputs_large, through the
harness's loader). Outputs (folder given): comparison.csv, comparison_paper_cases.tex, summary.md.

Usage (from the matrix root):
    python -B benchmark/baseline_achayuthakan/compare.py <output dir> [case ...]
    python -B benchmark/baseline_achayuthakan/compare.py --memory <case> <method>   (one run in this process; prints the peak)
"""
import csv
import json
import os
import re
import resource
import subprocess
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(MATRIX, "benchmark", "timing"))
sys.path.insert(0, os.path.join(MATRIX, "python"))
import numpy as np    # noqa: E402
import pandas as pd   # noqa: E402

CASES = ["case1", "case9", "case14", "case30", "case118", "case_ACTIVSg500", "case_ACTIVSg2000"]
PAPER = ["case1", "case9", "case30", "case118"]
METHODS = ["apfc", "net-inverse", "net-sparse", "gross-inverse", "gross-sparse"]
g = lambda x: "" if x is None else f"{x:.17g}"


def load(name):
    import cases as C
    return C.load(name, C.input_dirs())


def apfc_bus_table(data):
    """APFC (port, checks off): (demand rows D_<bus> with generator columns summed per bus, info, full table)."""
    from pfc_matrix import Case, Options, read_branch_params, apfc
    case = Case.from_frames(data["bus"], data["gen"], data["branch"], data["case"], params=read_branch_params(data["path"]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t, info = apfc(case, Options(checks=True))
    dem = t.loc[[r for r in t.index if re.fullmatch(r"D_\d+(_\d+)?", r)]]
    base = [re.match(r"(D_\d+)", r).group(1) for r in dem.index]
    multi = len(base) != len(set(base))
    dem = dem.groupby(base, sort=False).sum()
    gcols = {}
    for c in dem.columns:
        m = re.fullmatch(r"(G_\d+)(_\d+)?", c)
        if m:
            gcols.setdefault(m.group(1), []).append(c)
    out = pd.DataFrame({k: dem[v].sum(axis=1) for k, v in gcols.items()}, index=dem.index)
    other = [c for c in dem.columns if not re.fullmatch(r"G_\d+(_\d+)?", c)]
    return out, info, t, multi, other


def run_method(data, method):
    import achayuthakan as A
    if method == "apfc":
        return apfc_bus_table(data)[0]
    kind, how = method.split("-")
    return A.contributions(data["bus"], data["gen"], data["branch"], kind=kind, method=how)


def memory_child(case, method):
    data = load(case)
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    run_method(data, method)
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({"before": before, "after": after, "unit": "bytes" if sys.platform == "darwin" else "KiB"}))


def peak_memory(case, method):
    r = subprocess.run([sys.executable, "-B", os.path.abspath(__file__), "--memory", case, method],
                       capture_output=True, text=True, timeout=1800)
    return json.loads(r.stdout.strip().splitlines()[-1])


def diff(base, ap, demand):
    """Differences on the common rows and generator columns (missing cells are 0)."""
    rows = list(ap.index); cols = sorted(set(ap.columns) | set(base.columns), key=lambda c: int(c[2:]))
    B = base.reindex(index=rows, columns=cols, fill_value=0.0).to_numpy()
    P = ap.reindex(index=rows, columns=cols, fill_value=0.0).to_numpy()
    d = np.abs(B - P)
    r, c = np.unravel_index(int(np.argmax(d)), d.shape)
    rel = d / demand[:, None]
    return {"max_abs": float(d.max()), "at": f"{cols[c]} -> {rows[r]}", "max_rel": float(rel.max()),
            "mean_abs": float(d.mean()), "rows_missing": [x for x in rows if x not in base.index]}


def main(outd, cases):
    import achayuthakan as A
    os.makedirs(outd, exist_ok=True)
    rows = []
    for c in cases:
        data = load(c)
        ap, info_ap, full, multi, other = apfc_bus_table(data)
        bus = data["bus"].set_index("BUS_I")
        demand = np.array([float(bus.loc[int(r[2:]), "PD"]) for r in ap.index])
        res_ap = float(np.abs(ap.sum(axis=1).to_numpy() - demand).max())
        mem = {m: peak_memory(c, m) for m in METHODS}
        rec = {"case": c, "apfc_demand_rows": len(ap), "apfc_bus_with_several_demand_rows": multi,
               "apfc_other_source_columns": " ".join(other), "apfc_residual": res_ap,
               "n_e": info_ap["n_e"], "apfc_nnz_A": info_ap["nnz_A"], "apfc_density_A": info_ap["nnz_A"] / info_ap["n_e"] ** 2,
               "apfc_condest_A_1norm": info_ap["condest_A"], "apfc_peak_rss": mem["apfc"]["after"],
               "rss_unit": mem["apfc"]["unit"]}
        for kind in ("net", "gross"):
            ts, info = A.contributions(data["bus"], data["gen"], data["branch"], kind=kind, method="sparse", return_info=True)
            ti = A.contributions(data["bus"], data["gen"], data["branch"], kind=kind, method="inverse")
            dm = ts.loc[[r for r in ts.index if r.startswith("D_")]]
            dd = diff(dm, ap, demand)
            rec.update({
                f"{kind}_max_abs_diff_MW": dd["max_abs"], f"{kind}_max_rel_diff": dd["max_rel"], f"{kind}_mean_abs_diff_MW": dd["mean_abs"],
                f"{kind}_at": dd["at"], f"{kind}_rows_missing": len(dd["rows_missing"]),
                f"{kind}_residual": float(np.abs(dm.reindex(ap.index).sum(axis=1).to_numpy() - demand).max()),
                f"{kind}_inverse_vs_sparse": float(np.abs(ti.to_numpy() - ts.to_numpy()).max()),
                f"{kind}_n_b": info["n_b"], f"{kind}_nnz_A": info["nnz_A"], f"{kind}_density_A": info["density"],
                f"{kind}_condest_A_1norm": A.condest_1norm(info["A"]),
                f"{kind}_peak_rss_inverse": mem[f"{kind}-inverse"]["after"], f"{kind}_peak_rss_sparse": mem[f"{kind}-sparse"]["after"],
                f"{kind}_branch_demands": info["branch_demands"]})
            if kind == "gross":
                # hypothesis 1 evidence: the excess of each demand over its demand is the loss the gross variant
                # allocates to it; summed over all demands (bus and terminal) it equals the network losses
                br = data["branch"]
                losses = float((br["PF"] - br["PT"]).sum())
                rec["network_losses_MW"] = losses
                rec["gross_allocated_losses_MW"] = float(ts.to_numpy().sum() - info["demand_values"].sum())
                rec["gross_max_excess_over_demand_MW"] = float((dm.sum(axis=1).to_numpy()
                                                                - np.array([float(bus.loc[int(r[2:]), "PD"]) for r in dm.index])).max())
        rows.append(rec)
        print(f"{c}: net {g(rec['net_max_abs_diff_MW'])}, gross {g(rec['gross_max_abs_diff_MW'])}", flush=True)
    with open(os.path.join(outd, "comparison.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({k: (g(v) if isinstance(v, float) else v) for k, v in r.items()})
    write_tex(rows, outd)
    write_summary(rows, outd)


def write_tex(rows, outd):
    by = {r["case"]: r for r in rows}
    sci = lambda x: "0" if x == 0 else f"{x:.1e}".replace("e-", "\\text{e-}").replace("e+", "\\text{e}")
    tex = ["% Generated by benchmark/baseline_achayuthakan/compare.py (do not edit). Demand-level comparison of the",
           "% Achayuthakan et al. (2010) baseline with APFC: largest absolute difference (MW) of a generator-to-demand",
           "% contribution, matrix orders n_b (baseline) and n_e (APFC), and 1-norm condition estimates.",
           "\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
           "Case & $n_b$ & $n_e$ & net--downstream (MW) & gross--upstream (MW) & cond$_1$ baseline (net) & cond$_1$ APFC \\\\", "\\midrule"]
    for c in PAPER:
        if c in by:
            r = by[c]
            tex.append(f"{c} & {r['net_n_b']} & {r['n_e']} & ${sci(r['net_max_abs_diff_MW'])}$ & ${sci(r['gross_max_abs_diff_MW'])}$ & "
                       f"{r['net_condest_A_1norm']:.1f} & {r['apfc_condest_A_1norm']:.1f} \\\\")
    tex += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(outd, "comparison_paper_cases.tex"), "w").write("\n".join(tex) + "\n")


def write_summary(rows, outd):
    L = ["# Achayuthakan et al. (2010) baseline against APFC, demand level", "",
         "Generated by `compare.py`. Differences in MW over the cells generator bus x bus demand (rows D_<bus>); residual =",
         "largest |row sum - demand|; condition estimates in the 1-norm; peak memory = ru_maxrss of one run in a fresh process.", "",
         "| case | n_b | n_e | net: max abs diff | net: max rel diff | gross: max abs diff | gross: max rel diff | gross: at | residual APFC | residual net | residual gross | network losses | losses allocated by gross |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['case']} | {r['net_n_b']} | {r['n_e']} | {g(r['net_max_abs_diff_MW'])} | {g(r['net_max_rel_diff'])} | {g(r['gross_max_abs_diff_MW'])} | "
                 f"{g(r['gross_max_rel_diff'])} | {r['gross_at']} | {g(r['apfc_residual'])} | {g(r['net_residual'])} | {g(r['gross_residual'])} | "
                 f"{g(r['network_losses_MW'])} | {g(r['gross_allocated_losses_MW'])} |")
    L += ["", "| case | matrix | order | nnz | density | cond_1 estimate | peak memory (bytes) |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['case']} | APFC A | {r['n_e']} | {r['apfc_nnz_A']} | {g(r['apfc_density_A'])} | {g(r['apfc_condest_A_1norm'])} | {r['apfc_peak_rss']} |")
        for k in ("net", "gross"):
            L.append(f"| {r['case']} | baseline {k} | {r[f'{k}_n_b']} | {r[f'{k}_nnz_A']} | {g(r[f'{k}_density_A'])} | {g(r[f'{k}_condest_A_1norm'])} | "
                     f"inverse {r[f'{k}_peak_rss_inverse']}, sparse {r[f'{k}_peak_rss_sparse']} |")
    open(os.path.join(outd, "summary.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "--memory":
        memory_child(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1], sys.argv[2:] or CASES)
