"""Figures of a timing results folder (after summarize.py): writes timing_loglog.pdf/.png and timing_ratio.pdf/.png into
the folder. No text for a manuscript is produced.

- timing_loglog: median time (s) against n_e, log-log, one panel per kind (APFC, RPFC), the port with checks off and
  the old method. A cell where the old method did not finish is drawn at the time-out of the run (environment.json,
  args.timeout) with a marker per failure kind (time-out, iteration cap, error) and a label of the kind; the y-position
  of these markers is the time-out, not a measured time.
- timing_ratio: the ratio old/port of the medians (checks off) against n_e, log-log, both kinds.
A smoke folder carries the label "smoke, not for publication" in every panel.

Usage: python -B benchmark/timing/plot_timing.py <results folder>
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

FAIL = {"time-out": ("X", "did not finish (time-out)"), "iteration cap": ("^", "did not finish (iteration cap)"),
        "error": ("s", "did not finish (error)")}


def fail_kind(status):
    for k in FAIL:
        if f"did not finish ({k}" in status:
            return k
    return None


def load(folder):
    env = json.load(open(os.path.join(folder, "environment.json")))
    rows = list(csv.DictReader(open(os.path.join(folder, "summary.csv"))))
    return env, rows


def plot(folder):
    env, rows = load(folder)
    smoke = bool(env.get("smoke"))
    timeout = float((env.get("args") or {}).get("timeout") or 0) or None
    title_suffix = " (smoke, not for publication)" if smoke else ""
    ok = lambda r: r["status"].split(" (slow)")[0] == "ok" and r["n_e"] and r["median_ns"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, kind in zip(axes, ("apfc", "rpfc")):
        for method, checks, mk, lab in (("port", "off", "o-", "matrix method (port), checks off"), ("old", "", "s--", "old iterative method")):
            pts = sorted((float(r["n_e"]), float(r["median_ns"]) / 1e9) for r in rows
                         if r["kind"] == kind and r["method"] == method and r["checks"] == checks and ok(r))
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, ms=4, lw=1, label=lab)
        if timeout:
            used = set()
            for r in rows:
                k = fail_kind(r["status"]) if r["kind"] == kind and r["method"] == "old" else None
                if k and r["n_e"]:
                    m, lab = FAIL[k]
                    ax.plot(float(r["n_e"]), timeout, m, color="tab:red", ms=7, ls="none",
                            label=None if k in used else lab + " (drawn at the time-out)")
                    used.add(k)
            ax.axhline(timeout, color="tab:red", lw=0.6, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("n_e (elements of the expanded network)")
        ax.set_title(kind.upper() + title_suffix, fontsize=9)
        ax.grid(True, which="both", lw=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("median time per run (s)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(folder, f"timing_loglog.{ext}"), dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for kind, mk in (("apfc", "o-"), ("rpfc", "s--")):
        pts = sorted((float(r["n_e"]), float(r["ratio_old_over_port_off"])) for r in rows
                     if r["kind"] == kind and r["method"] == "old" and r.get("ratio_old_over_port_off") and r["n_e"])
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], mk, ms=4, lw=1, label=kind.upper())
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("n_e (elements of the expanded network)")
    ax.set_ylabel("ratio of medians, old / port (checks off)")
    ax.set_title("ratio old/port" + title_suffix, fontsize=9)
    ax.grid(True, which="both", lw=0.3)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(folder, f"timing_ratio.{ext}"), dpi=200)
    plt.close(fig)
    return [os.path.join(folder, f"timing_{n}.{e}") for n in ("loglog", "ratio") for e in ("pdf", "png")]


if __name__ == "__main__":
    for f in plot(sys.argv[1]):
        print("wrote", f)
