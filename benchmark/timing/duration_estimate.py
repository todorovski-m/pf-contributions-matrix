"""T14f step D: bounds of the duration of the final series from its settings, not from a measurement.
Cells: per combination of cases.COMBOS, the port with checks off and on, and the old method if it can finish (by the
predictor; for the small cases the first series' outcome is used: every small combination finishes except tap_loop
RPFC). Settings (run_timing.sh defaults): budget 30 s per cell (at least 5 repetitions, 2 warm-up runs; a slow cell has
3 repetitions), quiet wait at most 300 s per cell, old-method time-out 1800 s per run, stall probe 120 s for predicted
stalls. Prints the counts and the bounds; the parts that depend on run times are named, not estimated.
Usage (from the matrix root): python -B benchmark/timing/duration_estimate.py"""
import sys
sys.path.insert(0, "benchmark/timing")
import cases as C, old_stall_predictor as P
BUDGET, WAIT, TIMEOUT, PROBE = 30, 300, 1800, 120
dirs = C.input_dirs()
n_comb = n_port = n_old = n_stall = n_risk = n_finish_large = 0
for c, kinds in C.COMBOS:
    data = None
    for k in kinds:
        n_comb += 1; n_port += 2
        if c in C.SMALL:
            if not (c == "tap_loop" and k == "rpfc"):
                n_old += 1
            continue
        data = data or C.load(c, dirs)
        r = P.predict(data, k, C.correct_qsf(data)[0])
        if r["label"] == "predicted stall":
            n_stall += 1
        elif r["division_by_zero_risk"]:
            n_risk += 1
        else:
            n_old += 1; n_finish_large += 1
cells = n_port + n_old
print("| quantity | value |\n|---|---|")
for q, v in [("combinations", n_comb), ("port cells (checks off and on)", n_port), ("old-method cells expected to finish", n_old),
             ("of them large", n_finish_large), ("large combinations predicted to stall (guard run: the probe)", n_stall),
             ("large combinations with only a division-by-zero risk (guard run ends with the error)", n_risk),
             ("timed cells", cells),
             ("lower bound, timed cells x budget (s)", cells * BUDGET),
             ("added at most by the quiet wait, timed cells x 300 s (s)", cells * WAIT),
             ("added at most by the probes, predicted stalls x 120 s (s)", n_stall * PROBE),
             ("before T14f: the same guard runs at most, predicted stalls x 1800 s (s)", n_stall * TIMEOUT)]:
    print(f"| {q} | {v} |")
