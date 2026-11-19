"""Replace the first comment block written by savecase in benchmark/reconstruction/cases/<case>.m by a header that
describes the change (brief T14g, part C). Usage: python add_header.py <case file> ..."""
import sys
HEAD = {
    "syn9_pst": ["Phase-shifter test case for the reconstruction check (brief T14g, part C), derived from case9.",
                 "Change: SHIFT = 5 degrees on branch 8 (8-9; R 0.032, X 0.161, B 0.306 pu), which carries 86.6 MW in",
                 "case9; TAP stays 0 (ratio 1). Everything else as in case9. Written with savecase (make_pst_cases.m)."],
    "syn9_pst_tap": ["Phase-shifter test case with an off-nominal tap (brief T14g, part C), derived from case9.",
                     "Change: SHIFT = 5 degrees and TAP = 1.05 on branch 8 (8-9; R 0.032, X 0.161, B 0.306 pu).",
                     "Everything else as in case9. Written with savecase (make_pst_cases.m)."],
}
for path in sys.argv[1:]:
    name = path.rsplit("/", 1)[-1][:-2]
    lines = open(path).read().splitlines()
    assert lines[0].startswith("function mpc = "), path
    i = 1
    while i < len(lines) and lines[i].startswith("%"):
        i += 1
    head = [f"%{name.upper()}  " + HEAD[name][0]] + ["%   " + l for l in HEAD[name][1:]]
    open(path, "w").write("\n".join([lines[0]] + head + lines[i:]) + "\n")
    print("header written:", path)
