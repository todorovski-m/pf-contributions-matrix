function make_pst_cases(outdir)
%MAKE_PST_CASES  Write the phase-shifter test cases syn9_pst.m and syn9_pst_tap.m into OUTDIR (brief T14g, part C).
%   All are case9 (MATPOWER 7.1) with a change on one line:
%     syn9_pst              branch 8-9 (R 0.032, X 0.161, B 0.306 pu; 86.6 MW in case9): SHIFT = 5 degrees, TAP 0
%     syn9_pst_tap          branch 8-9: SHIFT = 5 degrees and TAP = 1.05
%     syn9_pst_carrier      branch 4-5 (R 0.017, X 0.092, B 0.158 pu): SHIFT = 5 degrees, TAP 0; the shifted branch is
%                           a reactive carrier (pattern (+,+)), brief T16 S4a
%     syn9_pst_carrier_tap  branch 9-4 (R 0.01, X 0.085, B 0.176 pu): SHIFT = 5 degrees and TAP = 1.05; the shifted
%                           branch is a reactive carrier after the flip (pattern (-,-)), brief T16 S4a
%   (found by search_pst_carrier.m). Everything else as in case9. Written with savecase; the header is added by
%   add_header.py.
define_constants;
for v = {{'syn9_pst', 0, 5, 8, 9}, {'syn9_pst_tap', 1.05, 5, 8, 9}, {'syn9_pst_carrier', 0, 5, 4, 5}, {'syn9_pst_carrier_tap', 1.05, 5, 9, 4}}
  name = v{1}{1}; tap = v{1}{2}; shift = v{1}{3};
  mpc = loadcase('case9');
  k = find(mpc.branch(:, F_BUS) == v{1}{4} & mpc.branch(:, T_BUS) == v{1}{5});
  assert(numel(k) == 1);
  mpc.branch(k, TAP) = tap;
  mpc.branch(k, SHIFT) = shift;
  savecase(fullfile(outdir, [name '.m']), mpc);
end
end
