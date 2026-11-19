function make_pst_cases(outdir)
%MAKE_PST_CASES  Write the phase-shifter test cases syn9_pst.m and syn9_pst_tap.m into OUTDIR (brief T14g, part C).
%   Both are case9 (MATPOWER 7.1) with a change on branch 8 (8-9; R 0.032, X 0.161, B 0.306 pu), which carries
%   86.6 MW in case9:
%     syn9_pst      SHIFT = 5 degrees (TAP stays 0, that is ratio 1)
%     syn9_pst_tap  SHIFT = 5 degrees and TAP = 1.05
%   Everything else as in case9. Written with savecase; the header is added by add_header.py.
define_constants;
for v = {{'syn9_pst', 0, 5}, {'syn9_pst_tap', 1.05, 5}}
  name = v{1}{1}; tap = v{1}{2}; shift = v{1}{3};
  mpc = loadcase('case9');
  k = find(mpc.branch(:, F_BUS) == 8 & mpc.branch(:, T_BUS) == 9);
  assert(numel(k) == 1);
  mpc.branch(k, TAP) = tap;
  mpc.branch(k, SHIFT) = shift;
  savecase(fullfile(outdir, [name '.m']), mpc);
end
end
