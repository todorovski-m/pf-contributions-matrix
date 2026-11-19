function search_pst_carrier()
%SEARCH_PST_CARRIER  Search case9 variants with a phase shift on one line in which the shifted branch is a reactive
%   carrier (brief T16, S4a): for every branch with line charging and shifts of -10 to 10 degrees (and tap 1 or 1.05),
%   run the power flow and print the variants whose shifted branch has the reactive pattern (+,+) or (-,-) of the
%   from-end flow QF and the delivered to-end flow -QT (tol 1e-6). Prints Markdown table rows.
define_constants;
base = loadcase('case9');
opt = mpoption('verbose', 0, 'out.all', 0);
fprintf('| branch | tap | shift (degrees) | converged | PF | QF | QT delivered | pattern |\n|---|---|---|---|---|---|---|---|\n');
for k = find(base.branch(:, BR_B) ~= 0)'
  for tap = [0 1.05]
    for sh = [-10 -5 -2 -1 1 2 5 10]
      mpc = base; mpc.branch(k, TAP) = tap; mpc.branch(k, SHIFT) = sh;
      r = runpf(mpc, opt);
      qf = r.branch(k, QF); qt = -r.branch(k, QT);
      s = @(x) char('+' * (x > 1e-6) + '-' * (x < -1e-6) + '0' * (abs(x) <= 1e-6));
      if r.success && ((qf > 1e-6 && qt > 1e-6) || (qf < -1e-6 && qt < -1e-6))
        fprintf('| %d-%d | %g | %g | %d | %.17g | %.17g | %.17g | (%s,%s) |\n', base.branch(k, F_BUS), base.branch(k, T_BUS), ...
                tap, sh, r.success, r.branch(k, PF), qf, qt, s(qf), s(qt));
      end
    end
  end
end
end
