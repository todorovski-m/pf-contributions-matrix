function mpc = pfc_crafted_branch(pf, pt, qf, qt)
%PFC_CRAFTED_BRANCH  Solved case9 with one extra branch 4-5 (parallel to branch 4-5, R = 0, X = 1e9, B = 0, so
%   that its series current is negligible) whose terminal flows are set to PF = pf, QF = qf and the delivered flows
%   pt, qt (stored with the MATPOWER sign: PT = -pt, QT = -qt). The extra branch is the last one (index 10).
%   For tests of the terminal classification with OPTS.solved = true.
define_constants;
mpc = runpf(loadcase('cases/case9'), mpoption('verbose', 0, 'out.all', 0));
r = mpc.branch(2, :);                        % branch 4-5
r([BR_R BR_X BR_B]) = [0 1e9 0];
r([PF QF PT QT]) = [pf qf -pt -qt];
mpc.branch(end + 1, :) = r;
end
