function rep = pfc_replace_both_ends(mpc, k)
%PFC_REPLACE_BOTH_ENDS  The equivalent of a solved case in which branch K draws active power from both of its
%   buses (pattern (+,-)): the branch is deleted and replaced by two loads, generators with negative active output
%   at its from bus (PG = -PF) and at its to bus (PG = -PT, MATPOWER sign: the power the branch draws there), with
%   QG = 0; every other flow and voltage is kept. For the APFC equivalence test of brief T14c.
define_constants;
g = repmat(mpc.gen(1, :), 2, 1);
g(:, GEN_BUS) = mpc.branch(k, [F_BUS; T_BUS]);
g(:, PG) = -[mpc.branch(k, PF); mpc.branch(k, PT)];
g(:, QG) = 0; g(:, GEN_STATUS) = 1;
rep = mpc;
rep.branch(k, :) = [];
rep.gen = [rep.gen; g];
if isfield(rep, 'gencost'), rep.gencost = [rep.gencost; repmat(rep.gencost(1, :), 2, 1)]; end
end
