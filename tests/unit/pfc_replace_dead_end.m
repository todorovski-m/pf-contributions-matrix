function rep = pfc_replace_dead_end(mpc, dead_bus, shunt)
%PFC_REPLACE_DEAD_END  The equivalent of a solved case with one dead-end bus DEAD_BUS: the dead-end bus and its
%   branch are deleted and replaced, at the other bus o of the branch, by the terminal injection the branch
%   represents; every other flow and voltage is kept. With P_o, Q_o the power that the branch draws from bus o
%   (MATPOWER sign: PF, QF if o is the from bus, PT, QT if o is the to bus):
%     a generator at bus o with PG = -P_o (a demand row G_o if P_o > 0);
%     SHUNT false: the same generator has QG = -Q_o (a source G_o if Q_o < 0, a demand row G_o if Q_o > 0);
%     SHUNT true (the dead-end bus carries a bus shunt): QG = 0, and bus o gets a bus shunt of power -Q_o
%       (BS(o) increased by -Q_o / VM(o)^2), the power that the model delivers at bus o.
define_constants;
k = find(mpc.branch(:, F_BUS) == dead_bus | mpc.branch(:, T_BUS) == dead_bus);
assert(numel(k) == 1, 'one branch at the dead-end bus expected');
if mpc.branch(k, F_BUS) == dead_bus
    o = mpc.branch(k, T_BUS); P_o = mpc.branch(k, PT); Q_o = mpc.branch(k, QT);
else
    o = mpc.branch(k, F_BUS); P_o = mpc.branch(k, PF); Q_o = mpc.branch(k, QF);
end
rep = mpc;
rep.branch(k, :) = [];
ib = find(rep.bus(:, BUS_I) == dead_bus);
assert(~any(rep.gen(:, GEN_BUS) == dead_bus), 'no generator at the dead-end bus expected');
rep.bus(ib, :) = [];
g = rep.gen(1, :);
g(GEN_BUS) = o; g(PG) = -P_o; g(GEN_STATUS) = 1;
io = find(rep.bus(:, BUS_I) == o);
if shunt
    g(QG) = 0;
    rep.bus(io, BS) = rep.bus(io, BS) - Q_o / rep.bus(io, VM)^2;
else
    g(QG) = -Q_o;
end
rep.gen(end + 1, :) = g;
if isfield(rep, 'gencost')
    rep.gencost(end + 1, :) = rep.gencost(1, :);
end
end
