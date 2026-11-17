function prepare_case(name, out_mat, map_csv)
%PREPARE_CASE  Prepare a MATPOWER case for the export (benchmark/timing/export_large.sh), in a scratch folder:
%   renumber the buses 1..nb in the order of the bus table (gen and branch bus numbers mapped accordingly), and
%   strip the matrices to the standard MATPOWER columns (bus 13, gen 21, branch 13). Saves the case to OUT_MAT
%   (savecase, MAT format: no loss of precision) and the mapping to MAP_CSV (old_bus,new_bus). No MATPOWER file is
%   changed. Physics-neutral; checked by neutrality_check.m.
define_constants;
mpc = loadcase(name);
old = mpc.bus(:, BUS_I);
nb = size(mpc.bus, 1);
e2n = sparse(old, ones(nb, 1), (1:nb)', max(old), 1);
mpc.bus(:, BUS_I) = (1:nb)';
mpc.gen(:, GEN_BUS) = full(e2n(mpc.gen(:, GEN_BUS)));
mpc.branch(:, F_BUS) = full(e2n(mpc.branch(:, F_BUS)));
mpc.branch(:, T_BUS) = full(e2n(mpc.branch(:, T_BUS)));
if any(mpc.gen(:, GEN_BUS) == 0) || any(mpc.branch(:, F_BUS) == 0) || any(mpc.branch(:, T_BUS) == 0)
  error('prepare_case:map', 'a generator or branch refers to a bus that is not in the bus table');
end
if isfield(mpc, 'areas') && ~isempty(mpc.areas)
  mpc.areas(:, 2) = full(e2n(mpc.areas(:, 2)));     % price reference bus of the area
end
mpc.bus = mpc.bus(:, 1:min(13, size(mpc.bus, 2)));
mpc.gen = mpc.gen(:, 1:min(21, size(mpc.gen, 2)));
mpc.branch = mpc.branch(:, 1:min(13, size(mpc.branch, 2)));
savecase(out_mat, mpc);
fid = fopen(map_csv, 'w'); fprintf(fid, 'old_bus,new_bus\n'); fprintf(fid, '%d,%d\n', [old'; 1:nb]); fclose(fid);
end
