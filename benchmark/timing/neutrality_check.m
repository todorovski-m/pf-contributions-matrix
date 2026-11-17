function neutrality_check(name, prepared_mat, map_csv, out_csv)
%NEUTRALITY_CHECK  Power flow of the original case NAME against the power flow of the prepared case (renumbered and
%   stripped, PREPARED_MAT), mapped back through MAP_CSV, both with runpf and the same options. Appends one line to
%   OUT_CSV: case, the maximum absolute differences of VM, VA (degrees), PF, QF, PT, QT, the number of iterations of
%   both power flows, both success flags, and whether the mapping is a bijection onto 1..nb.
define_constants;
opt = mpoption('verbose', 0, 'out.all', 0);
r0 = runpf(loadcase(name), opt);
r1 = runpf(loadcase(prepared_mat), opt);
m = dlmread(map_csv, ',', 1, 0);
nb = size(r0.bus, 1);
bij = isequal(sort(m(:, 2))', 1:nb) && numel(unique(m(:, 1))) == nb;
% rows are in the same order in both (renumbering keeps the table order); check that the bus numbers map back
back = isequal(m(r1.bus(:, BUS_I), 1), r0.bus(:, BUS_I)) && isequal(m(r1.branch(:, F_BUS), 1), r0.branch(:, F_BUS)) ...
       && isequal(m(r1.branch(:, T_BUS), 1), r0.branch(:, T_BUS)) && isequal(m(r1.gen(:, GEN_BUS), 1), r0.gen(:, GEN_BUS));
d = @(c, M) max(abs(r0.(M)(:, c) - r1.(M)(:, c)));
fid = fopen(out_csv, 'a');
fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%d\n', name, d(VM, 'bus'), d(VA, 'bus'), d(PF, 'branch'), ...
        d(QF, 'branch'), d(PT, 'branch'), d(QT, 'branch'), r0.iterations, r1.iterations, r0.success, r1.success, bij, back);
fclose(fid);
end
