% apfc.m or rpfc.m (KIND) on the exported flows of case CNAME (benchmark/inputs_large/<case>.mat, the
% output of my_runpf.m with the delivered PT and QT): the solved struct is built from the exported data with PT and
% QT restored to the MATPOWER sign convention, and run with opts.solved = true. Writes results/<case>_<kind>.csv
% (writetable of tests/octave_shims, %.17g) in the current folder. CNAME, KIND and INPUT are set by the caller.
define_constants;
r = load(INPUT);
mpc = struct('version', '2', 'baseMVA', r.baseMVA, 'bus', r.bus, 'gen', r.gen, 'branch', r.branch, 'success', r.success);
mpc.branch(:, [PT QT]) = -mpc.branch(:, [PT QT]);
t0 = tic;
feval(KIND, mpc, struct('solved', true, 'write', true, 'name', ['results/' CNAME]));
printf('%s %s finished in Octave\n', CNAME, KIND);
