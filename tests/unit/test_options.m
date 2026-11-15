% Unit tests for the options and the input forms of apfc.m and rpfc.m (brief T12b, Step 1).
% Run with: zsh tests/octave_in_scratch.sh test_options
define_constants;
o = pfc_checks('options', struct());
pfc_report('options: defaults', o.tol == 1e-6 && o.tol_div == 1e-9 && o.write && ~o.solved && isempty(o.name));
[ok, msg] = pfc_expect_error(@() pfc_checks('options', struct('tolerance', 1)), 'pfc:options:unknownField');
pfc_report('options: unknown field is an error', ok, msg);
[ok, msg] = pfc_expect_error(@() pfc_checks('options', struct('tol', -1)), 'pfc:options:value');
pfc_report('options: negative tol is an error', ok, msg);
[ok, msg] = pfc_expect_error(@() apfc(42), 'pfc:input:casedata');
pfc_report('input: numeric casedata is an error', ok, msg);
[ok, msg] = pfc_expect_error(@() apfc(loadcase('cases/case9'), struct('solved', true)), 'pfc:input:notSolved');
pfc_report('input: unsolved struct with opts.solved is an error', ok, msg);

% Acceptance 4: opts.solved = true on the output of one runpf call gives the same tables as the default call.
cases = {'case1', 'case9', 'case14', 'case30', 'case33', 'case118', 'tap_loop'};
for k = 1:numel(cases)
    c = cases{k};
    r = runpf(loadcase(['cases/' c]), mpoption('verbose', 0, 'out.all', 0));
    [ta, ~, da] = apfc(['cases/' c], struct('write', false));
    [tb, ~, db] = apfc(r, struct('solved', true));
    [tc] = apfc(loadcase(['cases/' c]));                 % unsolved struct, no name: no file
    pfc_report(sprintf('solved struct: apfc %s identical tables', c), isequal(ta, tb) && isequal(ta, tc));
    [ta, ~, da] = rpfc(['cases/' c], struct('write', false));
    [tb, ~, db] = rpfc(r, struct('solved', true));
    [tc] = rpfc(loadcase(['cases/' c]));
    pfc_report(sprintf('solved struct: rpfc %s identical tables', c), isequal(ta, tb) && isequal(ta, tc));
end
d = dir('results/*.csv');
pfc_report('write: no file written for opts.write = false or for a struct without a name', isempty(d), ...
    sprintf('%d files written', numel(d)));
apfc(loadcase('cases/case9'), struct('name', 'cases/case9_named'));
pfc_report('write: struct input with opts.name writes results/<name>_apfc.csv', exist('results/case9_named_apfc.csv', 'file') == 2);
% diag fields
[~, ~, d1] = rpfc('cases/case9', struct('write', false));
[~, ~, d2] = rpfc('cases/case9', struct('write', false));
need = {'n_bus', 'n_branch', 'n_demand', 'n_source', 'n_e', 'nnz_N', 'nnz_A', 'issparse', 'condest_A', ...
        'attenuation_range', 'removed', 'flipped_branches', 'converted', 'opts'};
pfc_report('diag: all fields present', all(isfield(d1, need)));
pfc_report('diag: condest of A is the same in two calls', d1.condest_A == d2.condest_A, ...
    sprintf('%.17g vs %.17g', d1.condest_A, d2.condest_A));
