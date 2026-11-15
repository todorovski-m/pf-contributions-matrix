% diag output on the shared cases, case24_ieee_rts, tap_loop and the synthetic cases (brief T12b, acceptance 2).
% diag.negative_series_flow of rpfc is checked against the carrier branches with a negative series flow
% computed by the old iterative method (tests/data/negative_series_flow_carriers.csv), for all twelve cases,
% and it must not produce a warning.
% Run with: zsh tests/octave_in_scratch.sh test_diag   (it copies tests/data/*.csv next to the tests)
cases = {'case1', 'case9', 'case14', 'case24_ieee_rts', 'case30', 'case33', 'case118', 'tap_loop', ...
         'syn9_neg_pd', 'syn9_neg_pg', 'syn9_neg_qd', 'syn9_neg_qg', 'syn9_case3_tap'};
% expected carrier sets
fn = 'unit/negative_series_flow_carriers.csv';
lines = strsplit(strtrim(fileread(fn)), sprintf('\n'));
expected = containers.Map(); covered = {};
for k = 2:numel(lines)
    p = strsplit(strtrim(lines{k}), ',');
    if numel(p) < 3 || isempty(p{2})
        covered{end + 1} = p{1}; %#ok<SAGROW>
        if ~isKey(expected, p{1}), expected(p{1}) = {}; end
    else
        if ~isKey(expected, p{1}), expected(p{1}) = {}; end
        v = expected(p{1}); v{end + 1} = [p{2} '-' p{3}]; expected(p{1}) = v;
    end
end
pfc_report(sprintf('negative series flow: expected data cover all %d cases', numel(cases)), all(ismember(cases, covered)));
warning('off', 'all');
warning('on', 'pfc:rpfc:negativeSeriesFlow');   % must not appear
fprintf('| case | function | n_e | nnz(N) | nnz(A) | issparse N, A, B | condest(A) | attenuation range | acyclic | spectral radius | warnings | negative series flow (count: branches) |\n');
fprintf('|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(cases)
    c = cases{k};
    for fnc = {'apfc', 'rpfc'}
        lastwarn('', '');
        [~, ~, d] = feval(fnc{1}, ['cases/' c], struct('write', false));
        [~, wid] = lastwarn;
        nsf = '';
        if strcmp(fnc{1}, 'rpfc')
            nsf = sprintf('%d: %s', d.negative_series_flow.count, strjoin(d.negative_series_flow.branches, '; '));
            got = regexprep(d.negative_series_flow.branches, '^\d+ \((\d+)-(\d+)\)$', '$1-$2');
            exp_set = expected(c);
            ok = d.negative_series_flow.count == numel(exp_set) && isempty(setxor(got, exp_set)) && ...
                 ~strcmp(wid, 'pfc:rpfc:negativeSeriesFlow') && ~any(strcmp(d.warnings, 'pfc:rpfc:negativeSeriesFlow'));
            pfc_report(sprintf('diag: %s rpfc negative_series_flow equals the old method''s carrier set (%d branches), no warning', ...
                c, numel(exp_set)), ok, sprintf('got {%s}, expected {%s}, last warning ''%s''', strjoin(got, ', '), strjoin(exp_set, ', '), wid));
        end
        fprintf('| %s | %s | %d | %d | %d | %d, %d, %d | %.17g | [%.17g, %.17g] | %d | %.17g | %s | %s |\n', c, fnc{1}, d.n_e, ...
            d.nnz_N, d.nnz_A, d.issparse.N, d.issparse.A, d.issparse.B, d.condest_A, d.attenuation_range(1), ...
            d.attenuation_range(2), d.acyclic, d.spectral_radius, strjoin(d.warnings, ', '), nsf);
        if strcmp(c, 'tap_loop') && strcmp(fnc{1}, 'rpfc')
            ok = abs(d.spectral_radius - 0.17937561102916488) <= 1e-12 && ~d.acyclic && ...
                 isequal(d.negative_series_flow.branches, {'1 (2-1)'});
            tap_line = sprintf('diag: tap_loop rpfc rho = %.17g, acyclic = 0, negative series flow in transformer A (branch 1, 2-1)', d.spectral_radius);
            tap_ok = ok;
        elseif any(strcmp(c, {'case1', 'case9', 'case14', 'case24_ieee_rts', 'case30', 'case33', 'case118'}))
            pfc_report(sprintf('diag: %s %s acyclic with spectral radius 0', c, fnc{1}), d.acyclic && d.spectral_radius == 0, ...
                sprintf('acyclic = %d, rho = %.17g', d.acyclic, d.spectral_radius));
        end
    end
end
pfc_report(tap_line, tap_ok);
