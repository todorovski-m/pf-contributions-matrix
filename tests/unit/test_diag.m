% diag output on the shared cases, tap_loop and the synthetic cases (brief T12b, acceptance 2).
% Run with: zsh tests/octave_in_scratch.sh test_diag
warning('off', 'all');
cases = {'case1', 'case9', 'case14', 'case30', 'case33', 'case118', 'tap_loop', ...
         'syn9_neg_pd', 'syn9_neg_pg', 'syn9_neg_qd', 'syn9_neg_qg'};
fprintf('| case | function | n_e | nnz(N) | nnz(A) | issparse N, A, B | condest(A) | attenuation range | acyclic | spectral radius | warnings | negative series flow branches |\n');
fprintf('|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(cases)
    c = cases{k};
    for fn = {'apfc', 'rpfc'}
        [~, ~, d] = feval(fn{1}, ['cases/' c], struct('write', false));
        nsf = '';
        if isfield(d, 'negative_series_flow_branches'), nsf = strjoin(d.negative_series_flow_branches, '; '); end
        fprintf('| %s | %s | %d | %d | %d | %d, %d, %d | %.17g | [%.17g, %.17g] | %d | %.17g | %s | %s |\n', c, fn{1}, d.n_e, ...
            d.nnz_N, d.nnz_A, d.issparse.N, d.issparse.A, d.issparse.B, d.condest_A, d.attenuation_range(1), ...
            d.attenuation_range(2), d.acyclic, d.spectral_radius, strjoin(d.warnings, ', '), nsf);
        if strcmp(c, 'tap_loop') && strcmp(fn{1}, 'rpfc')
            ok = abs(d.spectral_radius - 0.17937561102916488) <= 1e-12 && ~d.acyclic && ...
                 any(strcmp(d.warnings, 'pfc:rpfc:negativeSeriesFlow')) && isequal(d.negative_series_flow_branches, {'1 (2-1)'});
            tl = sprintf('PASS diag: tap_loop rpfc rho = %.17g, acyclic = 0, negative series flow in transformer A (branch 1, 2-1)', d.spectral_radius);
            if ~ok, tl = ['FAIL' tl(5:end)]; end
            tap_line = tl;
        elseif any(strcmp(c, {'case1', 'case9', 'case14', 'case30', 'case33', 'case118'}))
            if ~(d.acyclic && d.spectral_radius == 0)
                fprintf('FAIL diag: %s %s acyclic = %d, rho = %.17g\n', c, fn{1}, d.acyclic, d.spectral_radius);
            end
        end
    end
end
disp(tap_line);
