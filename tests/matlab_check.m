function matlab_check(outdir, refdir)
%MATLAB_CHECK  Run apfc and rpfc in Matlab and compare with the Octave results.
%
%   MATLAB_CHECK(OUTDIR) runs APFC and RPFC on case9, case118 and
%   case24_ieee_rts, writes <case>_apfc.csv and <case>_rpfc.csv with
%   Matlab's writetable to the folder OUTDIR (created if needed), and prints,
%   for each file, whether the row labels and column names agree with the
%   Octave results and the maximum absolute and relative difference of the
%   values. MATLAB_CHECK(OUTDIR, REFDIR) uses the Octave results in REFDIR
%   (default: tests/octave_results, see the README there).
%
%   Run it from the root of the repository, with MATPOWER on the Matlab
%   path and WITHOUT tests/octave_shims on the path (the built-in table
%   functions must be used). Expected: identical labels; value differences
%   of about 1e-12 or less (Matlab and Octave use different linear-algebra
%   libraries, and writetable writes 15 significant digits while the Octave
%   results have 17). It also prints issparse for N, A and B (expected 1).
%
%   This script has not been run: Matlab was not available when it was
%   written. It uses only functions available in Matlab R2018b.

if nargin < 1 || isempty(outdir)
    outdir = fullfile(tempdir, 'pfc_matlab_check');
end
if nargin < 2 || isempty(refdir)
    refdir = fullfile('tests', 'octave_results');
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
cases = {'case9', 'case118', 'case24_ieee_rts'};
fprintf('%-28s %-8s %-8s %-24s %-24s\n', 'file', 'labels', 'names', 'max abs diff', 'max rel diff');
for k = 1:numel(cases)
    c = cases{k};
    [ta, ~, da] = apfc(fullfile('cases', c), struct('write', false));
    [tr, ~, dr] = rpfc(fullfile('cases', c), struct('write', false));
    writetable(ta, fullfile(outdir, [c '_apfc.csv']));
    writetable(tr, fullfile(outdir, [c '_rpfc.csv']));
    fprintf('%s: issparse apfc N, A, B = %d %d %d; rpfc N, A, B = %d %d %d\n', c, ...
        da.issparse.N, da.issparse.A, da.issparse.B, dr.issparse.N, dr.issparse.A, dr.issparse.B);
    for kind = {'apfc', 'rpfc'}
        f = [c '_' kind{1} '.csv'];
        [hm, rm, vm] = read_csv(fullfile(outdir, f));
        [ho, ro, vo] = read_csv(fullfile(refdir, f));
        same_rows = isequal(rm, ro);
        same_names = isequal(hm, ho);
        if isequal(size(vm), size(vo))
            d = abs(vm - vo);
            big = abs(vo) >= 1e-9;
            mx = max(d(:));
            mr = max(d(big) ./ abs(vo(big)));
            if isempty(mr), mr = 0; end
            fprintf('%-28s %-8d %-8d %-24.17g %-24.17g\n', f, same_rows, same_names, mx, mr);
        else
            fprintf('%-28s %-8d %-8d size differs: %dx%d (Matlab) vs %dx%d (Octave)\n', f, same_rows, ...
                same_names, size(vm, 1), size(vm, 2), size(vo, 1), size(vo, 2));
        end
    end
end
fprintf('Matlab outputs written to %s\n', outdir);
end

function [header, rows, values] = read_csv(fn)
%READ_CSV  Read a CSV file with a header line and a first text column.
lines = strsplit(strtrim(fileread(fn)), sprintf('\n'));
lines = strtrim(lines);
header = strsplit(lines{1}, ',');
n = numel(lines) - 1;
rows = cell(n, 1);
values = zeros(n, numel(header) - 1);
for i = 1:n
    p = strsplit(lines{i + 1}, ',');
    rows{i} = p{1};
    values(i, :) = str2double(p(2:end));
end
end
