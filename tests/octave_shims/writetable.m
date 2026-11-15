function writetable(T, fname)
% T0 stand-in for Matlab's writetable (scratch use only): writes CSV with %.17g.
fid = fopen(fname, 'w');
fprintf(fid, 'Line'); fprintf(fid, ',%s', T.names{:}); fprintf(fid, '\n');
for i = 1:size(T.data, 1)
  fprintf(fid, '%s', T.rows{i}); fprintf(fid, ',%.17g', full(T.data(i, :))); fprintf(fid, '\n');
end
fclose(fid);
end
