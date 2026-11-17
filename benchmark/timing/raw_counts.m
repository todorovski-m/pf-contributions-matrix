% Raw branch counts of the cases in CASES (set before calling) as loadcase sees them (before my_runpf.m), printed as
% CSV lines case,branches_raw. Used by scan.py to count the branches that my_runpf.m deleted before the export.
printf('case,branches_raw\n');
for k = 1:numel(CASES)
  m = loadcase(CASES{k});
  printf('%s,%d\n', CASES{k}, size(m.branch, 1));
end
