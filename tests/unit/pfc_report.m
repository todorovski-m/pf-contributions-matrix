function pfc_report(name, ok, msg)
%PFC_REPORT  Print one PASS or FAIL line for a unit test.
if nargin < 3, msg = ''; end
if ok
    fprintf('PASS %s\n', name);
else
    fprintf('FAIL %s: %s\n', name, msg);
end
end
