casedata = {'case1', 'case9' 'case14', 'case24_ieee_rts', 'case30', 'case33', 'case118','tap_loop'};
for i = 1:numel(casedata)
    disp(casedata{i});
    apfc(['cases/' casedata{i}]);
    rpfc(['cases/' casedata{i}]);
end
