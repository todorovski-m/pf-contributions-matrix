function T = array2table(A, varargin)
% T0 stand-in for Matlab's array2table (scratch use only): keeps data and variable names.
T.data = A; T.names = {}; T.rows = {};
for k = 1:2:numel(varargin)
  if strcmpi(varargin{k}, 'VariableNames'), T.names = varargin{k+1}; end
end
end
