function varargout = pfc_checks(action, varargin)
%PFC_CHECKS  Options and validity checks shared by APFC and RPFC.
%
%   OPTS = PFC_CHECKS('options', OPTS) fills in the default options and
%   validates them. Fields (all optional):
%     tol      threshold for deletion and classification decisions
%              (default 1e-6, in the units of the MATPOWER case: MW, Mvar)
%     tol_div  guard against division by (near) zero and threshold for
%              negligible self-contribution terms (default 1e-9)
%     write    write results/<case>_apfc.csv or _rpfc.csv (default true;
%              skipped for struct input without a name)
%     solved   CASEDATA already holds a solved power flow, so RUNPF is not
%              called (default false)
%     name     case name used for the output file when CASEDATA is a struct
%              (default '': no file is written)
%
%   The actions are dispatched by name so that each check can be unit tested
%   with crafted inputs (tests/unit/).

switch action
    case 'options'
        varargout{1} = check_options(varargin{:});
    otherwise
        error('pfc:checks:unknownAction', 'pfc_checks: unknown action ''%s''.', action);
end
end

function opts = check_options(opts)
if nargin < 1 || isempty(opts)
    opts = struct();
end
if ~isstruct(opts)
    error('pfc:options:type', 'opts must be a struct.');
end
defaults = struct('tol', 1e-6, 'tol_div', 1e-9, 'write', true, 'solved', false, 'name', '');
known = fieldnames(defaults);
given = fieldnames(opts);
for k = 1:numel(given)
    if ~any(strcmp(given{k}, known))
        error('pfc:options:unknownField', 'Unknown option ''%s''. Known options: %s.', ...
            given{k}, strjoin(known', ', '));
    end
end
for k = 1:numel(known)
    if ~isfield(opts, known{k})
        opts.(known{k}) = defaults.(known{k});
    end
end
for k = {'tol', 'tol_div'}
    v = opts.(k{1});
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && v > 0 && isfinite(v))
        error('pfc:options:value', 'opts.%s must be a positive finite scalar.', k{1});
    end
end
if ~ischar(opts.name)
    error('pfc:options:value', 'opts.name must be a character vector.');
end
opts.write = logical(opts.write);
opts.solved = logical(opts.solved);
end
