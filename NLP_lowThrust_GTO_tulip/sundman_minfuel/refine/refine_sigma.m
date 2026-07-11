function [sigmaNew, isNew, nDropped] = refine_sigma(sigma, score, opts)
% REFINE_SIGMA  Bisect the top-scoring collocation intervals, preserving nodes.
%
% Inserts the midpoint of each selected sigma-interval, keeping every
% original node exactly (the no-resample discipline). Guards: never bisect
% an interval narrower than hFloor; cap inserted nodes per call at maxAdd.
% Dropped selections are counted (never silently truncated).
%
% INPUTS:
%   sigma - current normalized nodes, 0->1 [(N+1)x1]
%   score - per-interval refinement score [1xN], >= 0
%   opts  - struct: K max intervals to bisect [default min(8,nnz(score>0))],
%           hFloor min sigma-interval width to bisect [default 1e-9],
%           maxAdd cap on inserted nodes [default 40]
%
% OUTPUTS:
%   sigmaNew - refined nodes, sorted, originals preserved [(N'+1)x1]
%   isNew    - logical mask of inserted nodes [(N'+1)x1]
%   nDropped - viable intervals (score>0, width>=hFloor) excluded by a guard:
%              hFloor-thin, or beyond the maxAdd cap when maxAdd<=K [scalar]
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md

sigma = sigma(:);   N = numel(sigma) - 1;   score = score(:).';
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'K'),      opts.K = min(8, nnz(score > 0)); end
if ~isfield(opts, 'hFloor'), opts.hFloor = 1e-9;             end
if ~isfield(opts, 'maxAdd'), opts.maxAdd = 40;               end

h = diff(sigma).';                             % [1xN] interval widths
[~, ord] = sort(score, 'descend');
sel = [];  nDropped = 0;
cap = min(opts.K, opts.maxAdd);
for q = 1:numel(ord)
    k = ord(q);
    if score(k) <= 0, break; end               % no more scored intervals
    if h(k) < opts.hFloor, nDropped = nDropped + 1; continue; end   % guard: too thin
    if numel(sel) >= cap
        % viable interval excluded by a cap. If maxAdd (a safety guard) is the
        % binding limit, count it as dropped; if the K selection quota binds
        % (K < maxAdd), the interval is simply not-selected by design.
        if opts.maxAdd <= opts.K, nDropped = nDropped + 1; end
        continue;                              % keep scanning to count remaining viable
    end
    sel(end+1) = k; %#ok<AGROW>
end
add = false(1, N);  add(sel) = true;

sigmaNew = zeros(N + 1 + numel(sel), 1);  isNew = false(size(sigmaNew));
w = 1;
for k = 1:N
    sigmaNew(w) = sigma(k);  w = w + 1;
    if add(k)
        sigmaNew(w) = 0.5*(sigma(k) + sigma(k+1));  isNew(w) = true;  w = w + 1;
    end
end
sigmaNew(w) = sigma(N + 1);
end
