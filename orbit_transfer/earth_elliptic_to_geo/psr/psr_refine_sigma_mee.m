function [sigmaNew, isNew, nDropped] = psr_refine_sigma_mee(sigma, score, opts)
% PSR_REFINE_SIGMA_MEE  Bisect the top-scoring collocation intervals,
% preserving every original node -- PORTED from
% NLP_lowThrust_GTO_tulip/PSR/lib/refine_sigma.m (identical mechanics: insert
% the midpoint of each selected interval, never bisect an interval narrower
% than hFloor, cap inserted nodes per call at maxAdd, count dropped
% selections rather than silently truncating). K/maxAdd DEFAULTS differ from
% the tulip original -- see psr_mee_refine.m's header for why (MEE's ~171
% switches vs the tulip campaign's ~25-40 make a top-8-per-round cap
% incompatible with a 4-round stabilization budget; default here is
% "refine every positively-scored interval", capped only by a generous
% maxAdd safety valve).
%
% Factored out as its own pure function (mirrors interp_warmstart.m's
% rationale) so the mesh-insertion logic is unit-testable without a solve.
%
% INPUTS:
%   sigma - current normalized nodes, 0->1 [(N+1)x1]
%   score - per-interval refinement score [1xN], >= 0
%   opts  - struct: .K max intervals to bisect [default Inf -> all
%           positively-scored intervals], .hFloor min sigma-interval width
%           to bisect [default 1e-9], .maxAdd cap on inserted nodes
%           [default 2000]
%
% OUTPUTS:
%   sigmaNew - refined nodes, sorted, originals preserved [(N'+1)x1]
%   isNew    - logical mask of inserted nodes [(N'+1)x1]
%   nDropped - viable intervals (score>0, width>=hFloor) excluded by a
%              guard: hFloor-thin, or beyond the maxAdd cap when
%              maxAdd<=K [scalar]
%
% REFERENCES:
%   [1] NLP_lowThrust_GTO_tulip/PSR/lib/refine_sigma.m (read-only source of
%       this port, mechanics unchanged).
%   [2] psr_mee_refine.m (K/maxAdd default-deviation rationale).
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'K') || isempty(opts.K),           opts.K = Inf;   end
if ~isfield(opts, 'hFloor') || isempty(opts.hFloor), opts.hFloor = 1e-9; end
if ~isfield(opts, 'maxAdd') || isempty(opts.maxAdd), opts.maxAdd = 2000;  end

sigma = sigma(:);   N = numel(sigma) - 1;   score = score(:).';
if isinf(opts.K), K = nnz(score > 0); else, K = opts.K; end

h = diff(sigma).';                             % [1xN] interval widths
[~, ord] = sort(score, 'descend');
sel = [];  nDropped = 0;
cap = min(K, opts.maxAdd);
for q = 1:numel(ord)
    k = ord(q);
    if score(k) <= 0, break; end               % no more scored intervals
    if h(k) < opts.hFloor, nDropped = nDropped + 1; continue; end   % guard: too thin
    if numel(sel) >= cap
        % viable interval excluded by a cap. If maxAdd (a safety guard) is
        % the binding limit, count it as dropped; if the K selection quota
        % binds (K < maxAdd), the interval is simply not-selected by design.
        if opts.maxAdd <= K, nDropped = nDropped + 1; end
        continue;                              % keep scanning to count remaining viable
    end
    sel(end + 1) = k; %#ok<AGROW>
end
add = false(1, N);  add(sel) = true;

sigmaNew = zeros(N + 1 + numel(sel), 1);  isNew = false(size(sigmaNew));
w = 1;
for k = 1:N
    sigmaNew(w) = sigma(k);  w = w + 1;
    if add(k)
        sigmaNew(w) = 0.5 * (sigma(k) + sigma(k + 1));  isNew(w) = true;  w = w + 1;
    end
end
sigmaNew(w) = sigma(N + 1);
end
