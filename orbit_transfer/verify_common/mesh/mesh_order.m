function o = mesh_order(vals, factors)
% MESH_ORDER  Observed order of accuracy and Richardson limit from a mesh ladder.
%
% Given a quantity evaluated at L successive refinement levels (coarse to
% fine) with a CONSTANT refinement ratio r, estimates
%
%   p    = log( (v2 - v1) / (v3 - v2) ) / log(r)
%   rich = v3 + (v3 - v2) / (r^p - 1)
%
% ASSUMPTION, AND WHERE IT BREAKS. Both formulas assume the error admits a
% smooth asymptotic expansion in h with a single leading term. That assumption
% is questionable across a control discontinuity: when switch locations MOVE
% between meshes, the leading coefficient oscillates with mesh phase and a
% three-level extrapolation can be badly misleading. The study's Richardson
% policy therefore forbids applying this to switch counts, raw switch times on
% an unaligned grid, raw endpoint defect multipliers, sign percentages, KKT
% residuals, and anything measured across a branch or topology change. This
% function cannot detect those situations -- the CALLER must respect the
% policy. What it does detect is the arithmetic symptom (differences that flip
% sign or stop shrinking), and it refuses to report an order when it sees one.
%
% REFUSING IS THE POINT. A confidently wrong p is worse than a NaN, because
% the study's conclusions are read directly off these numbers. Sign flip, zero
% last difference, or fewer than three levels all yield NaN rather than a
% plausible-looking value.
%
% INPUTS:
%   vals    - the quantity at each level, coarse to fine [1xL], L >= 2
%   factors - node-count multipliers, coarse to fine [1xL], e.g. [1 2 4 8].
%             The ratio between successive factors must be CONSTANT.
%
% OUTPUTS:
%   o - struct:
%       .p                observed order from the finest three levels, NaN if
%                         L < 3 or the sequence is degenerate [scalar]
%       .rich             Richardson-extrapolated limit, NaN under the same
%                         conditions [scalar]
%       .dLast            |v_L - v_{L-1}| [scalar]
%       .rel              .dLast relative to |v_L| [scalar]
%       .monotone         true if successive differences share a sign and are
%                         shrinking in magnitude [logical]
%       .pWindows         sliding three-level slopes, [1x(L-2)]; empty if L < 3
%       .windowsConsistent  true if all .pWindows agree within 25% of their
%                         mean -- the evidence that an asymptotic regime was
%                         actually reached, which a single slope cannot give
%       .r                the (constant) refinement ratio [scalar]
%
% REFERENCES:
%   [1] Roache, "Verification and Validation in Computational Science and
%       Engineering," 1998 -- observed order and grid-convergence index.
%   [2] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 2
%       and the "Richardson policy (restricted)" section.

vals    = vals(:).';
factors = factors(:).';
L = numel(vals);
assert(L >= 2, 'mesh_order:levels', 'need at least 2 levels, got %d', L);
assert(numel(factors) == L, 'mesh_order:size', ...
    'factors has %d entries, vals has %d', numel(factors), L);

ratios = factors(2:end) ./ factors(1:end-1);
r = ratios(1);
assert(max(abs(ratios - r)) < 1e-12 * max(1, abs(r)), 'mesh_order:ratio', ...
    ['refinement ratio must be CONSTANT across the ladder (got [%s]); ' ...
     'p = log(.)/log(r) has no single r otherwise'], num2str(ratios, '%g '));
assert(r > 1, 'mesh_order:ratio', 'factors must increase coarse-to-fine (r = %g)', r);

d = diff(vals);                       % L-1 successive differences

o = struct('p', NaN, 'rich', NaN, ...
           'dLast', abs(d(end)), 'rel', abs(d(end)) / max(abs(vals(end)), realmin), ...
           'monotone', false, 'pWindows', [], 'windowsConsistent', false, 'r', r);

% Monotone: same sign throughout AND shrinking in magnitude. Either failure
% means the sequence is not settling, whatever a slope formula would say.
if all(d > 0) || all(d < 0)
    o.monotone = all(diff(abs(d)) < 0);   % magnitudes strictly shrinking
end

if L < 3
    return
end

% Sliding three-level slopes. With four levels these are two independent
% estimates whose agreement is the evidence of an asymptotic regime; the plan
% requires four levels precisely so this check exists.
nw = L - 2;
pw = nan(1, nw);
for k = 1:nw
    d1 = d(k);  d2 = d(k+1);
    pw(k) = local_slope(d1, d2, r);
end
o.pWindows = pw;

if all(isfinite(pw))
    mu = mean(pw);
    o.windowsConsistent = abs(mu) > 0 && max(abs(pw - mu)) <= 0.25 * abs(mu);
end

o.p = pw(end);                        % the finest window is the reported order
if isfinite(o.p)
    denom = r^o.p - 1;
    if abs(denom) > eps
        o.rich = vals(end) + d(end) / denom;
    end
end
end

% ---------------------------------------------------------------------------
function p = local_slope(d1, d2, r)
% LOCAL_SLOPE  One three-level observed order, or NaN if the pair is degenerate.
%
% Degenerate means: either difference is zero (the sequence stalled, so the
% ratio is 0 or infinite), or they have opposite signs (the sequence turned
% around, so no single power law fits). Both return NaN by design -- see the
% "refusing is the point" note in the main header.
%
% INPUTS:  d1, d2 - successive differences [scalar]; r - refinement ratio
% OUTPUTS: p - observed order [scalar], NaN when undefined
if d1 == 0 || d2 == 0 || sign(d1) ~= sign(d2) || ~isfinite(d1) || ~isfinite(d2)
    p = NaN;
    return
end
p = log(abs(d1 / d2)) / log(r);
end
