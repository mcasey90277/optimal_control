function m = mesh_match_switches(tA, tB, tol)
% MESH_MATCH_SWITCHES  Pair bang-bang switches across two meshes by physical time.
%
% Comparing switch structure across refinement levels cannot be done by index:
% the meshes have different node counts, and a finer mesh may resolve a short
% arc the coarser one missed entirely, so the switch COUNTS differ too. The
% only stable identity a switch has across meshes is its location in physical
% time, which is what this pairs on.
%
% AN UNMATCHED SWITCH IS A FINDING, NOT A NUISANCE. A switch present at one
% refinement level and absent at another is not a robust feature of the
% continuous solution -- it is evidence that the topology has not stabilized,
% which the study must report rather than absorb. This function therefore
% returns the unmatched sets explicitly and never silently drops or merges.
%
% MATCHING IS GLOBAL-NEAREST, NOT SEQUENTIAL. Repeatedly taking the closest
% surviving pair (and striking out both) is what prevents an early switch from
% claiming a partner that belongs to a later one. A sequential left-to-right
% walk gets this wrong whenever a switch is missing from one side: every
% subsequent pairing is shifted by one and the reported dt values become
% meaningless.
%
% INPUTS:
%   tA  - switch times on mesh A, physical units, monotonic increasing [1xp]
%   tB  - switch times on mesh B, physical units, monotonic increasing [1xq]
%         (by convention A is the coarser mesh and B the finer, which fixes
%          the sign of m.dt below; nothing enforces it)
%   tol - matching window in the same physical units [scalar]. A pair further
%         apart than this is not formed.
%
% OUTPUTS:
%   m - struct:
%       .pairs       matched index pairs, one row per pair [kx2] as [iA, iB],
%                    sorted by iA
%       .unmatchedA  indices into tA with no partner [1xp']
%       .unmatchedB  indices into tB with no partner [1xq']
%       .dt          signed differences tB(iB) - tA(iA) for each pair [1xk].
%                    Signed, not absolute: a systematic drift in one direction
%                    is different evidence from symmetric scatter.
%       .maxAbsDt    max(|dt|), or NaN when no pair formed [scalar]
%       .n           number of pairs [scalar]
%
% REFERENCES:
%   [1] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 1.
%   [2] mesh_switch_times.m -- produces the sub-grid times this consumes.
%       Do NOT feed this function thresholded node indices; see that header.

tA = tA(:).';
tB = tB(:).';
assert(isscalar(tol) && isfinite(tol) && tol > 0, ...
    'mesh_match_switches:tol', 'tol must be a positive finite scalar');

p = numel(tA);  q = numel(tB);
m = struct('pairs', zeros(0,2), 'unmatchedA', [], 'unmatchedB', [], ...
           'dt', [], 'maxAbsDt', NaN, 'n', 0);

if p == 0 || q == 0
    m.unmatchedA = 1:p;
    m.unmatchedB = 1:q;
    if p == 0, m.unmatchedA = zeros(1,0); end
    if q == 0, m.unmatchedB = zeros(1,0); end
    return
end

% Full cost matrix, then repeated global minimum. p and q are switch counts
% (25-360 here), so the dense matrix is trivially small.
D = abs(tA(:) - tB(:).');          % [p x q]
freeA = true(1,p);  freeB = true(1,q);
pairs = zeros(0,2);

while true
    Dm = D;
    Dm(~freeA, :) = inf;
    Dm(:, ~freeB) = inf;
    [best, lin] = min(Dm(:));
    if ~isfinite(best) || best > tol
        break                       % nothing left inside the window
    end
    [ia, ib] = ind2sub(size(Dm), lin);
    pairs(end+1, :) = [ia, ib];     %#ok<AGROW> -- bounded by min(p,q)
    freeA(ia) = false;
    freeB(ib) = false;
end

if ~isempty(pairs)
    pairs = sortrows(pairs, 1);
    m.pairs = pairs;
    m.dt    = tB(pairs(:,2)) - tA(pairs(:,1));
    m.maxAbsDt = max(abs(m.dt));
    m.n = size(pairs,1);
end
m.unmatchedA = find(freeA);
m.unmatchedB = find(freeB);
end
