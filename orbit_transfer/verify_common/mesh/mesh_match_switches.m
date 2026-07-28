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
%   tol - matching window in the same physical units. Either a scalar, or a
%         PER-SWITCH vector [1xp] giving each tA(i) its own window. The vector
%         form is what the study uses: the physical step of a mesh that is
%         uniform in longitude or in a Sundman variable varies by 36-39x across
%         a transfer (measured, earth 10 N), so a single global window is at
%         once too loose where the mesh is fine in time and too tight where it
%         is coarse. A median-based scalar window was measured to spuriously
%         refuse a genuine pair whose switch had moved 0.382 while the median
%         step was 0.156. A pair further apart than ITS OWN window is not
%         formed.
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
tol = tol(:).';
p = numel(tA);  q = numel(tB);
assert(~isempty(tol) && all(isfinite(tol)) && all(tol > 0), ...
    'mesh_match_switches:tol', 'tol entries must be positive and finite');
assert(isscalar(tol) || numel(tol) == p, 'mesh_match_switches:tol', ...
    'tol must be scalar or have one entry per tA switch (%d), got %d', p, numel(tol));
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

% Mask out every pair that exceeds ITS OWN window up front. Doing this before
% the greedy loop (rather than testing the running minimum against a scalar) is
% what makes a per-switch window correct: with a vector tol, the globally
% smallest remaining distance may be inadmissible while some larger distance is
% still inside its own, larger, window -- so the loop must not stop at the
% first inadmissible minimum.
if isscalar(tol)
    D(D > tol) = inf;
else
    D(D > repmat(tol(:), 1, q)) = inf;
end

freeA = true(1,p);  freeB = true(1,q);
pairs = zeros(0,2);

while true
    Dm = D;
    Dm(~freeA, :) = inf;
    Dm(:, ~freeB) = inf;
    [best, lin] = min(Dm(:));
    if ~isfinite(best)
        break                       % nothing admissible left
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
