function s = mesh_switch_times(sg, tNode, burn, Shat, thr)
% MESH_SWITCH_TIMES  Locate bang-bang switches to SUB-GRID accuracy, in physical time.
%
% WHY THIS EXISTS. A mesh-convergence study cannot measure the order of switch
% times using a switch location read off the grid. The repo's existing detector
% (foc_check.m, rep.nSwitches) locates a switch as the node index where the
% thresholded throttle changes, i.e. swI = find(diff(burn) ~= 0). That index is
% quantized to the mesh by construction, so an order fitted to it returns p ~ 1
% for ANY solution -- it measures the grid, not the trajectory. This function
% resolves the crossing inside its bracketing interval instead.
%
% PREFERRED ESTIMATOR: the zero of the switching function. At a regular switch
% the switching function crosses zero transversally (Sdot ~= 0), so it is
% smooth through the crossing and linear interpolation of it is a genuine
% sub-grid locator. The throttle is a step there and carries no such
% information, which is why the throttle crossing is only the FALLBACK -- and
% why the estimator actually used is always reported in s.method, so a report
% can never silently mix the two.
%
% The de-weighted switching function (foc_check's rep.Sdeweighted) is the right
% input: the raw dual-derived S carries the mesh's nodal quadrature weights,
% which vary by orders of magnitude on a Sundman mesh and would bias the
% interpolation toward the wider side of the bracket.
%
% INPUTS:
%   sg    - independent-variable grid, monotonic increasing [N1x1 or 1xN1]
%   tNode - physical time at each node [N1x1 or 1xN1]. Pass sg itself when the
%           independent variable IS physical time. This is the carried time
%           state for Sundman/longitude-domain transcriptions.
%   burn  - throttle state at each node, logical or 0/1 [1xN1]
%   Shat  - de-weighted switching function at each node [1xN1], or [] to force
%           the throttle fallback
%   thr   - (optional) throttle in [0,1] at each node [1xN1]; required only for
%           the fallback path, where the crossing of thr = 0.5 is located
%
% OUTPUTS:
%   s - struct:
%       .tSw         switch times in the units of tNode [1xk]
%       .sgSw        switch locations in the sg domain [1xk]
%       .n           number of switches [scalar]
%       .method      'switchfn' | 'throttle' | 'mixed' [char]
%       .bracket     bracketing node indices, one row per switch [kx2]
%       .subGridFrac fractional position of the crossing within its bracket,
%                    in [0,1] [1xk]. THE DIAGNOSTIC THAT MATTERS: values piled
%                    at 0 or 1 mean the estimator has collapsed back onto
%                    nodes and the extraction is quantized after all.
%       .fellBack    per-switch logical [1xk], true where the switching
%                    function was unusable and the throttle was used instead
%
% REFERENCES:
%   [1] Agamawi, Hager & Rao, "Mesh Refinement Method for Solving Bang-Bang
%       Optimal Control Problems Using Direct Collocation," arXiv:1905.11895 --
%       switch location is a dedicated event-location problem, not a quantity
%       to be read off nodes.
%   [2] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Global
%       Constraint on switch extraction, and Task 1a.

if nargin < 5, thr = []; end
sg    = sg(:).';
tNode = tNode(:).';
burn  = logical(burn(:).');
N1    = numel(sg);
assert(numel(tNode) == N1, 'mesh_switch_times:size', ...
    'tNode has %d entries, sg has %d', numel(tNode), N1);
assert(numel(burn) == N1, 'mesh_switch_times:size', ...
    'burn has %d entries, sg has %d', numel(burn), N1);
if ~isempty(Shat), Shat = Shat(:).'; end
if ~isempty(thr),  thr  = thr(:).';  end

s = struct('tSw', [], 'sgSw', [], 'n', 0, 'method', 'switchfn', ...
           'bracket', zeros(0,2), 'subGridFrac', [], 'fellBack', []);

% Bracket every switch by the throttle transition. The bracket is where the
% switch IS; the estimators below only decide WHERE INSIDE it.
kSw = find(diff(burn) ~= 0);
s.n = numel(kSw);
if s.n == 0
    if isempty(Shat), s.method = 'throttle'; end
    return
end

sgSw     = nan(1, s.n);
frac     = nan(1, s.n);
fellBack = false(1, s.n);

for q = 1:s.n
    ka = kSw(q);  kb = ka + 1;
    useSwitchFn = false;
    if ~isempty(Shat)
        Sa = Shat(ka);  Sb = Shat(kb);
        % A usable bracket needs a genuine sign change with a nonzero spread;
        % without one, linear interpolation is meaningless (or divides by ~0).
        useSwitchFn = isfinite(Sa) && isfinite(Sb) && (Sa ~= Sb) && ...
                      (sign(Sa) ~= sign(Sb)) && (Sa ~= 0 || Sb ~= 0);
    end

    if useSwitchFn
        f = Sa / (Sa - Sb);                 % zero of the linear interpolant
    else
        % Fallback: the throttle's own 0.5 crossing. Guarded rather than
        % silently returning a node index -- a quantized answer is exactly the
        % failure mode this function exists to prevent.
        fellBack(q) = true;
        if isempty(thr)
            error('mesh_switch_times:noEstimator', ...
                ['switch %d has no usable switching-function sign change and no ' ...
                 'throttle was supplied; cannot locate it sub-grid'], q);
        end
        ta = thr(ka);  tb = thr(kb);
        assert(ta ~= tb, 'mesh_switch_times:flatThrottle', ...
            'switch %d brackets a flat throttle (%.3g); cannot locate', q, ta);
        f = (ta - 0.5) / (ta - tb);
    end

    f       = min(max(f, 0), 1);            % keep it inside its own bracket
    frac(q) = f;
    sgSw(q) = sg(ka) + f * (sg(kb) - sg(ka));
end

% Map to physical time through the SAME bracket and fraction. Interpolating
% tNode rather than assuming t = sg is what makes this correct on Sundman and
% longitude-domain meshes, where the physical step varies by orders of
% magnitude across the transfer.
tSw = nan(1, s.n);
for q = 1:s.n
    ka = kSw(q);  kb = ka + 1;
    tSw(q) = tNode(ka) + frac(q) * (tNode(kb) - tNode(ka));
end

s.tSw         = tSw;
s.sgSw        = sgSw;
s.bracket     = [kSw(:), kSw(:) + 1];
s.subGridFrac = frac;
s.fellBack    = fellBack;
if all(fellBack)
    s.method = 'throttle';
elseif any(fellBack)
    s.method = 'mixed';
else
    s.method = 'switchfn';
end
end
