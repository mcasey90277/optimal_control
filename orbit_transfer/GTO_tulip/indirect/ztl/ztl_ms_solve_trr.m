function [z, out] = ztl_ms_solve_trr(z0, prob, opts)
% ZTL_MS_SOLVE_TRR  Solve the multiple-shooting system with lsqnonlin's
% trust-region-reflective method and the EXACT block Jacobian.
%
% The MS Jacobian is ill-conditioned (~1e11, un-improvable by node placement
% -- perigee amplification is intrinsic) and the residual valley is curved,
% so a crude diagonal-LM crawls. lsqnonlin TRR is a mature trust-region
% least-squares (2-D subspace minimization, sparse-aware) that handles both
% -- given the exact Jacobian, which ztl_ms_residual supplies. This is the
% same descent SS could not achieve; the question Z1 answers is whether the
% MS structure lets it reach 1e-8.
%
% INPUTS:
%   z0   - initial unknown vector [14M-7 x 1]
%   prob - problem struct for ztl_ms_residual
%   opts - (optional): .tolR [1e-10] .maxIter [400] .verbose [true]
%
% OUTPUTS:
%   z   - solution;  out: .resNorm .iters .funcCount .flag (lsqnonlin
%         exitflag) .firstOrd .termErr .maxCont
%
% REFERENCES: Coleman & Li, SIAM J. Optim. 6(2), 1996 (TRR); lsqnonlin.

if nargin < 3, opts = struct(); end
g = @(f,d) getdef(opts, f, d);
tolR = g('tolR', 1e-10);  maxIter = g('maxIter', 400);
if g('verbose', true), disp = 'iter'; else, disp = 'off'; end

lo = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'Display', disp, ...
    'FunctionTolerance', tolR^2, ...
    'OptimalityTolerance', 1e-14, ...
    'StepTolerance', 1e-15, ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', 20*maxIter);

[z, resn2, R, flag, output] = lsqnonlin(@(zz) resjac(zz, prob), z0(:), [], [], lo);
[~, ~, ri] = ztl_ms_residual(z, prob, false);

out = struct('resNorm', sqrt(resn2), 'iters', output.iterations, ...
    'funcCount', output.funcCount, 'flag', flag, ...
    'firstOrd', output.firstorderopt, 'termErr', ri.termErr, ...
    'maxCont', ri.maxCont, 'grazed', ri.grazed);
if norm(R) == 0, end  %#ok<*NASGU>
end

% ---------------------------------------------------------------------------
function [R, J] = resjac(z, prob)
if nargout > 1
    [R, J] = ztl_ms_residual(z, prob, true);
else
    R = ztl_ms_residual(z, prob, false);
end
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
