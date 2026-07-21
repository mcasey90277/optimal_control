function [z, prob, info] = ztl_ms_seed(lam0, rv0, rvf, tf, P, M, method)
% ZTL_MS_SEED  Build a DYNAMICALLY CONSISTENT multiple-shooting seed by
% chopping a single-shooting trajectory into arcs.
%
% Integrates the ramp-family flow from [rv0; 1; lam0] SEQUENTIALLY through M
% arcs (node times from ztl_ms_nodes), recording each node state. Because
% every interior node is the integrated endpoint of the previous arc, the
% continuity residual of the seed is ZERO (to integrator tolerance) by
% construction -- only the terminal BC is off. This is exactly the property
% Rung A could never get from mesh-accuracy duals: the seed is on the flow,
% not near it.
%
% INPUTS:
%   lam0   - initial costates (BE convention) [7x1]
%   rv0    - initial position/velocity [1x6]
%   rvf    - target position/velocity [1x6]
%   tf     - final time (ND) [scalar]
%   P      - ztl problem struct (.muStar .c .Tmax .eps + tolerances)
%   M      - number of arcs [scalar]
%   method - node placement 'uniform' | 'amplification' [default 'amplification']
%
% OUTPUTS:
%   z    - MS unknown vector [14M-7 x 1] = [lam0; Y_2; ...; Y_M]
%   prob - problem struct for ztl_ms_residual (.rv0 .rvf .tNodes .M .P)
%   info - struct: .maxContSeed (continuity residual of the seed, should be
%          ~integrator tol) .termErrSeed (terminal BC of the seed)
%
% REFERENCES: ztl_flow.m, ztl_ms_residual.m, ztl_ms_nodes.m; Z0_BUILD.md.

% Default 'uniform': amplification-equidistribution did NOT lower cond(J)
% (perigee amplification is intrinsic, un-splittable; probe_cond_amp) -- kept
% as an option, not the default.
if nargin < 7 || isempty(method), method = 'uniform'; end
tNodes = ztl_ms_nodes(lam0, rv0, tf, P, M, method);
prob = struct('rv0', rv0(:).', 'rvf', rvf(:).', 'tNodes', tNodes, 'M', M, 'P', P);

Y = cell(1, M);
Y{1} = [rv0(:); 1; lam0(:)];
for k = 1:M-1
    o = ztl_flow(Y{k}, [tNodes(k) tNodes(k+1)], P, false);
    Y{k+1} = o.yf;
end

z = zeros(14*M - 7, 1);
z(1:7) = lam0(:);
for k = 2:M
    z(7 + 14*(k-2) + (1:14)) = Y{k};
end

% verify the seed is on the flow (continuity ~0)
[~, ~, ri] = ztl_ms_residual(z, prob, false);
info = struct('maxContSeed', ri.maxCont, 'termErrSeed', ri.termErr);
end
