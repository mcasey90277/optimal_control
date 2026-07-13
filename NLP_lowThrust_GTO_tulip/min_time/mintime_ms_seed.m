function [z, prob, info] = mintime_ms_seed(lam0, tf, rv0, rvf, Tmax, c, muStar, M)
% MINTIME_MS_SEED  Build a multiple-shooting seed for MINTIME_MS_RESIDUAL by
% forward-propagating the single-shooting min-time arc and chopping it into M
% pieces. Nodes are set to the arc endpoints, so continuity is EXACTLY zero at
% the seed and the only initial residual is the terminal (= the single-shooting
% residual, ~1e-3); the MS solve then drives that to machine precision.
%
% INPUTS:
%   lam0   - single-shooting initial costates [7x1]
%   tf     - single-shooting final time (ND) [scalar]
%   rv0,rvf- endpoints (ND rotating) [1x6]
%   Tmax,c,muStar - dynamics constants
%   M      - number of shooting arcs [scalar, >=2]
%
% OUTPUTS:
%   z    - seed decision vector [ (14(M-1)+8) x 1 ]
%   prob - struct for mintime_ms_residual (.rv0 .rvf .sig .M .Tmax .c .muStar)
%   info - .maxCont (=~0) .termErr .tf

nY = 14;  assert(M >= 2, 'need M>=2 arcs');
sig = linspace(0, 1, M+1)';
Ynode = zeros(nY, M+1);
Ynode(:,1) = [rv0(:); 1; lam0(:)];
for k = 1:M
    dt = (sig(k+1) - sig(k)) * tf;
    [~, y] = pumpkyn.cr3bp.tfMinProp(dt, Ynode(:,k), Tmax, c, muStar);
    Ynode(:,k+1) = y(end, 1:nY).';
end
z = [lam0(:); reshape(Ynode(:, 2:M), [], 1); tf];
prob = struct('rv0', rv0(:).', 'rvf', rvf(:).', 'sig', sig, 'M', M, ...
              'Tmax', Tmax, 'c', c, 'muStar', muStar);
[~, ~, ri] = mintime_ms_residual(z, prob, false);
info = struct('maxCont', ri.maxCont, 'termErr', ri.termErr, 'tf', tf);
end
