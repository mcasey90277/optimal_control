function tNodes = ztl_ms_nodes(lam0, rv0, tf, P, M, method)
% ZTL_MS_NODES  Multiple-shooting node times, uniform or amplification-
% equidistributing.
%
% The MS conditioning is dominated by a FEW violent early-perigee arcs (at
% M=52 uniform, arc-STM norms span 4 .. 1.85e5, a 12000x ratio -- probe_cond_
% vs_nodes). Uniform-in-time nodes waste themselves in the benign outer-orbit
% coast and leave those perigee arcs under-resolved, so cond(J) plateaus at
% ~1e11 regardless of M. The cure is to place nodes so each arc carries
% roughly EQUAL log-amplification: equidistribute the ACCUMULATED STM growth
% g(t) = log||Phi(t,0)|| (a staircase that rises at each perigee passage).
% NOTE: int||A||dt is the WRONG monitor -- it is dominated by the
% instantaneous 1/r^3 perigee spike and packs nodes into the perigee instant
% while starving the straddling arcs (made cond WORSE, ~1e16). log||Phi(t,0)||
% is the true amplification each arc must resolve.
%
% INPUTS:
%   lam0   - initial costates (BE convention) [7x1]
%   rv0    - initial position/velocity [1x6]
%   tf     - final time (ND) [scalar]
%   P      - ztl problem struct
%   M      - number of arcs [scalar]
%   method - 'uniform' | 'amplification' [default 'amplification']
%
% OUTPUTS:
%   tNodes - node times [1x(M+1)], tNodes(1)=0, tNodes(end)=tf, increasing
%
% REFERENCES: node-equidistribution for MS BVPs (Ascher, Mattheij, Russell,
%   "Numerical Solution of BVPs for ODEs", ch.9); probe_cond_vs_nodes.m.

if nargin < 6 || isempty(method), method = 'amplification'; end

if strcmpi(method, 'uniform')
    tNodes = linspace(0, tf, M+1);
    return
end

% --- amplification clock g(t) = log||Phi(t,0)|| from an STM monitor pass -----
% Integrate [y; Phi(:)] over [0,tf] on a fine grid (eps=1 interior throttle
% is single-regime 'medium'; the monitor is for PLACEMENT, so a fixed regime
% is fine) and record log||Phi(t,0)|| at each grid point.
tGrid = linspace(0, tf, 400);
z0 = [rv0(:); 1; lam0(:); reshape(eye(14), [], 1)];
opts = odeset('RelTol', getdef(P,'odeRelTol',1e-13), 'AbsTol', getdef(P,'odeAbsTol',1e-15));
[~, Z] = ode89(@(tt,zz) monitor_rhs(zz, P), tGrid, z0, opts);

g = zeros(numel(tGrid), 1);
for i = 1:numel(tGrid)
    Phi = reshape(Z(i, 15:end), 14, 14);
    g(i) = log(norm(Phi));
end
g = cummax(g);                                   % monotone amplification clock
if g(end) - g(1) < 1e-6
    tNodes = linspace(0, tf, M+1);  return       % no growth -> uniform
end

gTargets = linspace(g(1), g(end), M+1);
[gu, ig] = unique(g);
tNodes = interp1(gu, tGrid(ig), gTargets, 'linear');
tNodes(1) = 0;  tNodes(end) = tf;               % pin exact endpoints
tNodes = sort(tNodes);                           % guard monotonicity
end

% ---------------------------------------------------------------------------
function dz = monitor_rhs(z, P)
% [y; Phi(:)] flow in the fixed 'medium' regime (placement monitor).
y = z(1:14);  Phi = reshape(z(15:end), 14, 14);
dz = [ztl_eom(y, P, 'medium'); reshape(ztl_A(y, P, 'medium')*Phi, [], 1)];
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
