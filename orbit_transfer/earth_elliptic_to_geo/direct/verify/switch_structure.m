function s = switch_structure(X, U, dL, sigma)
% SWITCH_STRUCTURE  Primal (mesh-dependent) bang-bang switch metrics for a
% converged MEE/sigma-domain low-thrust solution. PURELY PRIMAL -- reads only
% the throttle row and the node longitudes, no costates -- so it is immune to
% the raw-dual/primer anomaly (process/DESIGN_dual_map.md) that corrupts the
% PMP switching-function reconstruction at high eccentricity. This is the
% quantity compared across mesh densities in the P0 switch-count convergence
% certification (process/P0_SWITCH_MESH_CONVERGENCE.md).
%
% A "switch" is a crossing of the throttle through 0.5 between adjacent nodes;
% its true-longitude location is found by linear interpolation of thr across
% the interval (near-exact for a bang-bang throttle). NOTE the count is
% intrinsically mesh-dependent: two switches inside one interval are merged,
% so at coarse density it is a LOWER BOUND (see the P0 note -- 8 nodes/rev
% under-counts the 0.2 N rung by ~5%).
%
% INPUTS:
%   X     - state [P;ex;ey;hx;hy;m;t] [7x(N+1)]
%   U     - control [beta(3);thr] [4x(N+1)]
%   dL    - total true-longitude span (DeltaL) [scalar]
%   sigma - node grid, [0,1] [(N+1)x1] or [1x(N+1)]
%
% OUTPUTS:
%   s - struct:
%       .nSwitch      - count of thr 0.5-crossings [scalar]
%       .revs         - dL/(2*pi) [scalar]
%       .duty         - sigma-weighted burn fraction (trapezoid of thr>0.5) [scalar]
%       .swL          - unwrapped true longitude at each switch [1 x nSwitch]
%       .swPhase      - orbital phase L mod 2*pi at each switch [deg, 1 x nSwitch]
%       .swT          - physical time (ND) at each switch [1 x nSwitch]
%       .nodesPerRev  - (N)/revs [scalar]
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/process/P0_SWITCH_MESH_CONVERGENCE.md (the
%       convergence certification this metric feeds).
%   [2] earth_elliptic_to_geo/verify/meshstudy_switch.m (driver that sweeps
%       node density and compares this metric).
sigma = sigma(:).';
thr   = U(4,:);
burn  = thr > 0.5;
N1    = numel(sigma);
L     = pi + sigma*dL;                     % per-node true longitude (x0's L = pi)

% switch locations: interpolate where thr crosses 0.5 across each transition
% interval, in both true longitude (L) and physical time (X(7,:))
t   = X(7,:);
swI = find(diff(double(burn)) ~= 0);       % interval index of each transition
swL = zeros(1, numel(swI));
swT = zeros(1, numel(swI));
for q = 1:numel(swI)
    k = swI(q);
    denom = thr(k+1) - thr(k);
    if abs(denom) < 1e-12
        frac = 0.5;
    else
        frac = min(max((0.5 - thr(k)) / denom, 0), 1);
    end
    swL(q) = L(k) + frac*(L(k+1) - L(k));
    swT(q) = t(k) + frac*(t(k+1) - t(k));
end

% sigma-weighted duty cycle (trapezoid of the burn indicator)
dsig = diff(sigma);
bmid = 0.5*(burn(1:end-1) + burn(2:end));
s.duty        = sum(dsig .* bmid) / sum(dsig);
s.nSwitch     = numel(swI);
s.revs        = dL/(2*pi);
s.swL         = swL;
s.swPhase     = mod(swL, 2*pi) * 180/pi;
s.swT         = swT;
s.nodesPerRev = (N1-1) / s.revs;
end
