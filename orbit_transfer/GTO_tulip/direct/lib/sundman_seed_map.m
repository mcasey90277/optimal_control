function [sigma, X0, U0, tauf0] = sundman_seed_map(Xseed, Useed, tf, sgNorm, pSund, muStar, rv0, rvf)
% SUNDMAN_SEED_MAP  Map a time-mesh trajectory into Sundman coordinates
% (NO-RESAMPLE) for warm-starting CASADI_MINFUEL_SUNDMAN.
%
% Given a collocation-feasible solution on a time mesh, change the independent
% variable time -> tau via dt/dtau = kappa = r1^pSund, carrying time as an 8th
% state. Crucially this uses the seed's OWN nodes (sigma = tau/tauf0) with no
% interpolation onto a uniform mesh: downsampling a ~40-rev oscillatory
% trajectory leaves an irreducible ~1e-2 collocation defect that pins IPOPT in
% the restoration phase, whereas the seed's own nodes make the only initial
% infeasibility the small time-trap vs Sundman-trap mismatch, which the solver
% closes to ~1e-14. Endpoints are pinned exactly.
%
% INPUTS:
%   Xseed  - seed states [7xM] ([r;v;m]) on a time mesh
%   Useed  - seed controls [4xM]; rows 1:3 a thrust direction (unit OR cone
%            [w] with |w|=s), row 4 the throttle s in [0,1]
%   tf     - transfer time (ND) [scalar]
%   sgNorm - seed node parameter [Mx1], AFFINE in physical seed time (any
%            affine range; normalized inside). An arbitrary monotone
%            parameter is not enough: tSeed = sgNorm*tf below assumes it.
%   pSund  - Sundman power [scalar]
%   muStar - CR3BP mass ratio [scalar]
%   rv0    - departure state to pin at node 1 [1x6]
%   rvf    - arrival state to pin at node end [1x6]
%
% OUTPUTS:
%   sigma - Sundman nodes tau/tauf0 [ (M') x1 ], 0->1 (M' <= M after dedup)
%   X0    - warm-start states [8x M'] ([r;v;m;t])
%   U0    - warm-start controls [4x M'] ([alpha;s], ||alpha||=1)
%   tauf0 - total regularized length [scalar]
%
% REFERENCES:
%   [1] casadi_minfuel_sundman.m (the consumer; its carried-time trapezoid is
%       inverted exactly here so the seed carries no time defect).
%   [2] ../../process/LOW_THRUST_MINFUEL_CAMPAIGN.md, "no-resample seed".

sgNorm = sgNorm(:);  sgNorm = (sgNorm - sgNorm(1))/(sgNorm(end) - sgNorm(1));
tSeed  = sgNorm * tf;                                  % physical time per node

w = Useed(1:3,:);  s = Useed(4,:);
% unit direction. A cone seed has w = 0 on coast nodes; the old guard returned
% alpha = 0 there, which is NOT a unit vector and makes the gradient of the
% ||alpha||=1 constraint vanish at the warm start (review 2026-09-06, E7). Use
% the nearest valid direction along the mesh, else a fixed unit fallback.
nw = sqrt(sum(w.^2, 1));  good = nw > 1e-12;
alpha = zeros(size(w));
alpha(:, good) = w(:, good) ./ nw(good);
if any(~good)
    if any(good)
        gi = find(good);
        for kb = find(~good)
            [~, jn] = min(abs(gi - kb));  alpha(:, kb) = alpha(:, gi(jn));
        end
    else
        alpha(:, ~good) = repmat([1; 0; 0], 1, nnz(~good));
    end
end

r1  = sqrt((Xseed(1,:)+muStar).^2 + Xseed(2,:).^2 + Xseed(3,:).^2).';
kap = r1.^pSund;
dt  = diff(tSeed);
% Invert the solver's OWN carried-time trapezoid so the seed carries no
% systematic initial time-defect. The solver enforces (casadi_minfuel_sundman)
% Dt = (dtau/2)*(kappa_L + kappa_R) on the t-state, so the exact discrete
% inverse is dtau = 2*dt/(kappa_L + kappa_R) -- the reciprocal of the mean,
% NOT the mean of reciprocals 0.5*(1/kL+1/kR) (they agree only when kL=kR, so
% the old form biased fast-varying perigee intervals).
dtau = 2*dt./(kap(1:end-1) + kap(2:end));              % discrete inverse of the t-defect trapezoid
tau  = [0; cumsum(dtau)];  tauf0 = tau(end);  sig = tau/tauf0;

[sigma, ku] = unique(sig, 'stable');
X0 = [Xseed(:,ku); tSeed(ku).'];
U0 = [alpha(:,ku); s(ku)];
X0(1:6,1)   = rv0(:);  X0(7,1)   = 1;  X0(8,1)   = 0;
X0(1:6,end) = rvf(:);  X0(8,end) = tf;
end
