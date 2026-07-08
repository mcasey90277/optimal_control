function [yDot, Ht, S, aThrust] = lt_pmp_eom(~, y, Tmax, c, muStar)
% LT_PMP_EOM  Min-time low-thrust CR3BP state + costate dynamics under PMP.
%
% Augmented equations of motion for the minimum-time low-thrust transfer in
% the Earth-Moon circular restricted three-body problem (rotating,
% nondimensional barycentric frame). The optimal control is evaluated
% pointwise from the costates via Pontryagin's Minimum Principle:
% thrust direction alpha = -lambda_v/||lambda_v|| (Lawden's primer
% direction) and throttle u in {0,1} from the sign of the switching
% function S = -||lambda_v||*c/m - lambda_m (u = 1 when S <= 0).
%
% INPUTS:
%   ~       - time (unused; autonomous system) [scalar]
%   y       - augmented state [14x1]:
%               y(1:3)   r        position (ND, rotating barycentric)
%               y(4:6)   v        velocity (ND)
%               y(7)     m        mass fraction m/m0 (ND)
%               y(8:10)  lambda_r position costates
%               y(11:13) lambda_v velocity costates (primer = -lambda_v)
%               y(14)    lambda_m mass costate
%   Tmax    - max thrust acceleration at m = 1 (ND) [scalar]
%   c       - exhaust velocity Isp*g0 (ND) [scalar]
%   muStar  - Earth-Moon mass ratio m2/(m1+m2) [scalar]
%
% OUTPUTS:
%   yDot    - d/dtau of the augmented state [14x1]
%   Ht      - min-time Hamiltonian H = 1 + lambda'*f [scalar]
%   S       - switching function [scalar]
%   aThrust - thrust acceleration vector u*Tmax/m*alpha [3x1]
%
% REFERENCES:
%   [1] Zhang, Topputo, Bernelli-Zazzera, Zhao, "Low-Thrust Minimum-Fuel
%       Optimization in the Circular Restricted Three-Body Problem,"
%       JGCD 38(8), 2015.
%   [2] pumpkyn.cr3bp.tfMinEoM (reference implementation with analytic STM).

r        = y(1:3);
v        = y(4:6);
m        = y(7);
lambda_r = y(8:10);
lambda_v = y(11:13);
lambda_m = y(14);

% Distances to the primaries (Earth at [-muStar,0,0], Moon at [1-muStar,0,0])
dd = [r(1) + muStar;     r(2); r(3)];   % Earth -> spacecraft
rr = [r(1) - 1 + muStar; r(2); r(3)];   % Moon  -> spacecraft
d3 = sqrt(sum(dd.^2))^3;
r3 = sqrt(sum(rr.^2))^3;

% g(r): centrifugal + gravity; h(v): Coriolis
gr = [r(1); r(2); 0] - (1 - muStar)*dd./d3 - muStar*rr./r3;
hv = [2*v(2); -2*v(1); 0];

% Gravity-gradient + centrifugal Jacobian G = d g / d r (symmetric)
d5 = sqrt(sum(dd.^2))^5;
r5 = sqrt(sum(rr.^2))^5;
G  = diag([1, 1, 0]) ...
     - (1 - muStar)*(eye(3)./d3 - 3*(dd*dd.')./d5) ...
     -      muStar *(eye(3)./r3 - 3*(rr*rr.')./r5);

% Coriolis Jacobian Hc = d h / d v (skew)
Hc = [0 2 0; -2 0 0; 0 0 0];

% PMP control law (min-time): primer direction + bang-bang throttle
lamvMag = sqrt(sum(lambda_v.^2));
alpha   = -lambda_v./lamvMag;
S       = -lamvMag*c/m - lambda_m;
u       = 1;
if real(S) > 0     % explicit real(): S is complex under complex-step calls
    u = 0;
end
aThrust = u*Tmax/m.*alpha;

% State dynamics
rDot = v;
vDot = gr + hv + aThrust;
mDot = -u*Tmax/c;

% Costate dynamics: lambda_dot = -dH/dx
lambda_rDot = -G.'*lambda_v;
lambda_vDot = -lambda_r - Hc.'*lambda_v;
lambda_mDot = -lamvMag*u*Tmax/m^2;

yDot = [rDot; vDot; mDot; lambda_rDot; lambda_vDot; lambda_mDot];

if nargout > 1
    Ht = 1 + lambda_r.'*rDot + lambda_v.'*vDot + lambda_m*mDot;
end
end
