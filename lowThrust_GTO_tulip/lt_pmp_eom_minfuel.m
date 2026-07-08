function [yDot, Ht, S, u] = lt_pmp_eom_minfuel(~, y, Tmax, c, muStar, epsSmooth)
% LT_PMP_EOM_MINFUEL  Min-fuel low-thrust CR3BP dynamics under smoothed PMP.
%
% Augmented state + costate dynamics for the MINIMUM-FUEL transfer
% (fixed final time, cost J = int (Tmax/c) u dt = propellant fraction).
% The optimal thrust direction is unchanged from min-time (Lawden's primer,
% alpha = -lambda_v/||lambda_v||); the throttle law becomes bang-bang on
% the min-fuel switching function
%     S = 1 - ||lambda_v||*c/m - lambda_m
% (u = 1 for S < 0, coast for S > 0). For numerical work the bang-bang law
% is smoothed with the Bertrand-Epenoy logarithmic barrier, whose interior
% minimizer has the closed form
%     u = 1/(1 + exp(S/eps)) = (1 - tanh(S/(2 eps)))/2,
% implemented in the tanh form (no exp overflow; complex-step safe).
% eps -> 0 recovers bang-bang; drive it down by continuation. Because u is
% the interior argmin of the smoothed Hamiltonian, the envelope theorem
% leaves the costate equations formally identical to the hard-throttle
% ones evaluated at the smoothed u.
%
% INPUTS:
%   ~         - time (unused; autonomous) [scalar]
%   y         - augmented state [14x1], same layout as LT_PMP_EOM
%   Tmax      - max thrust acceleration at m = 1 (ND) [scalar]
%   c         - exhaust velocity (ND) [scalar]
%   muStar    - Earth-Moon mass ratio [scalar]
%   epsSmooth - throttle smoothing parameter (continuation: 1 -> ~1e-3)
%               [scalar]
%
% OUTPUTS:
%   yDot - d/dtau of the augmented state [14x1]
%   Ht   - smoothed min-fuel Hamiltonian (running cost included) [scalar]
%   S    - min-fuel switching function [scalar]
%   u    - smoothed throttle in (0,1) [scalar]
%
% REFERENCES:
%   [1] Bertrand, Epenoy, "New smoothing techniques for solving bang-bang
%       optimal control problems," OCAM 23(4), 2002.
%   [2] Zhang et al., "Low-Thrust Minimum-Fuel Optimization in the
%       Circular Restricted Three-Body Problem," JGCD 38(8), 2015.

r        = y(1:3);
v        = y(4:6);
m        = y(7);
lambda_r = y(8:10);
lambda_v = y(11:13);
lambda_m = y(14);

dd = [r(1) + muStar;     r(2); r(3)];
rr = [r(1) - 1 + muStar; r(2); r(3)];
d3 = sqrt(sum(dd.^2))^3;
r3 = sqrt(sum(rr.^2))^3;

gr = [r(1); r(2); 0] - (1 - muStar)*dd./d3 - muStar*rr./r3;
hv = [2*v(2); -2*v(1); 0];

d5 = sqrt(sum(dd.^2))^5;
r5 = sqrt(sum(rr.^2))^5;
G  = diag([1, 1, 0]) ...
     - (1 - muStar)*(eye(3)./d3 - 3*(dd*dd.')./d5) ...
     -      muStar *(eye(3)./r3 - 3*(rr*rr.')./r5);
Hc = [0 2 0; -2 0 0; 0 0 0];

% Min-fuel control law (smoothed)
lamvMag = sqrt(sum(lambda_v.^2));
alpha   = -lambda_v./lamvMag;
S       = 1 - lamvMag*c/m - lambda_m;
u       = (1 - tanh(S/(2*epsSmooth)))/2;

rDot = v;
vDot = gr + hv + u*Tmax/m.*alpha;
mDot = -u*Tmax/c;

lambda_rDot = -G.'*lambda_v;
lambda_vDot = -lambda_r - Hc.'*lambda_v;
lambda_mDot = -lamvMag*u*Tmax/m^2;

yDot = [rDot; vDot; mDot; lambda_rDot; lambda_vDot; lambda_mDot];

if nargout > 1
    Ht = (Tmax/c)*u + lambda_r.'*rDot + lambda_v.'*vDot + lambda_m*mDot;
end
end
