function [yDot, Ht, Se, u] = lt_pmp_eom_energy(~, y, Tmax, c, muStar)
% LT_PMP_EOM_ENERGY  Min-ENERGY low-thrust CR3BP dynamics under smooth PMP.
%
% Augmented state + costate dynamics for the MINIMUM-ENERGY transfer
% (fixed final time, cost J = int (1/2) u^2 dt with throttle u in [0,1]).
% Because the running cost is QUADRATIC in u, the optimal throttle is a
% CONTINUOUS saturated ramp, not the bang-bang of min-fuel:
%     dH/du = u - S_e = 0  ->  u* = sat(S_e, 0, 1),
%     S_e = Tmax*( ||lambda_v||/m + lambda_m/c ).
% The thrust DIRECTION is unchanged (Lawden's primer,
% alpha = -lambda_v/||lambda_v||). No smoothing homotopy is needed: u is
% already Lipschitz in the costates (only the two saturation corners are
% nonsmooth, a measure-zero set), so the shooting residual is C^0 and the
% convergence basin is far larger than min-fuel's.
%
% The costate ODEs are IDENTICAL to min-time / min-fuel: by the envelope
% theorem they depend only on the APPLIED u (and the primer direction),
% not on the objective. Only the control LAW differs between the three.
%
% INPUTS:
%   ~      - time (unused; autonomous) [scalar]
%   y      - augmented state [14x1], same layout as LT_PMP_EOM
%   Tmax   - max thrust acceleration at m = 1 (ND) [scalar]
%   c      - exhaust velocity (ND) [scalar]
%   muStar - Earth-Moon mass ratio [scalar]
%
% OUTPUTS:
%   yDot - d/dtau of the augmented state [14x1]
%   Ht   - min-energy Hamiltonian (running cost included) [scalar]
%   Se   - unconstrained energy minimizer u* = S_e (pre-saturation) [scalar]
%   u    - saturated throttle in [0,1] [scalar]
%
% REFERENCES:
%   [1] Bertrand & Epenoy, "New smoothing techniques...," OCAM 23(4), 2002.
%   [2] Caillau, Gergaud, Noailles, "3D transfer... energy minimization,"
%       JOTA 2003 (energy-to-fuel continuation).

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

% Min-energy control law (smooth saturation; complex-step safe)
lamvMag = sqrt(sum(lambda_v.^2));
alpha   = -lambda_v./lamvMag;
Se      = Tmax*(lamvMag/m + lambda_m/c);
% u = sat(Se,0,1). Branch on real(Se) so the complex step passes through
% the LINEAR interior (u = Se) and yields exactly zero sensitivity in the
% saturated regions (u const) -- the correct derivative of a clamp.
if real(Se) <= 0
    u = 0*Se;
elseif real(Se) >= 1
    u = 1 + 0*Se;
else
    u = Se;
end

rDot = v;
vDot = gr + hv + u*Tmax/m.*alpha;
mDot = -u*Tmax/c;

lambda_rDot = -G.'*lambda_v;
lambda_vDot = -lambda_r - Hc.'*lambda_v;
lambda_mDot = -lamvMag*u*Tmax/m^2;

yDot = [rDot; vDot; mDot; lambda_rDot; lambda_vDot; lambda_mDot];

if nargout > 1
    Ht = 0.5*u^2 + lambda_r.'*rDot + lambda_v.'*vDot + lambda_m*mDot;
end
end
