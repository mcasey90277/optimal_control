function [xdot, A, B] = pdg_dynamics(x, T, P)
% PDG_DYNAMICS  3-DOF powered-descent dynamics with analytic Jacobians.
%
%   rdot = v
%   vdot = g + T/m + aD          aD = -(1/2) rho(z) Cd A |v| v / m  (opt-in)
%   mdot = -|T| / (Isp g0)
%
% Complex-step safe: magnitudes via sqrt(sum(.^2) + 1e-12), no norm/abs/
% max. The +1e-12 is a division guard shared with the CasADi twin -- see
% the SMOOTHED MAGNITUDES note in the body.
%
% INPUTS:
%   x - state [r(3); v(3); m] [7x1]  (SI: m, m/s, kg)
%   T - thrust vector [3x1] [N]
%   P - booster_params struct (uses gvec, Isp, g0, drag.*)
%
% OUTPUTS:
%   xdot - state derivative [7x1]
%   A    - d(xdot)/dx [7x7]   (only when requested)
%   B    - d(xdot)/dT [7x3]   (only when requested)
%
% REFERENCES:
%   [1] Acikmese & Ploen, "Convex Programming Approach to Powered Descent
%       Guidance for Mars Landing," JGCD 2007. (vacuum model)

% SMOOTHED MAGNITUDES (external code review, 2026-08-09). |T| and |v| are
% both computed as sqrt(sum(.^2) + magEps), matching the CasADi twin
% pdg_rhs_casadi in solve_pdg_colloc.m, which has always carried the same
% +1e-12. Two reasons, one per magnitude:
%   * Tmag divides B(7,:) = -T'/(Tmag Isp g0) below. At exactly T=0 the
%     raw form is 0/0 -> NaN, which silently destroys the Jacobian for any
%     caller that asks for one at a zero-thrust state.
%   * vmag divides the drag Jacobian's (v v')/vmag term. v=0 is not a
%     hypothetical here: sim_closed_loop's vertical-arrest mode passes
%     through vz=0 by construction, and tvlqr_design calls this function
%     for A and B along the whole trajectory.
% Both smoothed forms have the correct limit (the drag Jacobian tends to
% -(kD/m)*sqrt(magEps)*I ~ 0 as v->0; B(7,:) tends to 0), and -- the
% reason this is done here rather than with max() or a branch -- the
% analytic Jacobians below remain EXACT derivatives of the smoothed xdot,
% so complex-step differentiation still validates them to machine
% precision. No abs/norm/max appears on any complex-step path.
% magEps is absolute in SI: |T| is O(1e5) N and |v| is O(1e2) m/s on this
% problem, so 1e-12 under the square root is ~1e-22 relative -- far below
% double precision, i.e. it changes nothing away from the singular point.
magEps = 1e-12;
v    = x(4:6);   m = x(7);
Tmag = sqrt(sum(T.^2) + magEps);
aT   = T / m;
aD   = zeros(3,1);
if P.drag.on
    rho  = P.drag.rho0 * exp(-x(3) / P.drag.H);
    vmag = sqrt(sum(v.^2) + magEps);
    kD   = 0.5 * rho * P.drag.Cd * P.drag.A;      % so aD = -kD |v| v / m
    aD   = -kD * vmag * v / m;
end
xdot = [v; P.gvec + aT + aD; -Tmag / (P.Isp * P.g0)];

if nargout > 1
    A          = zeros(7,7);
    A(1:3,4:6) = eye(3);
    A(4:6,7)   = -T / m^2;
    if P.drag.on
        % d(aD)/dv = -(kD/m)(|v| I + v v'/|v|); d/dz via drho/dz = -rho/H;
        % d/dm = +kD |v| v / m^2
        A(4:6,4:6) = -(kD/m) * (vmag*eye(3) + (v*v.')/vmag);
        A(4:6,3)   = A(4:6,3) + (kD*vmag/(P.drag.H*m)) * v;
        A(4:6,7)   = A(4:6,7) + kD * vmag * v / m^2;
    end
    B        = zeros(7,3);
    B(4:6,:) = eye(3) / m;
    B(7,:)   = -T.' / (Tmag * P.Isp * P.g0);
end
end
