function f = cr3bp_field(x, muStar)
%% Purpose:
%
%   The BALLISTIC circular restricted three-body field in the rotating
%   frame, ND units: xdot = v, vdot = g(r) + h(v), no thrust. Written out
%   rather than calling pumpkyn.cr3bp.eom because that routine's
%   dimension-argument convention is not documented for a single state, and
%   this derivative feeds a continuation tangent where an ambiguity would be
%   silent. Six lines, no dependencies, matched to the g(r) in
%   doc/algorithms_orbit_transfer.tex section 2.1.
%
%% Inputs:
%
%  x                        [6 x 1]                 [r; v]
%  muStar                   double                  mass ratio
%
%% Outputs:
%
%  f                        [6 x 1]                 dx/dt
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

r = x(1:3);  v = x(4:6);  mu = muStar;
d  = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2);
g = [r(1) - (1-mu)*(r(1)+mu)/d^3 - mu*(r(1)-1+mu)/rr^3;
     r(2) - (1-mu)*r(2)/d^3       - mu*r(2)/rr^3;
          - (1-mu)*r(3)/d^3       - mu*r(3)/rr^3];
h = [2*v(2); -2*v(1); 0];
f = [v; g + h];
end
