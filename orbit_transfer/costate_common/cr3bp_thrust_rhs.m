function dz = cr3bp_thrust_rhs(z, u, muStar, Tmax, c)
%% Purpose:
%
%   CR3BP dynamics with thrust and mass flow -- the single shared RHS for
%   flown-control verification (flown_control_error, true_min_altitude).
%   Extracted verbatim from certify_dro_mintime/local_f, which itself
%   mirrors dro_residual/local_rhs. One home (migration #4).
%
%% Inputs:
%
%  z                        [7 x 1]                 State [r; v; m]
%
%  u                        [4 x 1]                 Control [dir(3);
%                                                   throttle]; direction is
%                                                   normalized, throttle
%                                                   clamped to [0, 1]
%
%  muStar                   double                  CR3BP mass ratio
%
%  Tmax                     double                  ND thrust acceleration
%                                                   at unit mass fraction
%
%  c                        double                  ND exhaust velocity
%
%% Outputs:
%
%  dz                       [7 x 1]                 State derivative
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

mu = muStar;
r = z(1:3);  v = z(4:6);  m = z(7);
al = u(1:3);  al = al/max(norm(al),eps);  th = min(max(u(4),0),1);
dd = sqrt((r(1)+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
rr = sqrt((r(1)-1+mu)^2 + r(2)^2 + r(3)^2 + 1e-12);
gr = [r(1) - (1-mu)*(r(1)+mu)/dd^3 - mu*(r(1)-1+mu)/rr^3;
      r(2) - (1-mu)*r(2)/dd^3      - mu*r(2)/rr^3;
           - (1-mu)*r(3)/dd^3      - mu*r(3)/rr^3];
dz = [v; gr + [2*v(2); -2*v(1); 0] + (th*Tmax/m)*al; -(Tmax/c)*th];
end
