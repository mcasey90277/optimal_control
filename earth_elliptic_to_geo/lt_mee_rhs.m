function [dXdL, Ldot] = lt_mee_rhs(X, U, par)
% LT_MEE_RHS  Gauss variational equations for low-thrust 2-body motion in
% Modified Equinoctial Elements (MEE), independent variable L (true longitude).
%
% State X = [P; ex; ey; hx; hy; m; t] evolves per the paper's Gauss-equation
% block (RTN thrust components), converted from d/dt to d/dL via Ldot. Written
% without norm/abs/max/if on state-dependent quantities and with trig taken on
% mod(L,2*pi), so it evaluates on both numeric doubles and CasADi MX symbolics.
%
% INPUTS:
%   X   - State [P; ex; ey; hx; hy; m; t] [7x1]
%   U   - Control [beta(3); thr] with beta a unit RTN thrust direction
%         (radial, transverse, normal) and thr in [0,1] throttle [4x1]
%   par - Struct from kepler_lt_params, with extra field par.L set to the
%         independent-variable value (true longitude) at this node [struct]
%
% OUTPUTS:
%   dXdL - d/dL of state = (dX/dt)/Ldot [7x1]
%   Ldot - dL/dt [scalar]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004, p.6 (Gauss equations).
%   [2] earth_elliptic_to_geo/DESIGN.md sec 2 (problem statement).

P  = X(1);  ex = X(2);  ey = X(3);  hx = X(4);  hy = X(5);  m = X(6);
L  = par.L;  Tm = par.Tmax;  c = par.c;  mu = par.mu;

thr = U(4);
q = thr*U(1);  s = thr*U(2);  w = thr*U(3);

cL = cos(mod(L,2*pi));  sL = sin(mod(L,2*pi));
Z  = 1 + ex*cL + ey*sL;
A1 = ex + (1+Z)*cL;
A2 = ey + (1+Z)*sL;
Xh = 1 + hx^2 + hy^2;
hterm = hx*sL - hy*cL;
sqPmu = sqrt(P/mu);

Pdot  = (2*Tm/m)*sqrt(P^3/mu) * (s/Z);
exdot = (Tm/m)*sqPmu*(1/Z)*( Z*sL*q + A1*s - ey*hterm*w );
eydot = (Tm/m)*sqPmu*(1/Z)*(-Z*cL*q + A2*s + ex*hterm*w );
hxdot = (Tm/(2*m))*sqPmu*(Xh/Z)*cL*w;
hydot = (Tm/(2*m))*sqPmu*(Xh/Z)*sL*w;
Ldot  = sqrt(mu/P^3)*Z^2 + (1/m)*sqPmu*(1/Z)*hterm*w;
mdot  = -(Tm/c)*thr;             % ||(q,s,w)|| = thr since ||beta||=1
tdot  = 1;

dXdt  = [Pdot; exdot; eydot; hxdot; hydot; mdot; tdot];
dXdL  = dXdt / Ldot;
end
