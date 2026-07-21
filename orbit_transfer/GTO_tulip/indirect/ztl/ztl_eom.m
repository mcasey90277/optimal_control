function [yDot, aux] = ztl_eom(y, P, regime)
% ZTL_EOM  CR3BP min-fuel PMP dynamics, Bertrand-Epenoy RAMP throttle family,
% regime-explicit (complex-step safe within a fixed regime).
%
% Cost family J_eps = (Tmax/c) int [u - eps*u*(1-u)] dt (eps=1 energy,
% eps=0 fuel). Switching function S = 1 - ||lam_v||c/m - lam_m; optimal
% throttle is the clamped ramp u* = clamp(1/2 - S/(2 eps), 0, 1) whose
% clamp state is passed in as REGIME (no branching on y in here -- that is
% what makes complex-step probing of this field exact):
%   'on'     u = 1
%   'medium' u = 1/2 - S/(2 eps)   (eps > 0 only)
%   'off'    u = 0
% Thrust direction is the primer alpha = -lam_v/||lam_v||. Costate equations
% are regime-independent in FORM (envelope theorem); u's y-dependence in
% 'medium' enters the Jacobian via ztl_A's complex step automatically.
%
% INPUTS:
%   y      - augmented state [14x1]: [r(3); v(3); m; lam_r(3); lam_v(3);
%            lam_m]; MAY be complex (CS probing)
%   P      - struct: .muStar .c .Tmax (ND accel at m=1) .eps
%   regime - 'on' | 'medium' | 'off'
%
% OUTPUTS:
%   yDot - d/dt of the augmented state [14x1]
%   aux  - (optional) struct: .S .Sdot .u .Ht  (diagnostic; real inputs)
%
% REFERENCES:
%   [1] Bertrand & Epenoy, OCAM 23(4), 2002 (ramp smoothing family).
%   [2] Zhang et al., JGCD 38(8), 2015 (indirect min-fuel CR3BP).
%   [3] lt_pmp_eom_minfuel.m (source of the gr/hv/G/Hc blocks).
%   [4] Z0_BUILD.md (this increment's spec).

r     = y(1:3);
v     = y(4:6);
m     = y(7);
lam_r = y(8:10);
lam_v = y(11:13);
lam_m = y(14);
muStar = P.muStar;  c = P.c;  Tmax = P.Tmax;

% --- CR3BP field blocks (identical to lt_pmp_eom_minfuel; CS-safe) ---------
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

% --- ramp throttle at fixed regime (NO branching on y) ----------------------
lamvMag = sqrt(sum(lam_v.^2));       % NOT norm() -- CS safety
alpha   = -lam_v./lamvMag;
S       = 1 - lamvMag*c/m - lam_m;
switch regime
    case 'on',     u = 1;
    case 'off',    u = 0;
    case 'medium'
        assert(P.eps > 0, 'ztl_eom: medium regime requires eps > 0');
        u = 0.5 - S/(2*P.eps);
    otherwise
        error('ztl_eom:regime', 'unknown regime %s', regime);
end

rDot  = v;
vDot  = gr + hv + u*(Tmax/m).*alpha;
mDot  = -u*Tmax/c;
lrDot = -G.'*lam_v;
lvDot = -lam_r - Hc.'*lam_v;
lmDot = -lamvMag*u*Tmax/m^2;

yDot = [rDot; vDot; mDot; lrDot; lvDot; lmDot];

if nargout > 1
    Sdot = -(c/m)*(lam_v.'*lvDot)/lamvMag + (c*lamvMag/m^2)*mDot - lmDot;
    Ht   = (Tmax/c)*(u - P.eps*u*(1-u)) + lam_r.'*rDot + lam_v.'*vDot + lam_m*mDot;
    aux  = struct('S', S, 'Sdot', Sdot, 'u', u, 'Ht', Ht);
end
end
