function [F, A, aux] = cr3bp_minenergy_pmp(y, Tmax, c, muStar)
%% Purpose:
%
%   The MINIMUM-ENERGY PMP vector field for the CR3BP low-thrust transfer,
%   with its exact 14x14 Jacobian -- the fixed-final-time sibling of
%   pumpkyn.cr3bp.tfMinEoM. Same 14-state PMP y = [r; v; m; lam_r; lam_v;
%   lam_m], same dynamics and costate ODEs, but the running cost is the
%   Bertrand-Epenoy energy endpoint L = s^2 (throttle squared, physical-time
%   measure -- the SAME convention GTO_tulip's energy->fuel homotopy uses at
%   eps = 1), so the throttle is smooth and interior:
%
%       H  = s^2 + lam_r'v + lam_v'(g + h + s T alpha/m) - lam_m s T/c
%       alpha* = -lam_v/|lam_v|                       (primer direction)
%       s*     = clip( (T/2)(|lam_v|/m + lam_m/c), 0, 1 )
%
%   At saturation (s* = 1) F and A coincide with tfMinEoM's (u = 1) to
%   round-off -- the unit test pins every convention to the reference. The
%   Jacobian is built ONCE by CasADi SX automatic differentiation (persistent
%   Function), so no hand-derived third derivative of the potential can be
%   wrong; the clip is a CasADi fmin/fmax and differentiates one-sided.
%
%  ASSUMPTIONS / NOTES:
%
% • pumpkyn has no min-energy EoM, so this is the one costate-pipeline
%   dynamics that cannot be a pumpkyn call. Nothing is written into pumpkyn.
% • Requires CasADi on the path (self-resolved from ~/casadi-3.7.0 if not).
% • lam_v = 0 exactly is a singular point of the primer (as in tfMinEoM);
%   |lam_v| is regularized by 1e-300 only to keep the AD graph finite.
%
%% Inputs:
%
%  y                        [14 x 1]                [r; v; m; lam_r; lam_v;
%                                                   lam_m], ND rotating
%                                                   barycentric
%
%  Tmax                     double                  ND thrust acceleration
%                                                   at unit mass fraction
%
%  c                        double                  ND exhaust velocity
%
%  muStar                   double                  CR3BP mass ratio
%
%% Outputs:
%
%  F                        [14 x 1]                dy/dt
%
%  A                        [14 x 14]               dF/dy (exact, AD)
%
%  aux                      struct                  .s throttle, .alpha
%                                                   [3x1] direction, .H
%                                                   Hamiltonian (with the
%                                                   s^2 running cost),
%                                                   .sRaw the unclipped
%                                                   switching quantity
%
%% Revision History:
%  M. Casey                                                   (c) 08/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a catalog-like state; print s*, H, and the saturation check.
     mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     [F_, A_, ax_] = cr3bp_minenergy_pmp(y_, T_, c_, mu_);
     fprintf('demo: s* = %.4f, H = %.6f, |F| = %.4f, cond(A) = %.2e\n', ...
             ax_.s, ax_.H, sqrt(sum(F_.^2)), cond(A_));
     return
end

persistent fF fA fAux
if isempty(fF)
    if isempty(which('casadi.SX'))
        addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
    end
    ys = casadi.SX.sym('y', 14);  Ts = casadi.SX.sym('T');
    cs = casadi.SX.sym('c');       mus = casadi.SX.sym('mu');
    r = ys(1:3);  v = ys(4:6);  m = ys(7);
    lr = ys(8:10);  lv = ys(11:13);  lm = ys(14);
    dd = sqrt((r(1)+mus)^2 + r(2)^2 + r(3)^2);
    rr = sqrt((r(1)-1+mus)^2 + r(2)^2 + r(3)^2);
    gr = [r(1) - (1-mus)*(r(1)+mus)/dd^3 - mus*(r(1)-1+mus)/rr^3;
          r(2) - (1-mus)*r(2)/dd^3       - mus*r(2)/rr^3;
               - (1-mus)*r(3)/dd^3       - mus*r(3)/rr^3];
    hv = [2*v(2); -2*v(1); 0];
    nlv = sqrt(lv(1)^2 + lv(2)^2 + lv(3)^2 + 1e-300);
    alpha = -lv/nlv;
    sRaw = (Ts/2)*(nlv/m + lm/cs);
    s = fmin(fmax(sRaw, 0), 1);          % SX methods (dispatch, no import)
    acc = gr + hv + (s*Ts/m)*alpha;
    mdot = -s*Ts/cs;
    f7 = [v; acc; mdot];
    % costate ODEs: lamdot = -dH/dx with H = s^2 + lam'f (the s(y) dependence
    % drops out at the interior optimum by stationarity, and is inert when
    % clipped; differentiate H holding s fixed).
    lam = [lr; lv; lm];
    x7 = [r; v; m];
    Hs = s^2 + dot(lam, f7);
    lamdot = -jacobian(Hs, x7)';
    F14 = [f7; lamdot];
    A14 = jacobian(F14, ys);
    fF   = casadi.Function('F', {ys, Ts, cs, mus}, {F14}).expand();
    fA   = casadi.Function('A', {ys, Ts, cs, mus}, {A14}).expand();
    fAux = casadi.Function('aux', {ys, Ts, cs, mus}, {s, alpha, Hs, sRaw}).expand();
end

F = full(fF(y, Tmax, c, muStar));
if nargout > 1
    A = full(fA(y, Tmax, c, muStar));
end
if nargout > 2
    [s_, al_, H_, sr_] = fAux(y, Tmax, c, muStar);
    aux = struct('s', full(s_), 'alpha', full(al_), 'H', full(H_), ...
                 'sRaw', full(sr_));
end
end
