function [F, A, aux] = cr3bp_minfuel_pmp(y, Tmax, c, muStar, smooth)
%% Purpose:
%
%   The SMOOTHED ENERGY->FUEL PMP vector field for the CR3BP low-thrust
%   transfer, with its exact 14x14 Jacobian -- the homotopy sibling of
%   cr3bp_minenergy_pmp. Same 14-state PMP y = [r; v; m; lam_r; lam_v;
%   lam_m], same dynamics and costate ODEs; only the running cost L(s)
%   changes, selected by the smoothing family:
%
%     'eps'  (Bertrand-Epenoy, the GTO_tulip direct-campaign convention):
%            L = s - p*s(1-s) = (1-p)s + p*s^2,  p in (0, 1]
%            p=1 -> s^2 (the min-energy field EXACTLY); p->0 -> s (fuel).
%            Throttle: s* = clip( (Q - (1-p)) / (2p), 0, 1 ) -- continuous
%            in Q for every p > 0, slope 1/(2p).
%
%     'huber' (the PLQ arm, Mike 2026-08-31): L = huber_p(s) =
%            s^2/(2p) for s <= p,  s - p/2 for s > p  (p = kappa knee).
%            Throttle: s* = p*Q for Q < 1, s* = 1 for Q > 1 -- the affine
%            tail makes H(s) DEGENERATE at Q = 1, so the law JUMPS from p
%            to 1 there for every p < 1 (recorded theory, tested). p->0 ->
%            fuel. Note also: huber NEVER COASTS exactly -- s* = p*Q > 0
%            below the switch (no dead zone; bang-limit convergence is
%            first-order in p), whereas 'eps' produces true s = 0 coast
%            arcs (Q <= 1-p) at every finite p. The race (step-5 Task 3)
%            measures whether these structural differences matter in
%            continuation practice.
%
%   with Q = Tmax(|lam_v|/m + lam_m/c) the switch quantity: the fuel
%   (p -> 0) bang-bang law is s = 1 where Q > 1, s = 0 where Q < 1.
%
%       H = L(s) + lam_r'v + lam_v'(g + h + s T alpha/m) - lam_m s T/c
%       alpha* = -lam_v/|lam_v|                       (primer direction)
%
%   The Jacobian is built ONCE per family by CasADi SX automatic
%   differentiation (persistent Functions; the smoothing parameter is a
%   Function INPUT, so one build serves a whole continuation walk).
%
%  ASSUMPTIONS / NOTES:
%
% • p = 0 exactly is the nonsmooth bang-bang limit -- not representable as
%   a smooth field; walk p down and certify the limit separately.
% • lam_v = 0 is regularized by 1e-300 exactly as in the energy field.
% • Requires CasADi (self-resolved from ~/casadi-3.7.0 if absent).
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
%  smooth                   struct                  .family 'eps'|'huber',
%                                                   .p the smoothing
%                                                   parameter (eps in
%                                                   (0,1], or the Huber
%                                                   knee kappa > 0)
%
%% Outputs:
%
%  F                        [14 x 1]                dy/dt
%
%  A                        [14 x 14]               dF/dy (exact, AD)
%
%  aux                      struct                  .s throttle, .alpha
%                                                   [3x1], .H (with L(s)),
%                                                   .Q switch quantity,
%                                                   .sRaw the family's
%                                                   unclipped stationary
%                                                   throttle
%
%% Revision History:
%  M. Casey                                                   (c) 09/01/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: one state, both families at a mid parameter:
     mu_ = 0.012150585609624;  T_ = 0.1756418;  c_ = 8.673746;
     y_ = [0.85; 0.05; 0.01; 0.05; 0.55; -0.02; 0.97; ...
           15.6; 32.9; -0.09; -0.10; 0.045; -0.00015; 0.13];
     for fam = {struct('family','eps','p',0.3), struct('family','huber','p',0.3)}
         [F_, ~, ax_] = cr3bp_minfuel_pmp(y_, T_, c_, mu_, fam{1});
         fprintf('demo %-5s p=0.3: s* = %.4f, Q = %.4f, H = %.6f, |F| = %.4f\n', ...
                 fam{1}.family, ax_.s, ax_.Q, ax_.H, sqrt(sum(F_.^2)));
     end
     return
end

assert(isstruct(smooth) && all(isfield(smooth, {'family','p'})) && ...
       smooth.p > 0, 'cr3bp_minfuel_pmp: smooth must be {family, p > 0}');

persistent fn                              % fn.(family) = {fF, fA, fAux}
if isempty(fn) || ~isfield(fn, smooth.family)
    if isempty(which('casadi.SX'))
        addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
    end
    ys = casadi.SX.sym('y', 14);  Ts = casadi.SX.sym('T');
    cs = casadi.SX.sym('c');      mus = casadi.SX.sym('mu');
    ps = casadi.SX.sym('p');
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
    Q = Ts*(nlv/m + lm/cs);
    switch smooth.family
        case 'eps'
            sRaw = (Q - (1 - ps)) / (2*ps);
            s = fmin(fmax(sRaw, 0), 1);
            Lc = (1 - ps)*s + ps*s^2;
        case 'huber'
            % argmin of huber_p(s) - s*Q on [0,1]: p*Q on the quadratic
            % core for Q < 1, saturation for Q > 1 (degenerate AT Q = 1;
            % the min-norm selection p*Q is taken there):
            sRaw = ps*Q;
            s = if_else(Q > 1, 1, fmin(fmax(sRaw, 0), ps));
            Lc = if_else(s <= ps, s^2/(2*ps), s - ps/2);
        otherwise
            error('cr3bp_minfuel_pmp:family', ...
                  'unknown smoothing family ''%s''', smooth.family);
    end
    acc = gr + hv + (s*Ts/m)*alpha;
    mdot = -s*Ts/cs;
    f7 = [v; acc; mdot];
    % costate ODEs: lamdot = -dH/dx holding s fixed (interior stationarity
    % kills the s(y) dependence; the clipped/saturated branches are inert)
    lam = [lr; lv; lm];
    x7 = [r; v; m];
    Hs = Lc + dot(lam, f7);
    lamdot = -jacobian(Hs, x7)';
    F14 = [f7; lamdot];
    A14 = jacobian(F14, ys);
    fn.(smooth.family) = { ...
        casadi.Function('F',   {ys, Ts, cs, mus, ps}, {F14}).expand(), ...
        casadi.Function('A',   {ys, Ts, cs, mus, ps}, {A14}).expand(), ...
        casadi.Function('aux', {ys, Ts, cs, mus, ps}, {s, alpha, Hs, Q, sRaw}).expand()};
end

fset = fn.(smooth.family);
F = full(fset{1}(y, Tmax, c, muStar, smooth.p));
if nargout > 1
    A = full(fset{2}(y, Tmax, c, muStar, smooth.p));
end
if nargout > 2
    [s_, al_, H_, Q_, sr_] = fset{3}(y, Tmax, c, muStar, smooth.p);
    aux = struct('s', full(s_), 'alpha', full(al_), 'H', full(H_), ...
                 'Q', full(Q_), 'sRaw', full(sr_));
end
end
