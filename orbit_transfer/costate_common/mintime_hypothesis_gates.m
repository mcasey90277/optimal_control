function g = mintime_hypothesis_gates(z8, rv0, Tmax, c, muStar, opts)
%% Purpose:
%
%   Per-entry checks of the hypotheses under which the free-time conjugate
%   test (ms_conjugate_test, BCT form) is a SUFFICIENT second-order
%   certificate for a min-time catalog entry -- the gates the audit
%   (doc/mintime_second_order_audit.tex, section 4) found missing:
%
%     H2  strong Legendre:  H_{alpha alpha}|_{T S^2} = (T/m)|lam_v| I > 0
%         <=> |lam_v(t)| > 0 on [0, tf].            -> minLamV, tMinLamV
%     H3  all-burn is the PMP control: s* = 1 <=> Q_mt = |lam_v|/m + lam_m/c
%         > 0 on [0, tf] (pumpkyn tfMinEoM's own switching function).
%                                                    -> minQmt, tMinQmt
%     H1' no abnormal lift of the SAME trajectory. Lifts form the linear
%         space S = { lam : lam' = -A(t)' lam along the flown x(t) with the
%         control FROZEN, lam_v(t) || alpha(t) for all t, lam_m(tf) = 0 }.
%         lam.f is a linear functional on S (conserved along each lift)
%         whose kernel is the abnormal lifts; our normal lift has
%         lam.f = -1. Hence: no abnormal lift <=> dim S = 1. Computed as
%         7 - rank(C), C = [ [alpha_k]_x Psi_v(t_k) ; e_m' Psi(tf) ] over
%         samples t_k, Psi the fixed-control adjoint fundamental matrix.
%                                                    -> dimS, svRatio, sv
%         Self-consistency: the accepted lam0 must be in S (nullResid) and
%         lam.f == -1 along the arc (Hresid).
%
%   The arc is flown with pumpkyn tfMinProp from z8 (as seed_from_z8 does),
%   so the samples are the propagator's own; gate values are sample minima
%   (junction resolution is NOT assumed -- the dense flight is used).
%
%  ASSUMPTIONS / NOTES:
%
% • The fixed-control Jacobian A(t) = df/dx|_{alpha(t), s = 1} is built by
%   CasADi from the same CR3BP expressions as cr3bp_minfuel_pmp; alpha(t)
%   is the FLOWN direction -lam_v/|lam_v| of the accepted lift.
% • dim S is a numerical rank: sv(7)/sv(6) of C reported; dimS = 1 when
%   sv(7) < rankTol*sv(1) and sv(6) >= rankTol*sv(1).
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf], tfMin
%                                                   convention
%
%  rv0                      [6 x 1]                 Departure state, ND
%
%  Tmax, c, muStar          double                  As tfMinProp
%
%  opts                     struct (optional)       .nSamp [200] constraint
%                                                   samples, .rankTol [1e-8]
%
%% Outputs:
%
%  g                        struct                  .minLamV .tMinLamV
%                                                   .minQmt .tMinQmt .dimS
%                                                   .sv [7x1] .svRatio
%                                                   .nullResid .Hresid
%                                                   .nSwitchFlown (samples
%                                                   with Q_mt <= 0)
%
%% Revision History:
%  M. Casey                                                   (c) 09/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
nSamp   = 200;   if isfield(opts, 'nSamp'),   nSamp   = opts.nSamp;   end
rankTol = 1e-8;  if isfield(opts, 'rankTol'), rankTol = opts.rankTol; end

z8 = z8(:);  tf = z8(8);
[tau, Y] = pumpkyn.cr3bp.tfMinProp(tf, [rv0(:); 1; z8(1:7)], Tmax, c, muStar);
[tau, iu] = unique(tau);  Y = Y(iu, :);

% --- H2, H3 on the dense flight --------------------------------------------
lamV = Y(:, 11:13);  rho = sqrt(sum(lamV.^2, 2));
Qmt  = rho ./ Y(:, 7) + Y(:, 14) / c;
[g.minLamV, i1] = min(rho);  g.tMinLamV = tau(i1);
[g.minQmt,  i2] = min(Qmt);  g.tMinQmt  = tau(i2);
g.nSwitchFlown = nnz(Qmt <= 0);

% --- normal-lift check: lam.f == -1 along the arc ---------------------------
F = zeros(size(Y, 1), 7);
for k = 1:size(Y, 1)
    F(k, :) = fixedField(Y(k, :)', Tmax, c, muStar)';
end
lamF = sum(Y(:, 8:14) .* F, 2);
g.Hresid = max(abs(lamF + 1));

% --- H1': dimension of the lift space S -------------------------------------
Afun = fixedJacobian();                       % persistent CasADi Function
xOf  = @(t) interp1(tau, Y(:, 1:7),   t, 'pchip')';
aOf  = @(t) -interp1(tau, lamV, t, 'pchip')';
psiRhs = @(t, psi) reshape(-full(Afun(xOf(t), unitv(aOf(t)), Tmax, c, muStar))' * reshape(psi, 7, 7), [], 1);
ts = linspace(0, tf, nSamp);
[~, PSI] = ode113(psiRhs, ts, reshape(eye(7), [], 1), odeset('RelTol', 1e-10, 'AbsTol', 1e-12));
C = zeros(3*nSamp + 1, 7);
for k = 1:nSamp
    Psi = reshape(PSI(k, :), 7, 7);
    a = unitv(aOf(ts(k)));
    C(3*k-2:3*k, :) = skew(a) * Psi(4:6, :);         % lam_v(t_k) x alpha_k = 0
end
PsiT = reshape(PSI(end, :), 7, 7);
C(end, :) = PsiT(7, :);                              % lam_m(tf) = 0
sv = svd(C);
g.sv = sv;
g.svRatio = sv(7) / max(sv(6), realmin);
g.dimS = nnz(sv < rankTol * sv(1));
g.nullResid = norm(C * z8(1:7)) / norm(z8(1:7));
end

% ------------------------------------------------------------------------
function f = fixedField(y, Tmax, c, muStar)
% FIXEDFIELD  State dynamics with the control frozen at the flown direction
% (s = 1).  INPUTS: y [14x1]; params.  OUTPUTS: f [7x1].
r = y(1:3);  v = y(4:6);  m = y(7);  lv = y(11:13);
al = -lv / sqrt(sum(lv.^2) + 1e-300);
dd = sqrt((r(1)+muStar)^2 + r(2)^2 + r(3)^2);
rr = sqrt((r(1)-1+muStar)^2 + r(2)^2 + r(3)^2);
gr = [r(1) - (1-muStar)*(r(1)+muStar)/dd^3 - muStar*(r(1)-1+muStar)/rr^3;
      r(2) - (1-muStar)*r(2)/dd^3           - muStar*r(2)/rr^3;
           - (1-muStar)*r(3)/dd^3           - muStar*r(3)/rr^3];
hv = [2*v(2); -2*v(1); 0];
f = [v; gr + hv + (Tmax/m)*al; -Tmax/c];
end

function Afun = fixedJacobian()
% FIXEDJACOBIAN  CasADi Function A(x, alpha, T, c, mu) = df/dx with the
% control FROZEN (alpha an input, s = 1).  INPUTS: none.  OUTPUTS: Afun.
persistent fA
if isempty(fA)
    if isempty(which('casadi.SX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end
    xs = casadi.SX.sym('x', 7);  as = casadi.SX.sym('a', 3);
    Ts = casadi.SX.sym('T');  cs = casadi.SX.sym('c');  mus = casadi.SX.sym('mu');
    r = xs(1:3);  v = xs(4:6);  m = xs(7);
    dd = sqrt((r(1)+mus)^2 + r(2)^2 + r(3)^2);
    rr = sqrt((r(1)-1+mus)^2 + r(2)^2 + r(3)^2);
    gr = [r(1) - (1-mus)*(r(1)+mus)/dd^3 - mus*(r(1)-1+mus)/rr^3;
          r(2) - (1-mus)*r(2)/dd^3       - mus*r(2)/rr^3;
               - (1-mus)*r(3)/dd^3       - mus*r(3)/rr^3];
    hv = [2*v(2); -2*v(1); 0];
    f7 = [v; gr + hv + (Ts/m)*as; -Ts/cs];
    fA = casadi.Function('Afix', {xs, as, Ts, cs, mus}, {jacobian(f7, xs)}).expand();
end
Afun = fA;
end

function u = unitv(x)
% UNITV  x/|x|.  INPUTS: x [3x1].  OUTPUTS: u [3x1].
u = x / sqrt(sum(x.^2) + 1e-300);
end

function S = skew(a)
% SKEW  Cross-product matrix [a]_x.  INPUTS: a [3x1].  OUTPUTS: S [3x3].
S = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
