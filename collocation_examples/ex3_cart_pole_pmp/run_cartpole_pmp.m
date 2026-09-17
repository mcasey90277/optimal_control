function out = run_cartpole_pmp(opts)
%% Purpose:
%
%   Solve the cart-pole minimum-effort swing-up INDIRECTLY: build the
%   Pontryagin boundary-value problem and shoot it with the shared multiple-
%   shooting engine (costate_common/ms_bvp, or oclib's oc.ms_bvp).
%
%   Four unknowns -- lam(0) -- against four terminal conditions, at fixed
%   final time. The seed comes from the committed direct solution's defect
%   multipliers, mapped to costates by oc.duals_to_costates with the
%   TRAPEZOID station rule (the CR3BP campaigns use Hermite-Simpson, so this
%   demo exercises a rule the orbit work never does). Its primer sign vote
%   does not apply to a scalar force, so the sign is resolved here against
%   the direct solve's own control.
%
%   This file is also the second TOP-LEVEL consumer of the oclib package,
%   which is what admits ms_bvp to it: nothing here is an orbit, a CR3BP
%   quantity, or a pumpkyn call.
%
%% Inputs:
%
%  opts                     struct (optional)       .K segments [8], .plot
%                                                   [true when nargout = 0],
%                                                   .engine solver handle
%                                                   [@ms_bvp]
%
%% Outputs:
%
%  out                      struct                  .lam0, .J, .missTerminal,
%                                                   .missFlown, .statMax,
%                                                   .t, .X, .U, .Lam, .info
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here, fullfile(root, 'oclib'), fullfile(root, 'orbit_transfer', 'costate_common'));
d = @(f, v) fieldd(opts, f, v);
K      = d('K', 8);
engine = d('engine', @ms_bvp);
doPlot = d('plot', nargout == 0);

R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
p = R.p;  tf = R.tf;
xf = [0; pi; 0; 0];

%% Seed: the direct solve's multipliers become costates
%  duals_to_costates resolves the GLOBAL SIGN by a primer vote over 3-vector
%  thrust directions. This problem's control is a scalar force, so there is
%  no such vote: the mapping is asked for the costates unsigned (no uDir),
%  and the sign is fixed here by the same principle in scalar form -- the
%  implied control u = -lam'G/2 must agree in sign with the control the
%  direct solve actually used.
[lamS, tS, dg] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.muDefect, ...
                                             'tNodes', R.tN, 'velRows', 3:4));
uImplied = zeros(1, numel(tS));
Xs = interp1(R.tN, R.X.', tS, 'pchip').';
for k = 1:numel(tS)
    [~, Gk] = cartpole_field(Xs(:,k), p);
    uImplied(k) = -(lamS(:,k).'*Gk)/2;
end
uDirect = interp1(R.tN, R.U, tS, 'pchip');
if sum(uImplied.*uDirect) < 0, lamS = -lamS; end
assert(sum(uImplied.*uDirect) ~= 0, 'run_cartpole_pmp:sign', ...
       'the implied and solved controls are orthogonal: the seed carries no sign information');
tGrid = linspace(0, tf, K+1);
Xg   = interp1(R.tN, R.X.', tGrid, 'pchip').';
Lg   = interp1(tS, lamS.', tGrid, 'pchip', 'extrap').';
seed = struct('tf', tf, 'tGrid', tGrid, 'Y', [Xg; Lg]);
seed.Y(1:4,1) = [0; 0; 0; 0];                  % the fixed departure state

%% The problem, as three closures
prob = struct('ny', 8, 'freeIdx0', 5:8, ...
    'prop', @(dt, y0, needSTM) cartpole_pmp_prop(dt, y0, needSTM, p), ...
    'rhs',  @(y) cartpole_pmp_rhs(y, p), ...
    'terminal', @(y, needJ) terminalFcn(y, xf, needJ));

%  tolR = 1e-10, not 1e-12: this residual is ABSOLUTE while the costates
%  are O(1e3), so 1e-10 already sits at a relative level of ~1e-13.
%  Tighter than that chases noise: the propagator's own RelTol (1e-12,
%  Task 5) floors the achievable multiple-shooting residual near 4e-12
%  regardless of iteration budget (measured: polishMax=50, maxIter=300
%  reproduce the identical 3.865e-12 floor as the default settings). The
%  accuracy claim this demo makes lives on the terminal-miss gate below
%  (1e-9), not on this internal residual.
[~, info] = engine(prob, seed, struct('fixedTf', true, 'tolR', 1e-10, 'maxIter', 60));

%% Report: fly the answer, measure what the gates measure
%  NOTE (deviation from the brief's literal 501-point grid, DO NOT
%  "optimise" this back down): the reported trajectory doubles as the
%  source of U(t) for the flown-control check below, which samples u only
%  at these output times and pchip-interpolates between them. At 501
%  points that interpolation error alone -- not the physics -- dominates
%  the flown miss. Measured against the SAME converged lam0: 3.6e-6 at
%  501 points (fails the 1e-6 gate), 1.1e-8 at 2001, 5.4e-11 at 50001 --
%  monotone convergence with grid density, i.e. a pure resampling
%  artifact. 2001 points removes it with a wide margin while changing no
%  other gate materially (J, control RMS): see the task report.
lam0 = info.Y(5:8, 1);
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), linspace(0, tf, 2001), ...
                [0; 0; 0; 0; lam0], odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
X = Y(:,1:4).';  Lam = Y(:,5:8).';
U = zeros(1, numel(t));  stat = zeros(1, numel(t));
for k = 1:numel(t)
    [~, G] = cartpole_field(X(:,k), p);
    U(k)    = -(Lam(:,k).'*G)/2;
    stat(k) = abs(2*U(k) + Lam(:,k).'*G);
end
J = trapz(t, U.^2);

uOf = @(tt) interp1(t, U, min(max(tt, t(1)), t(end)), 'pchip');
zEnd = oc.fly_control([0; 0; 0; 0], [0 tf], ...
    @(tt, x) flyRhs(x, uOf(tt), p), struct('mode', 'span', 'solver', @ode113));

out = struct('lam0', lam0, 'J', J, ...
    'missTerminal', max(abs(X(:,end) - xf)), ...
    'missFlown', max(abs(zEnd - xf)), 'statMax', max(stat), ...
    't', t.', 'X', X, 'U', U, 'Lam', Lam, 'info', info);

if doPlot, plotAgainstDirect(out, R); end
if nargout == 0
    fprintf(['cart-pole PMP-BVP: J = %.6f (direct %.6f), terminal miss %.2e, ', ...
             'flown miss %.2e, stationarity %.2e\n'], ...
            out.J, R.J, out.missTerminal, out.missFlown, out.statMax);
    clear out
end
end

% ------------------------------------------------------------------------
function [g, dgdy] = terminalFcn(y, xf, needJ)
%% Purpose:
%
%   The four terminal conditions x(tf) = xf, and their Jacobian.
%
g = y(1:4) - xf;
dgdy = [];
if needJ, dgdy = [eye(4), zeros(4)]; end
end

function dx = flyRhs(x, u, p)
%% Purpose:
%
%   The TRUE dynamics with a supplied control -- what the flown check flies.
%
[F, G] = cartpole_field(x, p);
dx = F + G*u;
end

function plotAgainstDirect(out, R)
%% Purpose:
%
%   The indirect solution against the direct one: states and control.
%
figure('color', [1 1 1]);
subplot(2,1,1); plot(out.t, out.X(2,:), 'k', R.tN, R.X(2,:), 'r--'); grid on
ylabel('q_2 (rad)'); legend({'indirect (PMP-BVP)', 'direct (collocation)'}, 'Location', 'best');
title('cart-pole swing-up: indirect vs direct');
subplot(2,1,2); plot(out.t, out.U, 'k', R.tN, R.U, 'r--'); grid on
xlabel('t (s)'); ylabel('u (N)');
end

function v = fieldd(s, f, v0)
%% Purpose:
%
%   s.(f) if present, else the default.
%
if isfield(s, f), v = s.(f); else, v = v0; end
end
