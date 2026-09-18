function out = run_cartpole_pmp(opts)
%% Purpose:
%
%   Solve the cart-pole minimum-effort swing-up INDIRECTLY: build the
%   Pontryagin boundary-value problem and shoot it with the shared multiple-
%   shooting engine oc.ms_bvp (oclib). This file needs nothing from
%   orbit_transfer: costate_common keeps its own delegate (ms_bvp.m) for
%   its own campaigns, but this demo calls oc.ms_bvp directly and adds only
%   oclib to the path.
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
%                                                   [@oc.ms_bvp], .tolR
%                                                   [1e-10], .polishMax [5];
%                                                   a caller
%                                                   who instead passes the
%                                                   costate_common delegate
%                                                   (ms_bvp) must put
%                                                   orbit_transfer/
%                                                   costate_common on the
%                                                   path themselves -- this
%                                                   file no longer does
%
%% Outputs:
%
%  out                      struct                  .ok and .why (THE front
%                                                   door's verdict: converged,
%                                                   reporting flight reached
%                                                   t_f, terminal miss < 1e-9,
%                                                   flown miss < 1e-6), .lam0,
%                                                   .J, .missTerminal (from
%                                                   a single re-flight of
%                                                   lam0), .missEngine (the
%                                                   shooting's own last-arc
%                                                   residual), .missFlown,
%                                                   .Hrel (spread
%                                                   of the Hamiltonian along
%                                                   the arc), .t, .X, .U,
%                                                   .Lam, .info, .seed
%                                                   (.signCorr .flipped
%                                                   .ampRatio)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  M. Casey  Astra review: own verdict, sign confidence, H     (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here, fullfile(root, 'oclib'), fullfile(fileparts(here), 'cartpole_common'));
d = @(f, v) fieldd(opts, f, v);
K      = d('K', 8);
engine = d('engine', @oc.ms_bvp);
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
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.muDefect, ...
                                         'tNodes', R.tN, 'velRows', 3:4));
tS = tS(:).';                                   % ROWS throughout: a column here
uImplied = zeros(1, numel(tS));                 % would turn the products below
Xs = interp1(R.tN, R.X.', tS, 'pchip').';       % into an outer product, silently
for k = 1:numel(tS)
    [~, Gk] = cartpole_field(Xs(:,k), p);
    uImplied(k) = -(lamS(:,k).'*Gk)/2;
end
uDirect = reshape(interp1(R.tN, R.U, tS, 'pchip'), 1, []);
assert(all(isfinite(uImplied)) && all(isfinite(uDirect)) && any(uImplied) && any(uDirect), ...
       'run_cartpole_pmp:sign', 'the seed or the direct control is non-finite or identically zero');
% NORMALISED correlation, with a bar: "not exactly zero" admits a vote of no
% meaning. The sign vote cannot see a wrong costate SCALE, so the amplitude
% ratio is reported beside it (1 = the mapping's magnitude is right too).
signCorr = sum(uImplied.*uDirect) / sqrt(sum(uImplied.^2)*sum(uDirect.^2));
assert(abs(signCorr) > 0.5, 'run_cartpole_pmp:sign', ...
       'implied vs direct control correlation %.3f: the seed carries no usable sign information', signCorr);
if signCorr < 0, lamS = -lamS;  uImplied = -uImplied; end
seedDiag = struct('signCorr', abs(signCorr), 'flipped', signCorr < 0, ...
                  'ampRatio', sqrt(sum(uImplied.^2)/sum(uDirect.^2)));
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
tolR     = d('tolR', 1e-10);
polishMax = d('polishMax', 5);
[~, info] = engine(prob, seed, struct('fixedTf', true, 'tolR', tolR, 'maxIter', 60, ...
                                      'polishMax', polishMax));

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
% RelTol 2.5e-14, not the propagator's 1e-12: this single shot re-flies the
% WHOLE 5 s horizon from lam0, so its own integration error lands directly
% in the terminal miss that is this demo's headline accuracy claim. Measured
% on the converged root, re-flying the identical lam0: miss 6.88e-09 at
% RelTol 1e-12, 8.20e-10 at 1e-13, and 1.475e-10 at the 2.5e-14 actually set
% below, while the ENGINE's own last-arc terminal residual is 7.268e-14 --
% i.e. the loose figure was the measurement, not the solution. The gate stays
% at 1e-9. (The 3.67e-11 this paragraph used to quote came from a probe at
% RelTol 1e-14 and was never re-taken at the shipped setting; re-measured
% 2026-09-17 alongside cartpole_minenergy_study, which prints the same
% numbers.)
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), linspace(0, tf, 2001), ...
                [0; 0; 0; 0; lam0], odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
% the reporting flight has no collapse protection of its own: a truncated
% return would otherwise be reported as "the terminal state"
reachedTf = (numel(t) == 2001) && (t(end) == tf) && all(isfinite(Y(:)));
X = Y(:,1:4).';  Lam = Y(:,5:8).';
U = zeros(1, numel(t));  H = zeros(1, numel(t));
for k = 1:numel(t)
    [Fk, G] = cartpole_field(X(:,k), p);
    U(k) = -(Lam(:,k).'*G)/2;
    H(k) = U(k)^2 + Lam(:,k).'*(Fk + G*U(k));
end
J = trapz(t, U.^2);
% H is CONSTANT (not zero: t_f is fixed) along any trajectory of this
% autonomous Hamiltonian flow, so its spread measures how well the costate
% equation was integrated -- not optimality, and not the seed's basin. The
% residual 2u + lam'G is NOT reported: U is built as -lam'G/2 on the line
% above, so it is zero by construction and measures nothing.
Hrel = (max(H) - min(H)) / max(abs(H));

% THE CONTROL AS A FUNCTION OF TIME, in closed form. Interpolating the
% sampled U put the resampling error straight into the flown miss (8.9e-07
% against a 1e-06 gate, and 3.6e-06 on a coarser grid): the interpolant, not
% the physics, was being measured. lam(t) comes from the same converged
% costate arc -- cubic-Hermite in the COSTATES, whose own derivative the flow
% supplies exactly -- and u = -lam'G/2 is then evaluated from the flown
% state, so the control the check flies is the PMP control at that instant.
ppL = pchip(t, Lam);
uOf = @(tt, x) closedFormU(ppL, min(max(tt, t(1)), t(end)), x, p);
zEnd = oc.fly_control([0; 0; 0; 0], [0 tf], ...
    @(tt, x) flyRhs(x, uOf(tt, x), p), struct('mode', 'span', 'solver', @ode113, ...
                                              'RelTol', 1e-13, 'AbsTol', 1e-15));

missTerminal = max(abs(X(:,end) - xf));
missFlown    = max(abs(zEnd - xf));
% the SHOOTING's own terminal accuracy, independent of the reporting flight:
% propagate the engine's last junction over its own interval
dtLast    = info.tGrid(end) - info.tGrid(end-1);
yLast     = cartpole_pmp_prop(dtLast, info.Y(:,end), false, p);
missEngine = max(abs(yLast(1:4) - xf));
% THE FRONT DOOR'S OWN VERDICT. It used to measure and leave judging to the
% tests; a caller now gets .ok and the reason, so an unconverged or truncated
% solve cannot be read as an answer.
why = '';
if ~info.converged,           why = sprintf('engine not converged (|R| = %.2e)', info.normR);
elseif ~reachedTf,            why = 'reporting flight truncated or non-finite';
elseif ~(missTerminal < 1e-9), why = sprintf('terminal miss %.2e', missTerminal);
elseif ~(missFlown < 1e-6),   why = sprintf('flown miss %.2e', missFlown);
end
out = struct('ok', isempty(why), 'why', why, 'lam0', lam0, 'J', J, ...
    'missTerminal', missTerminal, 'missEngine', missEngine, 'missFlown', missFlown, 'Hrel', Hrel, ...
    't', t.', 'X', X, 'U', U, 'Lam', Lam, 'info', info, 'seed', seedDiag);

if doPlot, plotAgainstDirect(out, R); end
if ~out.ok, warning('run_cartpole_pmp:notOk', 'NOT a solution: %s', why); end
if nargout == 0
    fprintf(['cart-pole PMP-BVP (%s): J = %.6f (direct %.6f), terminal miss %.2e, ', ...
             'flown miss %.2e, H spread %.2e rel\n'], ...
            tern(out.ok, 'ok', ['NOT OK: ' why]), out.J, R.J, out.missTerminal, out.missFlown, out.Hrel);
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

function u = closedFormU(ppL, tt, x, p)
%% Purpose:
%
%   The PMP control at time tt: the costate from the converged arc, G from
%   the FLOWN state, u = -lam'G/2. No interpolation of u itself.
%
[~, G] = cartpole_field(x, p);
u = -(ppval(ppL, tt).'*G)/2;
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

function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, v0)
%% Purpose:
%
%   s.(f) if present, else the default.
%
if isfield(s, f), v = s.(f); else, v = v0; end
end
