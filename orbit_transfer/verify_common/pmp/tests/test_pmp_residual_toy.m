% TEST_PMP_RESIDUAL_TOY  Validate pmp_residual against an ANALYTIC extremal.
%
% This is Phase 3 Step P3.1 and it GATES everything downstream: if the checker
% cannot certify a solution whose exact answer is known in closed form, no
% residual it reports on a campaign row is interpretable.
%
% THE PROBLEM (the same double integrator as test_foc_check_toy):
%   min int_0^T u dt,  xdot = v,  vdot = u,  u in [0,1],
%   x(0) = v(0) = 0,   x(T) = 1,  v(T) FREE,  T = 2.
%
% THE ANALYTIC EXTREMAL, derived rather than fitted:
%   H = u + lam_x*v + lam_v*u  =>  lam_x_dot = 0, lam_v_dot = -lam_x.
%   So lam_x = c1 (const) and lam_v(t) = lam_v(T) + c1*(T-t).
%   v(T) is free with no terminal cost, so transversality gives lam_v(T) = 0
%   and hence lam_v(t) = c1*(T-t).
%   Switching function S = 1 + lam_v; the minimum condition gives u = 1 where
%   S < 0 and u = 0 where S > 0, i.e. burn first, then coast.
%   With u = 1 on [0,ts]: v = t, and thereafter v = ts, so
%   x(T) = ts^2/2 + ts*(T-ts) = ts*T - ts^2/2 = 1.
%   At T = 2 this is ts^2 - 4*ts + 2 = 0, so ts = 2 - sqrt(2).
%   S(ts) = 0 gives c1 = 1/(ts - T) = -1/sqrt(2).
%
%   Hamiltonian check: on the burn arc H = 1 + c1*T; on the coast arc
%   H = c1*ts. Both equal 1 - sqrt(2) = -0.414214, so H is constant across the
%   switch, as it must be for an autonomous fixed-horizon problem.
%
% TWO PARTS, because a checker must do two things:
%   (a) CORRECTNESS -- on a grid whose nodes INCLUDE the switch time, so no
%       interval straddles the control discontinuity, the residuals of the
%       exact solution must sit at integrator tolerance. Anything larger is a
%       bug in the checker itself.
%   (b) SENSITIVITY -- on a grid that does NOT contain the switch, exactly one
%       interval straddles it, and the checker must produce a residual spike
%       THERE and nowhere else. This is the phenomenon the whole study is
%       about: an unaligned mesh carries the wrong control over part of one
%       interval. A checker that reports (a) cleanly but misses (b) would be
%       silently blind to the only error that matters.
root = fileparts(fileparts(mfilename('fullpath')));  cd(root);  addpath(root);
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
import casadi.*

T   = 2;
ts  = 2 - sqrt(2);            % switch time
c1  = 1/(ts - T);             % = -1/sqrt(2)
Hex = 1 + c1*T;               % constant Hamiltonian, = 1 - sqrt(2)
fprintf('analytic: ts = %.12f  c1 = %.12f  H = %.12f\n', ts, c1, Hex);

% --- the transcription's dynamics and running cost, as CasADi Functions -----
xs = MX.sym('x',2);  us = MX.sym('u',1);
fFun = Function('f', {xs,us}, {[xs(2); us]});
LFun = Function('L', {xs,us}, {us});
uPMP = @(x,l) double( (1 + l(2)) < 0 );      % minimum condition: u = 1 iff S<0

% --- exact solution sampled on a grid ---------------------------------------
exact = @(t) local_exact(t, ts, c1, T);

% ============================ (a) SWITCH-ALIGNED ============================
% Evaluated PER ARC, not on one grid spanning the switch. At a bang-bang switch
% the control is genuinely double-valued -- u = 1 as the burn arc ends and
% u = 0 as the coast arc begins -- so a single shared node cannot carry both,
% and a grid that merely CONTAINS the switch time still interpolates 1 -> 0
% across the first coast interval. (Measured: doing exactly that gives a
% spurious 1.77e-02 residual, which the self-gate correctly refuses.) A true
% retained breakpoint is a DUPLICATED node, which is what a multi-phase
% formulation provides and what per-arc evaluation reproduces here.
arcs = { [0 ts], [ts T] };
for a = 1:2
    sgA = linspace(arcs{a}(1), arcs{a}(2), 41);
    [XA, UA, LA] = local_sample(sgA, exact, a == 1);
    RA = pmp_residual(sgA, XA, UA, LA, fFun, LFun, ...
            struct('uPMP', uPMP, 'variant', 'both', 'verbose', true));

    assert(RA.gatePass, 'a%d: self-gate must pass on the exact solution (defect %.3e)', ...
        a, RA.selfDefect);
    assert(max(RA.Rx_interp) < 1e-10, ...
        'a%d: primal residual must be at integrator tolerance, got %.3e', a, max(RA.Rx_interp));
    assert(max(RA.Rlam_interp) < 1e-10, ...
        'a%d: COSTATE residual must be at integrator tolerance, got %.3e', a, max(RA.Rlam_interp));
    assert(max(RA.Rx_pmp) < 1e-10 && max(RA.Rlam_pmp) < 1e-10, ...
        'a%d: the PMP-control variant must also be tight (%.3e / %.3e)', ...
        a, max(RA.Rx_pmp), max(RA.Rlam_pmp));
    assert(abs(mean(RA.H) - Hex) < 1e-10, ...
        'a%d: Hamiltonian should be %.12f, got %.12f', a, Hex, mean(RA.H));
    assert(RA.Hvar < 1e-10, 'a%d: Hamiltonian must be CONSTANT, variation %.3e', a, RA.Hvar);
    fprintf('  (a%d) arc [%.4f %.4f]: residuals at tolerance, H = %.12f constant. PASS\n', ...
            a, arcs{a}(1), arcs{a}(2), mean(RA.H));
end

% ============================ (b) UNALIGNED =================================
% A uniform grid over [0,T]; ts = 0.5857864 falls strictly inside one interval.
sgB = linspace(0, T, 81);
kSw = find(sgB(1:end-1) < ts & sgB(2:end) > ts);
assert(isscalar(kSw), 'b: expected exactly one straddling interval, got %d', numel(kSw));
[XB, UB, LB] = local_sample(sgB, exact);
% The self-gate is EXPECTED to fire here and that is correct behaviour, not a
% failure: on an unaligned grid the trapezoid genuinely cannot represent this
% exact bang-bang solution, so a large recomputed defect is the true state of
% affairs. (For an NLP solution the gate should always pass, because the NLP
% drove those same defects to zero -- which is why it is a wiring check there.)
% defectTol is raised only so the run proceeds to the residual map.
RB = pmp_residual(sgB, XB, UB, LB, fFun, LFun, ...
        struct('uPMP', uPMP, 'variant', 'interp', 'verbose', true, 'defectTol', 1));
assert(~RB.gatePass || RB.selfDefect > 1e-10, ...
    'b: an unaligned grid MUST show a non-trivial recomputed defect');

rx = RB.Rx_interp;
others = rx([1:kSw-1, kSw+1:end]);
fprintf('  (b) straddling interval %d: R_x = %.3e ; max elsewhere = %.3e\n', ...
        kSw, rx(kSw), max(others));
assert(rx(kSw) > 1e3 * max(max(others), 1e-16), ...
    'b: the straddling interval must dominate -- got %.3e vs %.3e elsewhere', ...
    rx(kSw), max(others));
assert(max(others) < 1e-10, ...
    'b: NON-straddling intervals must still be at tolerance, got %.3e', max(others));
[~, kMax] = max(rx);
assert(kMax == kSw, 'b: the residual peak must be AT the straddling interval (%d), got %d', ...
    kSw, kMax);
fprintf('  (b) unaligned grid: spike localized at the switch, elsewhere clean. PASS\n');

fprintf('test_pmp_residual_toy: ALL PASS\n');

% ===========================================================================
function [X, U, lam] = local_sample(sg, exact, burnArc)
% LOCAL_SAMPLE  Evaluate the analytic extremal at the grid nodes.
%
% burnArc (optional) forces the control on the ARC being sampled, which is what
% makes the switch endpoint unambiguous: the same instant t = ts carries u = 1
% when it closes the burn arc and u = 0 when it opens the coast arc.
% INPUTS:  sg [1xN1]; exact - handle t -> [x; v; lam_x; lam_v; u]; burnArc
% OUTPUTS: X [2xN1]; U [1xN1]; lam [2xN1]
n = numel(sg);
X = zeros(2,n);  U = zeros(1,n);  lam = zeros(2,n);
for k = 1:n
    z = exact(sg(k));
    X(:,k)   = z(1:2);
    lam(:,k) = z(3:4);
    U(k)     = z(5);
end
if nargin >= 3 && ~isempty(burnArc)
    U(:) = double(burnArc);
end
end

% ---------------------------------------------------------------------------
function z = local_exact(t, ts, c1, T)
% LOCAL_EXACT  The closed-form extremal at time t.
%
% AT the switch the control is taken as the BURN value. The choice is
% arbitrary at a measure-zero point and does not affect the propagated
% residuals, since a node value only enters as one endpoint of the linear
% interpolant on each side.
% INPUTS: t; ts switch time; c1 costate slope; T horizon
% OUTPUTS: z = [x; v; lam_x; lam_v; u] [5x1]
if t <= ts
    v = t;             x = t^2/2;        u = 1;
else
    v = ts;            x = ts^2/2 + ts*(t - ts);  u = 0;
end
z = [x; v; c1; c1*(T - t); u];
end
