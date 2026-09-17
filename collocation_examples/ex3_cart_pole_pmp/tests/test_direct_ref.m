function ok = test_direct_ref()
%% Purpose:
%
%   The committed direct-solution fixture is a real solution of the problem
%   the indirect demo solves: right shape, its own trapezoidal defects
%   small, and a control that stayed inside the relaxed bound (so it is the
%   UNCONSTRAINED optimum the PMP solve is compared to). The trapezoidal
%   defect check is the check that measures dynamics accuracy: it recomputes
%   the derivative field from the independently authored
%   ex2_cart_pole_swing_up/try2 helpers (cart_accel, pendulum_accel), not
%   from a copy of gen_direct_ref.m's own dynamics, so an error shared by
%   both would not cancel out. The initial- and terminal-state checks below
%   measure something different: those nodes are pinned by equality
%   constraints in gen_direct_ref.m's NLP, not propagated by the dynamics,
%   so a near-exact match there is evidence the boundary rows were wired up
%   correctly, not evidence of integration accuracy.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(here), 'ex2_cart_pole_swing_up', 'try2'));
R = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));

ok = chk(ok, isequal(size(R.X), [4 201]) && isequal(size(R.U), [1 201]) ...
             && isequal(size(R.muDefect), [4 200]) && numel(R.tN) == 201, ...
         'shapes: X [4 x 201], U [1 x 201], muDefect [4 x 200]');
ok = chk(ok, R.tN(1) == 0 && abs(R.tN(end) - 5) < 1e-12 && abs(R.tf - 5) < 1e-12, ...
         'the grid spans the fixed 5 s horizon');
ok = chk(ok, max(abs(R.X(:,1) - [0;0;0;0])) < 1e-8, ...
         'initial state matches the pinned target [0;0;0;0] (an equality-constraint node, not propagated)');
ok = chk(ok, max(abs(R.X(:,end) - [0;pi;0;0])) < 1e-6, ...
         sprintf(['terminal state matches the pinned target [0;pi;0;0] (an equality-' ...
                  'constraint node, not propagated; miss %.1e)'], ...
                 max(abs(R.X(:,end) - [0;pi;0;0]))));
ok = chk(ok, max(abs(R.U)) < 0.9*R.uMax, ...
         sprintf('the control stayed off the relaxed bound: max |u| = %.1f N of %.0f N', ...
                 max(abs(R.U)), R.uMax));

% its own trapezoidal defects, recomputed here from the stored nodes using
% the example's OWN, independently authored dynamics helpers (not a copy
% of gen_direct_ref.m's internal ref_field) so a shared error cannot cancel
h = diff(R.tN);
F = zeros(4, 201);
for k = 1:201
    q2 = R.X(2,k);  q1dot = R.X(3,k);  q2dot = R.X(4,k);  u = R.U(k);
    F(:,k) = [q1dot; q2dot; ...
              cart_accel(q2, q2dot, u, R.p.L, R.p.m1, R.p.m2, R.p.g); ...
              pendulum_accel(q2, q2dot, u, R.p.L, R.p.m1, R.p.m2, R.p.g)];
end
d = R.X(:,2:end) - R.X(:,1:end-1) - (h/2).*(F(:,2:end) + F(:,1:end-1));
% Gate 1e-12, not the original 1e-6: the committed fixture measures 6.5e-14
% here, so 1e-6 would have accepted a fixture eight orders of magnitude
% worse than the one actually shipped. 1e-12 keeps ~2 orders of margin
% above the measured value (cross-version/BLAS noise) while refusing
% anything that regressed materially.
ok = chk(ok, max(abs(d(:))) < 1e-12, sprintf('trapezoidal defects small: %.1e', max(abs(d(:)))));
ok = chk(ok, R.J > 0 && isfinite(R.J), sprintf('a finite positive cost: J = %.6f', R.J));

if ok, fprintf('TEST_DIRECT_REF: ALL PASS\n');
else,  fprintf('TEST_DIRECT_REF: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
