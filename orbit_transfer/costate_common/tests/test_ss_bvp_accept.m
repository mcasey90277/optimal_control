function test_ss_bvp_accept()
% TEST_SS_BVP_ACCEPT  The generic single-shooting acceptance gate.
%
% The catalog pipeline's third gate is "feed the refined entry to the
% SINGLE-shooting solver; PASS = it comes back unchanged". For min-time that
% solver is pumpkyn tfMin; for every other cost there is no pumpkyn twin, so
% ss_bvp_accept plays the role generically on the same three-closure prob.
%
% THREE TESTS on the oscillator y'' = -y, fixed tf = pi/2, y(0) = 0, hit
% y(tf) = 1 (exact y'(0) = 1):
%   1. A ROOT IS ACCEPTED: seed y'(0) = 1 -> residual at the seed ~ 1e-15,
%      |dz| < 1e-6, accepted = true.
%   2. A NEAR-MISS IS NOT: seed y'(0) = 1 + 1e-3 converges to the root but
%      MOVED by 1e-3 -> accepted = false, and the returned z is the root.
%   3. FREE-TF FORM WORKS TOO: the free-tf oscillator (hit y = 1, y' = 0)
%      from its exact (y'(0), tf) = (1, pi/2) is accepted.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per test; errors on failure)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

A = [0 1; -1 0];
prop = @(dt, y0, nS) demo_prop(A, dt, y0, nS);
prF = struct('ny', 2, 'freeIdx0', 2, 'prop', prop, 'rhs', @(y) A*y, ...
             'terminal', @(y, nJ) deal(y(1) - 1, [1 0]));
tf = pi/2;

% --- 1: exact root ---------------------------------------------------------
[z, info] = ss_bvp_accept(prF, [0; 1], tf, struct('fixedTf', true));
check('root: residual at seed ~ 0', info.normR0 < 1e-12, sprintf('%.1e', info.normR0));
check('root: |dz| < 1e-6', info.dz < 1e-6, sprintf('%.1e', info.dz));
check('root: accepted', info.accepted);
check('root: z returned = y1(free)', numel(z) == 1 && abs(z - 1) < 1e-8);

% --- 2: near miss ----------------------------------------------------------
[z2, info2] = ss_bvp_accept(prF, [0; 1 + 1e-3], tf, struct('fixedTf', true));
check('miss: converged to root', info2.converged && abs(z2 - 1) < 1e-8);
check('miss: moved ~1e-3', abs(info2.dz - 1e-3) < 1e-6, sprintf('dz %.2e', info2.dz));
check('miss: NOT accepted', ~info2.accepted);

% --- 3: free tf ------------------------------------------------------------
prT = prF;  prT.terminal = @(y, nJ) deal([y(1)-1; y(2)], eye(2));
[z3, info3] = ss_bvp_accept(prT, [0; 1], pi/2, struct());
check('free: accepted', info3.accepted, sprintf('dz %.1e', info3.dz));
check('free: z = [y''(0); tf]', numel(z3) == 2 && abs(z3(2) - pi/2) < 1e-8);
fprintf('test_ss_bvp_accept: ALL PASS\n');
end

function check(label, cond, detail)
if nargin < 3, detail = ''; end
if cond
    fprintf('  [PASS] %s %s\n', label, detail);
else
    fprintf('  [FAIL] %s %s\n', label, detail);
    error('test_ss_bvp_accept:fail', '%s', label);
end
end

function [yh, PHI] = demo_prop(A, dt, y0, needSTM)
E = expm(A*dt);
yh = E*y0;
if needSTM, PHI = E; else, PHI = []; end
end
