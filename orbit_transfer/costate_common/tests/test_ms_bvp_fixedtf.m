function test_ms_bvp_fixedtf()
% TEST_MS_BVP_FIXEDTF  The fixed-final-time variant of the ms_bvp engine.
%
% THREE TESTS on the linear oscillator y'' = -y (exact flow available):
%   1. FIXED TF SOLVES. y(0)=0 fixed, y'(0) free, tf = pi/2 FIXED, hit
%      y(tf) = 1 (ONE terminal condition for ONE unknown). Exact: y'(0) = 1.
%      The unknown vector must carry no tf entry, and info.tGrid must end
%      at the fixed tf.
%   2. FIXED TF IS NOT A GUARD. Feed a seed whose .tf is deliberately wrong
%      for the free problem (tf = 1.0, where y(1) = sin(1) != 1); with
%      fixedTf the engine must return the fixed-tf answer y'(0) = 1/sin(1),
%      not move tf.
%   3. DEFAULT UNCHANGED. Without the flag the free-tf demo still returns
%      (y'(0), tf) = (1, pi/2) -- the two paths must not interfere.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per test; errors on failure)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

A = [0 1; -1 0];
prop = @(dt, y0, nS) demo_prop(A, dt, y0, nS);
K = 4;

% --- 1: fixed tf = pi/2, terminal y(tf) = 1 only ------------------------
pr = struct('ny', 2, 'freeIdx0', 2, 'prop', prop, 'rhs', @(y) A*y, ...
            'terminal', @(y, nJ) deal(y(1) - 1, [1 0]));
tf = pi/2;  tg = linspace(0, tf, K+1);
sd = struct('tf', tf, 'tGrid', tg, 'Y', [0.9*sin(tg); 0.9*cos(tg)]);
[p, info] = ms_bvp(pr, sd, struct('tolR', 1e-10, 'fixedTf', true));
check('fixed tf: unknown count nf + ny(K-1)', numel(p) == 1 + 2*(K-1));
check('fixed tf: y''(0) = 1', abs(p(1) - 1) < 1e-8, sprintf('got %.10f', p(1)));
check('fixed tf: converged', info.converged);
check('fixed tf: tGrid ends at tf', abs(info.tGrid(end) - tf) < 1e-15);

% --- 2: wrong seed.tf stays put ------------------------------------------
tf2 = 1.0;  tg2 = linspace(0, tf2, K+1);
sd2 = struct('tf', tf2, 'tGrid', tg2, 'Y', [sin(tg2); cos(tg2)]);
[p2, info2] = ms_bvp(pr, sd2, struct('tolR', 1e-10, 'fixedTf', true));
check('fixed tf: y''(0) = 1/sin(1)', abs(p2(1) - 1/sin(1)) < 1e-8, ...
      sprintf('got %.10f', p2(1)));
check('fixed tf: tf untouched', abs(info2.tGrid(end) - tf2) < 1e-15);

% --- 3: default free-tf path unchanged -----------------------------------
pr3 = struct('ny', 2, 'freeIdx0', 2, 'prop', prop, 'rhs', @(y) A*y, ...
             'terminal', @(y, nJ) deal([y(1)-1; y(2)], eye(2)));
tg3 = linspace(0, 1.2, K+1);
sd3 = struct('tf', 1.2, 'tGrid', tg3, 'Y', [sin(tg3); cos(tg3)]);
[p3, info3] = ms_bvp(pr3, sd3, struct('tolR', 1e-8));
check('free tf: y''(0) = 1', abs(p3(1) - 1) < 1e-7);
check('free tf: tf = pi/2', abs(p3(end) - pi/2) < 1e-7);
check('free tf: converged', info3.converged);
fprintf('test_ms_bvp_fixedtf: ALL PASS\n');
end

function check(label, cond, detail)
if nargin < 3, detail = ''; end
if cond
    fprintf('  [PASS] %s\n', label);
else
    fprintf('  [FAIL] %s %s\n', label, detail);
    error('test_ms_bvp_fixedtf:fail', '%s', label);
end
end

function [yh, PHI] = demo_prop(A, dt, y0, needSTM)
E = expm(A*dt);
yh = E*y0;
if needSTM, PHI = E; else, PHI = []; end
end
