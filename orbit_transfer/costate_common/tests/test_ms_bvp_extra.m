function ok = test_ms_bvp_extra()
%% Purpose:
%
%   Tests ms_bvp's EXTRA SCALAR UNKNOWNS (prob.nExtra): unknowns appended
%   after t_f that are not junction states, with their own equations
%   (prob.extraEq) and an extra Jacobian block in the terminal condition
%   (prob.terminal's third output). Needed for the homogeneous PMP
%   normalisation rho^2 + |lam_0|^2 = 1 (FINDINGS 33 correction): the
%   objective multiplier rho is exactly such an unknown -- it appears only in
%   the terminal H = rho + lam'f = 0 and in the normalisation equation.
%
%   Fixture: the engine's own harmonic-oscillator demo BVP, with one extra
%   unknown sigma and the equation sigma - 2 = 0, UNCOUPLED from the
%   dynamics. The solution must equal the nExtra = 0 solution exactly, with
%   sigma = 2; and the packing (t_f before the extras, extras last) must be
%   what the header promises.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

% the demo BVP: y'' = -y, y(0) = 0, y(tf) = 1, y'(tf) = 0 (free tf -> pi/2)
A  = [0 1; -1 0];
pr = struct('ny', 2, 'freeIdx0', 2, ...
    'prop', @(dt, y0, needSTM) demoProp(A, dt, y0, needSTM), ...
    'rhs',  @(y) A*y, ...
    'terminal', @(y, needJ) deal([y(1) - 1; y(2)], [1 0; 0 1]));
K = 4;  tg = linspace(0, 1.2, K+1);
Yg = [sin(tg); cos(tg)];
sd = struct('tf', 1.2, 'tGrid', tg, 'Y', Yg);
[p0, i0] = ms_bvp(pr, sd, struct('tolR', 1e-10));
ok = chk(ok, i0.converged && abs(p0(end) - pi/2) < 1e-7, ...
         sprintf('baseline: tf = %.9f (pi/2), converged = %d', p0(end), i0.converged));

% same problem + one extra unknown sigma with the equation sigma - 2 = 0
prx = pr;
prx.nExtra   = 1;
prx.terminal = @(y, needJ, x) deal([y(1) - 1; y(2)], [1 0; 0 1], zeros(2, 1));
prx.extraEq  = @(p1, x) deal(x - 2, zeros(1, numel(p1)), 1);
sdx = sd;  sdx.extra = 0.5;                      % seed for sigma
[px, ix] = ms_bvp(prx, sdx, struct('tolR', 1e-10));
ok = chk(ok, ix.converged, sprintf('with nExtra = 1: converged = %d, normR = %.1e', ix.converged, ix.normR));
ok = chk(ok, numel(px) == numel(p0) + 1, sprintf('packing: n = %d (baseline %d + 1)', numel(px), numel(p0)));
ok = chk(ok, abs(px(end) - 2) < 1e-9, sprintf('extra unknown solved: sigma = %.9f (want 2)', px(end)));
ok = chk(ok, abs(px(end-1) - pi/2) < 1e-7, sprintf('t_f sits BEFORE the extras: p(end-1) = %.9f', px(end-1)));
ok = chk(ok, norm(px(1:end-1) - p0) < 1e-8, ...
         sprintf('uncoupled extra leaves the BVP solution unchanged: |dp| = %.1e', norm(px(1:end-1) - p0)));

% assembleOnly must expose the extras too
[~, ia] = ms_bvp(prx, sdx, struct('assembleOnly', true));
ok = chk(ok, isfield(ia, 'nExtra') && ia.nExtra == 1 && numel(ia.p) == numel(p0) + 1, ...
         'assembleOnly reports nExtra and the full unknown vector');

if ok, fprintf('TEST_MS_BVP_EXTRA: ALL PASS\n');
else,  fprintf('TEST_MS_BVP_EXTRA: FAILURE (see lines above)\n');
end
end

function [yh, PHI] = demoProp(A, dt, y0, needSTM)
% DEMOPROP  Exact linear propagation with STM.  INPUTS: A; dt; y0; needSTM.
% OUTPUTS: yh; PHI.
PHI = expm(A*dt);  yh = PHI*y0;
if ~needSTM, PHI = []; end
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
