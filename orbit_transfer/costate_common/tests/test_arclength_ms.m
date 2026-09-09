function ok = test_arclength_ms()
%% Purpose:
%
%   Tests arclength_ms -- the GENERIC pseudo-arclength continuation engine
%   (any residual factory, any parameter) -- on an algebraic problem with a
%   KNOWN fold, so that the tangent, corrector, step control, fold
%   classification and grid-level crossing detection are each checked
%   against exact answers before the engine ever touches a transfer.
%
%   Fixture: R(x, q) = [x1^2 + x2^2 - 1;  x2 - q] with x in R^2, q the
%   parameter. The solution curve is the unit circle, x = (+-sqrt(1-q^2), q):
%   a FOLD in q at (0, 1) and (0, -1), where the fixed-q Jacobian
%   [2x1 2x2; 0 1] loses rank (x1 = 0) while the augmented
%   [2x1 2x2 0; 0 1 -1] keeps full row rank. Starting at (1, 0), q = 0,
%   heading toward increasing q, the arc must: reach q = 1 (tangent's
%   q-component -> 0 then changes sign), continue onto the x1 < 0 branch,
%   and record crossings of every requested level with x1 of the right
%   sign on each pass.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

resFactory = @(q) @(x) circleRes(x, q);
dRdq       = @(x, q) [0; -1];
p0 = [1; 0];  q0 = 0;
A = arclength_ms(resFactory, dRdq, p0, q0, struct( ...
    'direction', +1, 'ds', 0.1, 'dsMax', 0.3, 'nStep', 200, 'Dx', [1; 1], ...
    'qStop', [-0.3 2], 'levels', [0.5 0.9 -0.5], 'newtonTol', 1e-12, ...
    'logFile', ''));

% every accepted point is a root on the circle
r = cellfun(@(x) abs(x(1)^2 + x(2)^2 - 1), A.p);
ok = chk(ok, max(r) < 1e-10, sprintf('all %d accepted points are roots: max |R| = %.1e', numel(A.p), max(r)));

% it rounds the fold at q = 1: q rises to ~1 and comes back down
% accepted roots straddle the fold by up to ds; the fold itself is checked
% by its LOCALIZED position below
ok = chk(ok, max(A.q) > 0.98 && A.q(end) < -0.3, ...
         sprintf('rounded the fold: max accepted q = %.6f, final q = %.3f', max(A.q), A.q(end)));

% the fold is classified as a fold (tangent q-component sign change with
% rank-one loss of R_x and regular augmented matrix)
ok = chk(ok, numel(A.folds) == 1, sprintf('%d fold(s) detected (want exactly 1 before qStop)', numel(A.folds)));
if ~isempty(A.folds)
    kf = A.folds(1).index;
    ok = chk(ok, abs(A.folds(1).q - 1) < 1e-6, ...
             sprintf('fold LOCALIZED at q = %.8f (exact 1)', A.folds(1).q));
    ok = chk(ok, A.folds(1).classified && A.folds(1).sminX < 1e-2*A.folds(1).sminAug, ...
             sprintf('rank-one loss of R_x with regular augmented: smin %.1e vs %.1e', ...
                     A.folds(1).sminX, A.folds(1).sminAug));
end

% grid crossings: level 0.5 must be crossed TWICE (x1 > 0 going up, x1 < 0
% coming down), 0.9 twice, -0.5 not at all (never reached)
c = A.crossings;
n05 = sum(abs([c.level] - 0.5) < 1e-12);  n09 = sum(abs([c.level] - 0.9) < 1e-12);
nm5 = sum(abs([c.level] + 0.5) < 1e-12);
ok = chk(ok, n05 == 2 && n09 == 2 && nm5 == 0, ...
         sprintf('crossings: level 0.5 x%d, 0.9 x%d, -0.5 x%d (want 2, 2, 0)', n05, n09, nm5));
i05 = find(abs([c.level] - 0.5) < 1e-12);
if numel(i05) == 2
    x1s = sort(cellfun(@(pp) pp(1), {c(i05).p}));
    ok = chk(ok, abs(x1s(1) + sqrt(0.75)) < 1e-9 && abs(x1s(2) - sqrt(0.75)) < 1e-9, ...
             sprintf('level-0.5 crossings converged to +-sqrt(3)/2: %.9f, %.9f', x1s(2), x1s(1)));
    ok = chk(ok, all([c(i05).converged]), 'level crossings are re-solved at the exact level');
end

% budget/stop reasons are reported, not silent
ok = chk(ok, ischar(A.stop) && ~isempty(A.stop), sprintf('stop reason recorded: "%s"', A.stop));

% ---- a level crossed TWICE INSIDE ONE STEP ------------------------------
% Bracketing a level only between consecutive accepted roots misses a level
% the arc steps over and back within a single step -- which is exactly what
% happens near a fold, the place a sheet gains its second candidate. Fixed
% step 0.25 from angle 1.4 rad: the arc runs q = 0.9855 -> 1 (the fold) ->
% 0.9969, so level 0.999 is crossed twice with BOTH step endpoints below it.
a0 = 1.4;
A2 = arclength_ms(resFactory, dRdq, [cos(a0); sin(a0)], sin(a0), struct( ...
    'direction', +1, 'ds', 0.25, 'dsMin', 0.25, 'dsMax', 0.25, 'nStep', 3, ...
    'Dx', [1; 1], 'qStop', [-0.3 2], 'levels', 0.999, 'newtonTol', 1e-12));
c2 = A2.crossings;
qq = A2.q;
ok = chk(ok, all(qq < 0.999), sprintf('fixture: every accepted root stays below the level (max q = %.6f)', max(qq)));
ok = chk(ok, numel(c2) == 2, sprintf('level crossed twice inside one step: %d crossing(s) found (want 2)', numel(c2)));
if numel(c2) == 2
    x1s = sort(cellfun(@(v) v(1), {c2.p}));
    ok = chk(ok, all([c2.converged]) && abs(x1s(1) + x1s(2)) < 1e-8 && abs(abs(x1s(1)) - sqrt(1 - 0.999^2)) < 1e-8, ...
             sprintf('both roots of the double crossing recovered: x1 = %+.6f, %+.6f', x1s(1), x1s(2)));
end

if ok, fprintf('TEST_ARCLENGTH_MS: ALL PASS\n');
else,  fprintf('TEST_ARCLENGTH_MS: FAILURE (see lines above)\n');
end
end

function [R, J] = circleRes(x, q)
% CIRCLERES  Unit circle intersected with x2 = q.  INPUTS: x [2x1]; q.
% OUTPUTS: R [2x1]; J [2x2].
R = [x(1)^2 + x(2)^2 - 1; x(2) - q];
J = [2*x(1), 2*x(2); 0, 1];
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
