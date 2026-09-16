function ok = test_newton_fixed_q()
%% Purpose:
%
%   Tests newton_fixed_q -- scaled Newton at a fixed continuation parameter
%   -- on a residual with a closed-form root:
%       R(p; q) = [p1^2 + p2^2 - q;  p1 - p2],   root p = sqrt(q/2) [1; 1].
%   Checks:
%     1. CONVERGES to the analytic root at the requested q, quadratically
%        (few calls), reporting |R|_inf < tol.
%     2. SCALING: a different Dx lands on the same root (Newton is affine
%        invariant, so the scaled step must not change the answer).
%     3. AT THE ROOT: one call, converged, p unchanged.
%     4. THE LEVEL IS HONOURED: the factory is asked for exactly q.
%     5. NON-FINITE RESIDUAL: stops and reports not converged.
%     6. NON-FINITE STEP (singular Jacobian): stops, not converged.
%     7. CAP: an iteration cap too small to converge reports not converged
%        with the final residual and nMax + 1 calls.
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
addpath(here);

q     = 3;
pStar = sqrt(q/2)*[1; 1];
fac   = @(qq) @(p) circleLine(p, qq);

%% 1. Converges to the analytic root:
[p, conv, nR, nC] = newton_fixed_q(fac, q, [2; 0.5], [1; 1], 1e-12, 20);
ok = chk(ok, conv && nR < 1e-12 && max(abs(p - pStar)) < 1e-12, ...
         sprintf('converged to sqrt(q/2)[1;1]: |p - p*| %.1e, |R| %.1e', max(abs(p - pStar)), nR));
ok = chk(ok, nC <= 8, sprintf('quadratically: %d residual calls', nC));

%% 2. Scaling does not change the root:
[p2, conv2] = newton_fixed_q(fac, q, [2; 0.5], [10; 0.1], 1e-12, 30);
ok = chk(ok, conv2 && max(abs(p2 - pStar)) < 1e-12, ...
         sprintf('Dx = [10; 0.1] lands on the same root (%.1e)', max(abs(p2 - pStar))));

%% 3. Starting at the root:
[p3, conv3, ~, nC3] = newton_fixed_q(fac, q, pStar, [1; 1], 1e-12, 20);
ok = chk(ok, conv3 && nC3 == 1 && isequal(p3, pStar), 'at the root: one call, converged, p unchanged');

%% 4. The level is honoured:
asked = [];
spy   = @(qq) recordLevel(qq);
[~, conv4] = newton_fixed_q(spy, 0.7, [1; 0.2], [1; 1], 1e-12, 20);
ok = chk(ok, conv4 && isequal(asked, 0.7), sprintf('the factory is asked for exactly q = %g', asked));

%% 5. Non-finite residual:
facNaN = @(qq) @(p) nanResidual(p);
[~, conv5, nR5, nC5] = newton_fixed_q(facNaN, q, [1; 1], [1; 1], 1e-12, 20);
ok = chk(ok, ~conv5 && nC5 == 1 && isnan(nR5), 'a NaN residual stops at once, not converged');

%% 6. Non-finite step:
facSing = @(qq) @(p) singularJacobian(p);
warnWas = warning('off', 'MATLAB:singularMatrix');
restore = onCleanup(@() warning(warnWas));
[~, conv6, ~, nC6] = newton_fixed_q(facSing, q, [1; 1], [1; 1], 1e-12, 20);
ok = chk(ok, ~conv6 && nC6 == 1, 'a singular Jacobian (non-finite step) stops, not converged');

%% 7. Iteration cap:
[~, conv7, nR7, nC7] = newton_fixed_q(fac, q, [50; -40], [1; 1], 1e-12, 2);
ok = chk(ok, ~conv7 && nC7 == 3 && nR7 > 1e-12, ...
         sprintf('cap of 2: not converged, |R| %.1e reported, %d calls', nR7, nC7));

if ok, fprintf('TEST_NEWTON_FIXED_Q: ALL PASS\n');
else,  fprintf('TEST_NEWTON_FIXED_Q: FAILURE (see lines above)\n');
end

    function h = recordLevel(qq)
    %% Purpose:
    %
    %   Residual factory that records the level it was asked for.
    %
    asked = qq;
    h = @(p) circleLine(p, qq);
    end
end

% ------------------------------------------------------------------------
function [R, J] = circleLine(p, q)
%% Purpose:
%
%   Circle of radius sqrt(q) intersected with the line p1 = p2.
%
R = [p(1)^2 + p(2)^2 - q;  p(1) - p(2)];
J = [2*p(1), 2*p(2);  1, -1];
end

function [R, J] = nanResidual(~)
%% Purpose:
%
%   A residual that has left the domain.
%
R = [NaN; 0];  J = eye(2);
end

function [R, J] = singularJacobian(~)
%% Purpose:
%
%   A nonzero residual with an all-zero Jacobian (non-finite Newton step).
%
R = [1; 1];  J = zeros(2);
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
