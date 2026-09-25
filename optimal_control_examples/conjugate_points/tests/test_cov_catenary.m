function ok = test_cov_catenary()
%% Purpose:
%
%   The minimal surface of revolution, F = y sqrt(1 + y'^2) on [-1/2, 1/2]
%   with y = 1 at both ends, against INDEPENDENT answers:
%
%     extremals   y = c cosh(t/c) with c cosh(1/(2c)) = 1: two roots,
%                 found here by fzero on that scalar equation alone
%     slopes      y'(-1/2) = -sinh(1/(2c))
%     functional  J = int c cosh^2(t/c) dt = c (1/2 + (c/2) sinh(1/c))
%     conjugate   LINDELOF'S CONSTRUCTION: the point conjugate to a on a
%     point       catenary is where the tangent there meets the t-axis at
%                 the same place as the tangent at a (a geometric fact,
%                 derived without any Jacobi equation)
%     verdicts    shallow catenary: no conjugate point, lambda_1 > 0,
%                 Delta J > 0; deep catenary: conjugate point inside,
%                 one negative mode, Delta J < 0
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 all checks passed
%
%% References:
%
%   [1] G. A. Bliss, "Calculus of Variations," Mathematical Association of
%       America, 1925 (the surface of revolution of minimum area,
%       Lindelof's construction).
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

addpath(fileparts(fileparts(mfilename('fullpath'))));
ok = true;
a = -0.5;  b = 0.5;
prob = cov_problem('y*sqrt(1+yp^2)', a, b, 1, 1);
E = cov_extremals(prob, [-8 2]);

% the oracle: the two catenary parameters
fc = @(c) c.*cosh(0.5./c) - 1;
cc = sort([fzero(fc, [0.1 0.4]), fzero(fc, [0.5 1.5])], 'descend');   % shallow, deep
pTrue = -sinh(0.5./cc);
JTrue = cc.*(0.5 + cc/2.*sinh(1./cc));
ok = rep(ok, numel(E) == 2, 'two extremals', sprintf('%d found', numel(E)));
if numel(E) ~= 2, fprintf('test_cov_catenary: FAIL\n'); return, end
% E is sorted by J: shallow (lower area) first
ok = rep(ok, all(abs([E.p] - pTrue) < 1e-8), 'slopes vs -sinh(1/2c)', ...
         sprintf('errors %s', mat2str(abs([E.p] - pTrue), 3)));
ok = rep(ok, all(abs([E.J] - JTrue) < 1e-8), 'functional vs closed form', ...
         sprintf('errors %s', mat2str(abs([E.J] - JTrue), 3)));

% Lindelof: tangent at a meets y = 0 at tStar; conjugate point tc has
% tc - y(tc)/y'(tc) = tStar on the same catenary (deep one: c = cc(2))
c = cc(2);
y  = @(t) c*cosh(t/c);   yp = @(t) sinh(t/c);
tStar = a - y(a)/yp(a);
tcTrue = fzero(@(t) t - y(t)/yp(t) - tStar, [0.01 0.49]);
S = cov_shoot(prob, E(2).p, b);
ok = rep(ok, numel(S.tConj) == 1 && abs(S.tConj - tcTrue) < 1e-7, ...
         'deep: conjugate point vs Lindelof', ...
         sprintf('tc = %s, Lindelof %.10f', mat2str(S.tConj, 10), tcTrue));
S1 = cov_shoot(prob, E(1).p, b);
ok = rep(ok, isempty(S1.tConj), 'shallow: no conjugate point in (a, b]', ...
         sprintf('zeros %s', mat2str(S1.tConj, 6)));

% second variation and Delta J
V1 = cov_second_variation(prob, E(1));
V2 = cov_second_variation(prob, E(2));
ok = rep(ok, V1.nNeg == 0 && V1.lambda(1) > 0 && all(V1.dJ > 0), ...
         'shallow: lambda_1 > 0, Delta J > 0', ...
         sprintf('lambda_1 %.4f, dJ %s', V1.lambda(1), mat2str(V1.dJ, 3)));
ok = rep(ok, V2.nNeg == 1 && all(V2.dJ < 0), 'deep: one negative mode, Delta J < 0', ...
         sprintf('nNeg %d, dJ %s', V2.nNeg, mat2str(V2.dJ, 3)));
ok = rep(ok, V1.Pmin > 0 && V2.Pmin > 0, 'Legendre: P > 0 on both', ...
         sprintf('min P %.3f, %.3f', V1.Pmin, V2.Pmin));
fprintf('test_cov_catenary: %s\n', pf(ok));
end

% ---------------------------------------------------------------------------
function ok = rep(ok, c, name, msg)
% REP  Print one check and fold it into the verdict.
% INPUTS: ok running verdict, c this check, name, msg. OUTPUTS: ok.
fprintf('  %-40s %s   (%s)\n', name, pf(c), msg);
ok = ok && c;
end

function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
