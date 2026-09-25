function chk = oscillator_checks(prob)
%% Purpose:
%
%   The oscillator checks as a FUNCTION OF THE PROBLEM STRUCT, so the same
%   checks can be run on the true problem (test_cov_oscillator: all pass)
%   and on a deliberately corrupted one (test_cov_mutation: the check aimed
%   at the corruption must fail). A check that cannot fail proves nothing.
%
%   prob must be the oscillator F = y'^2 - y^2 with a = 0, ya = 0, yb = 1
%   (the analytic answers below assume it); only its handles may be
%   mutated.
%
%% Inputs:
%
%  prob                     struct                  cov_problem output
%
%% Outputs:
%
%  chk                      struct                  one field per check,
%                                                   each .ok (logical) .msg
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

b = prob.b;
chk = struct();

% extremal by shooting vs sin(t)/sin(b)
E = cov_extremals(prob, [-12 12]);
if numel(E) ~= 1
    chk.extremal = res(false, sprintf('%d extremals found, expected 1', numel(E)));
    return
end
tt = linspace(0, b, 201);
z = deval(E.S.sol, tt);
err = max(abs(z(1,:) - sin(tt)/sin(b)));
chk.extremal = res(abs(E.p - 1/sin(b)) < 1e-8 && err < 1e-8, ...
                   sprintf('slope err %.1e, curve err %.1e', abs(E.p - 1/sin(b)), err));

% the functional vs cot(b)
chk.functional = res(abs(E.J - cot(b)) < 1e-8*max(1, abs(cot(b))), ...
                     sprintf('J = %.10f, cot(b) = %.10f', E.J, cot(b)));

% Jacobi field vs sin(t), and its zeros vs k pi
tEnd = 1.5*b;
S = cov_shoot(prob, E.p, tEnd);
tt = linspace(0, tEnd, 301);
z = deval(S.sol, tt);
hErr = max(abs(z(3,:) - sin(tt)));
kpi = pi*(1:floor(tEnd/pi));
zerosOk = numel(S.tConj) == numel(kpi) && all(abs(S.tConj - kpi) < 1e-8);
chk.jacobi = res(hErr < 1e-8 && zerosOk, ...
                 sprintf('max|h - sin t| = %.1e, zeros %s', hErr, mat2str(S.tConj, 10)));

% neighbours cross at k pi for a LARGE delta too (linear problem)
P = cov_perturbed(prob, E, 1.0, tEnd);
crossOk = numel(P.tCross) == numel(kpi) && all(abs(P.tCross - kpi) < 1e-7);
chk.neighbours = res(crossOk, sprintf('crossings %s', mat2str(P.tCross, 10)));

% second variation: eigenvalues and the Morse count
V = cov_second_variation(prob, E, 400);
lamTrue = 2*(((1:3)*pi/b).^2 - 1);
relErr = max(abs(V.lambda(1:3).' - lamTrue)./max(1, abs(lamTrue)));
nNegTrue = nnz((1:10)*pi < b);
chk.second_variation = res(relErr < 1e-3 && V.nNeg == nNegTrue, ...
    sprintf('lambda(1:3) rel err %.1e, nNeg %d (expected %d)', relErr, V.nNeg, nNegTrue));

% Delta J: exactly quadratic, and its sign follows the lowest eigenvalue
quadErr = max(abs(V.dJ - V.dJquad))/max(abs(V.dJquad));
signOk = all(sign(V.dJ) == sign(V.lambda(1)));
chk.delta_J = res(quadErr < 1e-6 && signOk, ...
                  sprintf('|dJ - quad|/|quad| = %.1e, sign(dJ) = sign(lambda1): %d', quadErr, signOk));
end

% ---------------------------------------------------------------------------
function r = res(ok, msg)
% RES  One check result. INPUTS: ok logical, msg char. OUTPUTS: r struct.
r = struct('ok', logical(ok), 'msg', msg);
end
