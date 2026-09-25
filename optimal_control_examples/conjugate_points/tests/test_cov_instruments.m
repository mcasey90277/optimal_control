function ok = test_cov_instruments()
%% Purpose:
%
%   Three statements that tie the instruments to each other on EVERY
%   extremal of EVERY preset (none of them is a closed form, but each pair
%   is computed by different code, so agreement is evidence):
%
%     1. MORSE: the number of negative second-variation eigenvalues equals
%        the number of zeros of the Jacobi field in (a, b).
%     2. The shooting function's slope at a root equals the Jacobi field at
%        b: r'(p) = h(b) (a central difference of r against the ODE value).
%     3. Delta J along the lowest mode has the sign of lambda_1, and at small
%        eps matches the quadratic prediction (eps^2/2) eta' K eta.
%
%   And the nonlinear limit, on the pendulum's second extremal: the
%   neighbour's crossing converges to the conjugate point as delta -> 0,
%   at FIRST order (the error halves when delta halves).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 all checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

addpath(fileparts(fileparts(mfilename('fullpath'))));
ok = true;
Pr = cov_presets;
for k = 1:numel(Pr)
    q = Pr(k);
    prob = cov_problem(q.F, q.a, q.b, q.ya, q.yb);
    E = cov_extremals(prob, q.pRange);
    ok = rep(ok, ~isempty(E), sprintf('[%d] extremal found', k), q.name);
    for e = 1:numel(E)
        tag = sprintf('[%d.%d]', k, e);
        V = cov_second_variation(prob, E(e), 400, 0.01*[-1 1]);
        nIn = nnz(E(e).S.tConj < q.b - 1e-9*(q.b - q.a));
        ok = rep(ok, V.nNeg == nIn, [tag ' Morse: nNeg = #conj in (a,b)'], ...
                 sprintf('%d vs %d', V.nNeg, nIn));

        dp = 1e-5*max(1, abs(E(e).p));
        Sp = cov_shoot(prob, E(e).p + dp, q.b);  Sm = cov_shoot(prob, E(e).p - dp, q.b);
        rPrime = (Sp.yEnd - Sm.yEnd)/(2*dp);
        ok = rep(ok, abs(rPrime - E(e).hb) < 1e-5*max(1, abs(E(e).hb)), ...
                 [tag ' r''(p) = h(b)'], sprintf('%.8f vs %.8f', rPrime, E(e).hb));

        relQ = max(abs(V.dJ - V.dJquad))/max(abs(V.dJquad));
        ok = rep(ok, all(sign(V.dJ) == sign(V.lambda(1))) && relQ < 0.05, ...
                 [tag ' Delta J: sign and quadratic'], ...
                 sprintf('lambda_1 %+.4f, dJ %s, rel %.1e', V.lambda(1), mat2str(V.dJ, 3), relQ));
    end
end

% the delta -> 0 limit on a nonlinear problem
q = Pr(strcmp({Pr.name}, 'Pendulum (nonlinear)'));
prob = cov_problem(q.F, q.a, q.b, q.ya, q.yb);
E = cov_extremals(prob, q.pRange);
[~, iw] = max(arrayfun(@(e) numel(e.S.tConj), E));      % the one with a conjugate point
S = cov_shoot(prob, E(iw).p, q.b);
tc = S.tConj(1);
d = 0.4*2.^-(0:4);
err = nan(size(d));
for k = 1:numel(d)
    P = cov_perturbed(prob, E(iw), d(k), q.b);
    if ~isempty(P.tCross), err(k) = abs(P.tCross(1) - tc); end
end
% first order is ASYMPTOTIC: at delta = 0.4 higher-order terms still show
% (ratio 1.67), so require decreasing errors, ratios rising toward 2, and
% the last two within 0.1 of 2
ratio = err(1:end-1)./err(2:end);
asym = all(diff(err) < 0) && all(diff(ratio) > 0) && all(abs(ratio(end-1:end) - 2) < 0.1);
ok = rep(ok, all(isfinite(err)) && asym, ...
         'pendulum: crossing -> t_c at first order', ...
         sprintf('t_c %.6f, errors %s, ratios %s', tc, mat2str(err, 3), mat2str(ratio, 3)));
fprintf('test_cov_instruments: %s\n', pf(ok));
end

% ---------------------------------------------------------------------------
function ok = rep(ok, c, name, msg)
% REP  Print one check and fold it into the verdict.
% INPUTS: ok running verdict, c this check, name, msg. OUTPUTS: ok.
fprintf('  %-44s %s   (%s)\n', name, pf(c), msg);
ok = ok && c;
end

function s = pf(ok)
% PF  PASS/FAIL text. INPUTS: ok logical. OUTPUTS: s char.
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
