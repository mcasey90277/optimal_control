function ok = test_sleigh_oracle()
%% Purpose:
%
%   The sleigh's hand-derived reduced equations against three independent
%   references: the symbolic Lagrange-d'Alembert builder evaluated on the
%   same solution (D1), the exact tanh solution (D2), energy conservation
%   and the constraint residual (D3, D4) -- for both signs of a, both signs
%   of v0, the knife edge a = 0, the straight ride w0 = 0 and rest. Also
%   the stability statement (D5) and the center-of-mass mapping.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));          % lagrangianHamiltonianMechanics
addpath(fileparts(here), fullfile(root, 'mechlib'));

odeOpts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
t = 0:0.01:8;
w0 = 2;

for a = [0.2, -0.2, 0]
    p = sleigh_params(1, 0.1, a);
    [q, qd, L, A] = sleigh_lagrangian(p);
    sys = mech.lagrange_dalembert(q, qd, L, A);
    for v0 = [1, -1]
        tag = sprintf('a = %+.1f, v0 = %+.0f', a, v0);
        [~, Yt] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; v0; w0], odeOpts);
        Y = Yt.';
        out = sleigh_quantities(t, Y, p);

        % D1: hand-derived derivative vs the builder, on the full state
        d1 = 0;
        for kk = 1:20:numel(t)
            d1 = max(d1, max(abs(sys.f(t(kk), out.Yfull(:, kk)) - out.Yfulld(:, kk))));
        end
        ok = check(ok, d1 <= 1e-9, sprintf('%s: D1 hand-derived vs builder %.2e (<= 1e-9)', tag, d1));

        % D2: exact solution
        d2v = max(abs(Y(4, :) - out.vExact));  d2w = max(abs(Y(5, :) - out.wExact));
        ok = check(ok, d2v <= 1e-8 && d2w <= 1e-8, sprintf('%s: D2 exact solution, v %.2e, w %.2e (<= 1e-8)', tag, d2v, d2w));

        % D3, D4
        chk = mech.conservation_check(t, out.Yfull, sys);
        ok = check(ok, chk.ok, sprintf('%s: D3 energy drift %.2e, D4 constraint %.2e (<= 1e-9) %s', tag, chk.driftRel, chk.resid, chk.why));

        % mapping: the center-of-mass velocity satisfies the constraint
        ok = check(ok, chk.resid <= 1e-12, sprintf('%s: center-of-mass velocity satisfies A qd = 0 (%.1e)', tag, chk.resid));

        % D5 and the sign structure
        ok = check(ok, all(sign(Y(5, :)) == sign(w0)), sprintf('%s: w never changes sign', tag));
        if a ~= 0
            ok = check(ok, abs(Y(4, end) - out.asymptote) < 1e-3, ...
                       sprintf('%s: v(end) = %.6f -> sign(a) sqrt(C) = %.6f (|diff| < 1e-3)', tag, Y(4, end), out.asymptote));
            vs = Y(4, Y(4, :) ~= 0);                  % ignore exact-zero samples: (-, 0, +) is ONE crossing
            nFlip = nnz(diff(sign(vs)) ~= 0);
            if a*v0 > 0
                ok = check(ok, nFlip == 0, sprintf('%s: center of mass leads, v keeps its sign (flips %d)', tag, nFlip));
            else
                ok = check(ok, nFlip == 1, sprintf('%s: center of mass trails, the signed speed crosses zero exactly once (flips %d)', tag, nFlip));
            end
            dth = abs(Y(3, end) - out.thetaInf);
            ok = check(ok, dth < 1e-4, sprintf('%s: total heading change %.6f vs exact sign(w0) (b/|a|) acos(sign(a) v0/S) = %.6f (|diff| %.1e < 1e-4)', ...
                       tag, Y(3, end) - Y(3, 1), out.thetaInf - Y(3, 1), dth));
        else
            cv = max(abs(Y(4, :) - v0));  cw = max(abs(Y(5, :) - w0));
            ok = check(ok, cv < 1e-10 && cw < 1e-10, sprintf('%s: knife edge, v and w constant (%.1e, %.1e)', tag, cv, cw));
            ok = check(ok, max(abs(out.lambda - p.m*Y(4, :).*Y(5, :))) < 1e-12, sprintf('%s: knife edge, lambda = m v w', tag));
            lamB = 0;
            for kk = 1:20:numel(t)
                lamB = max(lamB, abs(out.lambda(kk) - sys.lambda(out.Yfull(1:3, kk), out.Yfull(4:6, kk))));
            end
            ok = check(ok, lamB < 1e-12, sprintf('%s: knife edge, lambda vs builder multiplier %.2e (< 1e-12)', tag, lamB));
            rad = sqrt(sum((out.rb - [v0/w0*(-sin(0)); v0/w0]).^2, 1));   % circle centre (0, v0/w0) for the start at the origin heading east
            ok = check(ok, max(abs(rad - abs(v0/w0))) < 1e-9, sprintf('%s: knife edge runs a circle of radius |v0/w0| (%.1e)', tag, max(abs(rad - abs(v0/w0)))));
        end
    end
end

% w0 = 0: straight line at speed v0 for a ~= 0; rest when v0 = w0 = 0
p = sleigh_params(1, 0.1, 0.2);
[~, Ys] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; 1; 0], odeOpts);
outS = sleigh_quantities(t, Ys.', p);
ok = check(ok, max(abs(Ys(:, 2))) < 1e-12 && max(abs(Ys(:, 3))) < 1e-12 && max(abs(Ys(:, 4) - 1)) < 1e-12, ...
           'w0 = 0: straight line, th constant, v constant');
ok = check(ok, max(abs(outS.vExact - 1)) < 1e-15 && all(outS.wExact == 0) && all(isfinite(outS.vExact)), ...
           'w0 = 0: exact solution is the constant, not NaN');
outR = sleigh_quantities(t, zeros(5, numel(t)), p);
ok = check(ok, all(isfinite(outR.vExact)) && all(isfinite(outR.wExact)) && outR.C == 0 && all(outR.E == 0), ...
           'v0 = w0 = 0: C = 0, exact solution finite (rest), E = 0');

% trailing launch with w0 = 0 is an (unstable) equilibrium: it never reverses
[~, Yu] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), t, [0; 0; 0; -1; 0], odeOpts);
ok = check(ok, max(abs(Yu(:, 4) + 1)) < 1e-12 && max(abs(Yu(:, 5))) < 1e-12, ...
           'a = +0.2, v0 = -1, w0 = 0: the trailing straight ride is an exact (unstable) equilibrium, v stays -1');

% tiny w0: the exact solution must stay finite (the atanh form would give NaN because C rounds to v0^2)
tL = 0:0.05:25;
% AbsTol must sit far below w0 = 1e-10, else the integrator itself is the error
[~, Yt] = ode113(@(tt, yy) sleigh_rhs(tt, yy, p), tL, [0; 0; 0; -1; 1e-10], odeset('RelTol', 1e-12, 'AbsTol', 1e-24));
outT = sleigh_quantities(tL, Yt.', p);
ok = check(ok, all(isfinite(outT.vExact)) && all(isfinite(outT.wExact)) && abs(outT.wExact(1) - 1e-10) < 1e-23, ...
           'w0 = 1e-10: exact solution finite, wExact(1) = w0');
kz = find(Yt(:, 4) > 0, 1);
ok = check(ok, ~isempty(kz) && abs(tL(kz) - 17.29) < 0.1 && max(abs(Yt(:, 4).' - outT.vExact)) < 1e-6, ...
           sprintf('w0 = 1e-10: the trailing sleigh reverses near 17.29 s (first v > 0 at %.2f s), exact solution tracks it', tL(max(kz, 1))));

ok = check(ok, throws(@() sleigh_params(0, 0.1, 0.2), 'sleigh_params:value') && ...
               throws(@() sleigh_params(1, -1, 0.2), 'sleigh_params:value') && ...
               throws(@() sleigh_params(1, 0.1, NaN), 'sleigh_params:value'), ...
           'sleigh_params refuses m <= 0, J <= 0, non-finite a with :value');
p0 = sleigh_params();
ok = check(ok, p0.m == 1 && p0.J == 0.1 && p0.a == 0.2 && abs(p0.Jb - 0.14) < 1e-15 && abs(p0.k - 0.2/0.14) < 1e-15, ...
           'sleigh_params defaults and derived constants');

if ok, fprintf('TEST_SLEIGH_ORACLE: ALL PASS\n'); else, fprintf('TEST_SLEIGH_ORACLE: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end
