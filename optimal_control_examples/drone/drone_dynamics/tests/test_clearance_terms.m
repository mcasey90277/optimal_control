function ok = test_clearance_terms()
%% Purpose:
%
%   clearance_terms against a plain loop written here; M = 1 and N > 1
%   shapes; complex-step safety (the imaginary part of h under a complex
%   perturbation equals the analytic directional derivative n . dx);
%   the ball term. hermite_resample reproduces a cubic exactly and returns
%   the node values at the node positions.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(5);  restore = onCleanup(@() rng(rngWas));

[~, ~, geom] = drone_layouts('sphere', 6);
M = 6;  N = 4;  pairs = pair_list(M);
X = randn(3*M, N);
[h, rho, n, hB] = clearance_terms(X, pairs, geom);
ok = chk(ok, isequal(size(h), [15 N]) && isequal(size(rho), [15 N]) && isequal(size(n), [3 15 N]) && isequal(size(hB), [M N]), ...
         'shapes: h, rho [P x N], n [3 x P x N], hBall [M x N]');
worst = 0;  worstB = 0;
for c = 1:N
    P = reshape(X(:,c), 3, M);
    for q = 1:15
        a = pairs(1,q);  b = pairs(2,q);
        dd = sqrt(sum((P(:,a) - P(:,b)).^2));
        worst = max(worst, abs(h(q,c) - (dd - 2*geom.r)) + abs(rho(q,c) - dd) + max(abs(n(:,q,c) - (P(:,a) - P(:,b))/dd)));
    end
    for a = 1:M
        worstB = max(worstB, abs(hB(a,c) - (sqrt(sum((P(:,a) - geom.ball.c).^2)) - geom.ball.R)));
    end
end
ok = chk(ok, worst < 1e-14, sprintf('pair terms vs a plain loop (%.1e)', worst));
ok = chk(ok, worstB < 1e-14, sprintf('ball terms vs a plain loop (%.1e)', worstB));

% complex step: d h_q / d x along a random direction equals n . (dx_a - dx_b)
dx = randn(3*M, 1);  hc = 1e-20;
hcs = clearance_terms(X(:,1) + 1i*hc*dx, pairs, geom);
worst = 0;
for q = 1:15
    a = pairs(1,q);  b = pairs(2,q);
    dxa = dx(3*(a-1)+(1:3));  dxb = dx(3*(b-1)+(1:3));
    worst = max(worst, abs(imag(hcs(q))/hc - n(:,q,1).'*(dxa - dxb)));
end
ok = chk(ok, worst < 1e-13, sprintf('complex step of h matches n . dx (%.1e)', worst));

[~, ~, g1] = drone_layouts('circle', 1);
[h1, ~, ~, hB1] = clearance_terms(randn(3, 3), pair_list(1), g1);
ok = chk(ok, isequal(size(h1), [0 3]) && isempty(hB1), 'M = 1: h is [0 x N], no ball term');

% hermite_resample: a cubic is reproduced exactly
tN = linspace(0, 2, 6);
c3 = randn(2, 4);
poly  = @(t) c3(:,1) + c3(:,2)*t + c3(:,3)*t.^2 + c3(:,4)*t.^3;
dpoly = @(t) c3(:,2) + 2*c3(:,3)*t + 3*c3(:,4)*t.^2;
[Xf, Vf, tf] = hermite_resample(poly(tN), dpoly(tN), tN, 10);
ok = chk(ok, numel(tf) == 51 && isequal(size(Xf), [2 51]) && isequal(size(Vf), [2 51]), 'hermite_resample: N nSub + 1 points');
ok = chk(ok, max(max(abs(Xf - poly(tf)))) < 1e-12 && max(max(abs(Vf - dpoly(tf)))) < 1e-11, 'hermite_resample reproduces a cubic and its derivative');
ok = chk(ok, max(abs(tf(1:10:end) - tN)) < 1e-15 && max(max(abs(Xf(:,1:10:end) - poly(tN)))) < 1e-15, 'nodes are returned exactly');

if ok, fprintf('TEST_CLEARANCE_TERMS: ALL PASS\n'); else, fprintf('TEST_CLEARANCE_TERMS: FAILURE\n'); end
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
