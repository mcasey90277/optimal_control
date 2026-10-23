function ok = test_barrier_terms()
%% Purpose:
%
%   The paper's C1 log barrier and its assembled derivatives:
%     1. barrier_shape is C1 at h = l: U and U' continuous, U'' continuous
%        too (both branches give eps/l^2 there), finite at h = 0 and h < 0;
%     2. U' and U'' vs complex step and central differences of U;
%     3. barrier_terms' gradient vs a complex step of U (complex-step safe);
%     4. barrier_terms' Hessian vs a complex step of the gradient, and
%        symmetric;
%     5. the ball term enters the i-block only;
%     6. M = 1: U = 0, zero gradient, empty-sparse Hessian;
%     7. U equals a plain sum over barrier_shape(h) with h from
%        clearance_terms (the assembly agrees with the pieces).
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
rngWas = rng(3);  restore = onCleanup(@() rng(rngWas));
epsB = 1e-5;  l = 1e-4;

%% 1. the shape
[Ua, U1a, U2a] = barrier_shape(l*(1 + 1e-12), epsB, l);
[Ub, U1b, U2b] = barrier_shape(l*(1 - 1e-12), epsB, l);
% the two probes are 2e-12 l apart on a slope eps/l, so U itself differs by 2e-12 eps; the C1 statement is the U' match below
ok = chk(ok, abs(Ua - Ub) < 1e-11*epsB && abs(U1a - U1b) < 1e-9*epsB/l && abs(U2a - U2b) < 1e-9*epsB/l^2, ...
         'barrier_shape is C1 at h = l (U, U'' and U'''' all continuous)');
[U0, U10, U20] = barrier_shape([0, -0.01], epsB, l);
ok = chk(ok, all(isfinite([U0 U10 U20])), 'finite at h = 0 and h < 0 (quadratic branch)');
ok = chk(ok, abs(barrier_shape(1, epsB, l)) < 1e-20 && abs(barrier_shape(0.5, epsB, l) + epsB*log(0.5)) < 1e-20, ...
         'U = -eps log h on the log branch');

%% 2. shape derivatives
hs = [2e-4, 1e-3, 0.05, 0.4, 5e-5, -2e-3];
hc = 1e-20;  hd = 1e-7;
worst1 = 0;  worst2 = 0;
for k = 1:numel(hs)
    [~, U1, U2] = barrier_shape(hs(k), epsB, l);
    Uc = barrier_shape(hs(k) + 1i*hc, epsB, l);
    worst1 = max(worst1, abs(imag(Uc)/hc - U1)/max(1, abs(U1)));
    [~, U1p] = barrier_shape(hs(k) + hd, epsB, l);
    [~, U1m] = barrier_shape(hs(k) - hd, epsB, l);
    worst2 = max(worst2, abs((U1p - U1m)/(2*hd) - U2)/max(1, abs(U2)));
end
ok = chk(ok, worst1 < 1e-13, sprintf('U'' vs complex step (%.1e rel)', worst1));
ok = chk(ok, worst2 < 1e-6, sprintf('U'''' vs central differences of U'' (%.1e rel)', worst2));

%% 3-4. assembled gradient and Hessian, with a ball
[x0, ~, geom] = drone_layouts('sphere', 5);
M = 5;  pairs = pair_list(M);
x = x0(:) + 0.05*randn(3*M, 1);
[U, gU, HU] = barrier_terms(x, pairs, geom);
ok = chk(ok, isscalar(U) && isequal(size(gU), [3*M 1]) && isequal(size(HU), [3*M 3*M]) && issparse(HU), 'shapes; HU sparse');
worstG = 0;  worstH = 0;
for k = 1:3*M
    xc = x;  xc(k) = xc(k) + 1i*hc;
    [Uc, gc] = barrier_terms(xc, pairs, geom);
    worstG = max(worstG, abs(imag(Uc)/hc - gU(k))/max(1, abs(gU(k))));
    worstH = max(worstH, max(abs(imag(gc)/hc - HU(:,k)))/max(1, max(abs(HU(:,k)))));
end
ok = chk(ok, worstG < 1e-12, sprintf('gradient vs complex step of U (%.1e rel)', worstG));
ok = chk(ok, worstH < 1e-10, sprintf('Hessian vs complex step of the gradient (%.1e rel)', worstH));
ok = chk(ok, max(max(abs(HU - HU.'))) < 1e-14*max(max(abs(HU))), 'Hessian symmetric');

%% 5. the ball enters only the diagonal agent block
geomNB = geom;  geomNB.ball = [];
[~, gNB, HNB] = barrier_terms(x, pairs, geomNB);
D = HU - HNB;
offBlock = 0;
for a = 1:M
    for b = 1:M
        if a ~= b, offBlock = max(offBlock, max(max(abs(D(3*(a-1)+(1:3), 3*(b-1)+(1:3)))))); end
    end
end
ok = chk(ok, offBlock == 0 && max(abs(gU - gNB)) > 0, 'ball term: i-block only, and it is not zero');

%% 6. M = 1
[~, ~, g1] = drone_layouts('circle', 1);
[U1, g1v, H1] = barrier_terms(randn(3,1), pair_list(1), g1);
ok = chk(ok, U1 == 0 && isequal(g1v, zeros(3,1)) && issparse(H1) && nnz(H1) == 0, 'M = 1: U = 0, zero gradient, empty Hessian');

%% 7. assembly agrees with the pieces
[h, ~, ~, hB] = clearance_terms(x, pairs, geom);
Uplain = sum(barrier_shape(h, geom.eps, geom.l)) + sum(barrier_shape(hB, geom.eps, geom.l));
ok = chk(ok, abs(U - Uplain) < 1e-15*max(1, abs(U)), 'U = plain sum of barrier_shape over clearance_terms');

if ok, fprintf('TEST_BARRIER_TERMS: ALL PASS\n'); else, fprintf('TEST_BARRIER_TERMS: FAILURE\n'); end
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
