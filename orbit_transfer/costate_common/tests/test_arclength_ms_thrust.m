function test_arclength_ms_thrust()
% TEST_ARCLENGTH_MS_THRUST  Regression: the GENERIC engine arclength_ms,
% fed the production thrust residual (normal chart ms_tfmin) with the same
% scaling arclength_thrust used, must reproduce the 2026-09-08 diagnostic
% arc from the 75.5 mN root -- the t_f(T) curve to 1e-7 relative over the
% first millinewton, every root converged, no fold. The reference curve
% lives in fixtures/thrust_arc_75mN.mat (built once by
% fixtures/build_thrust_arc_fixture.m from the archived run).
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per check; errors on any FAIL)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
F = load(fullfile(here, 'fixtures', 'thrust_arc_75mN.mat'));

lStar = 389703.264829278;  tStar = 382981.289129055;
ndT = @(TN) (TN/F.m0kg)*tStar^2/(lStar*1000);
n = numel(F.p0);  K = (n - 8)/14 + 1;
Y1 = [F.rv0(:); 1; F.p0(1:7)];
seed0 = struct('tf', F.p0(end), 'tGrid', linspace(0, F.p0(end), K+1), ...
               'Y', [Y1, reshape(F.p0(8:end-1), 14, K-1), zeros(14,1)]);
seed0.Y(:, end) = seed0.Y(:, end-1);

mkRes = @(TN) resHandle(F.rv0, F.rvf, seed0, ndT(TN), F.cnd, F.muStar);
sT = F.TN0;
dRdT = @(p, TN) fdRT(mkRes, p, TN);
Dx = max(abs(F.p0), 1e-2);                  % exactly arclength_thrust's choice

T1 = tic;
A = arclength_ms(mkRes, dRdT, F.p0, F.TN0, struct( ...
    'direction', -1, 'Dx', Dx, 'sq', sT, 'ds', 0.05, 'dsMin', 1e-4, 'dsMax', 3.0, ...
    'nStep', 60, 'qStop', [F.TN0 - 1.0e-3, inf], 'newtonTol', 1e-9));
fprintf('  arc: %d roots, %.1f -> %.3f mN, stop = %s, %.0f s\n', ...
    numel(A.q), A.q(1)*1e3, A.q(end)*1e3, A.stop, toc(T1));

nFail = 0;
function chk(ok, msg)
    if ok, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); nFail = nFail + 1; end
end

chk(numel(A.q) >= 10 && A.q(end) < F.TN0 - 0.9e-3, ...
    sprintf('arc covered the first millinewton (%d roots, to %.3f mN)', numel(A.q), A.q(end)*1e3));
chk(max(A.normR) < 1e-8, sprintf('every root converged: max |R| = %.1e', max(A.normR)));
chk(isempty(A.folds), 'no fold on the monotone 75.5 -> 74.5 mN stretch');

% t_f(T) against the archived arc (both monotone in T here)
tf = cellfun(@(p) p(end), A.p);
msk = F.refT_N >= min(A.q) & F.refT_N <= max(A.q);
tfRef = interp1(A.q, tf, F.refT_N(msk), 'pchip');
relErr = max(abs(tfRef - F.reftf_nd(msk)) ./ F.reftf_nd(msk));
chk(relErr < 1e-7, sprintf('t_f(T) matches the archived arc: max rel err %.1e over %d reference roots', relErr, nnz(msk)));

% the costate norm trend (the chart running toward abnormality) is reproduced
lam0 = cellfun(@(p) norm(p(1:7)), A.p);
lamRef = interp1(F.refT_N(msk), F.reflam0(msk), A.q(end), 'pchip');
chk(abs(lam0(end) - lamRef)/lamRef < 1e-5, ...
    sprintf('|lam0| at the arc end matches: %.4f vs archived %.4f', lam0(end), lamRef));

if nFail > 0, error('TEST_ARCLENGTH_MS_THRUST: %d FAIL', nFail); end
fprintf('TEST_ARCLENGTH_MS_THRUST: ALL PASS\n');
end

function h = resHandle(rv0, rvf, seed, Tnd, cnd, mu)
% RESHANDLE  Production normal-chart ms residual for one thrust value.
% INPUTS: rv0; rvf; seed; Tnd; cnd; mu.  OUTPUTS: h (p -> [R, J]).
sd = seed;  sd.Y(1:7, 1) = [rv0(:); 1];
[~, inf_] = ms_tfmin(rv0(:), rvf(:), sd, Tnd, cnd, mu, ...
    struct('assembleOnly', true, 'handleOnly', true, 'tfLo', 0.05, 'tfHi', 20, 'pMax', 1e7));
h = inf_.residual;
end

function Rt = fdRT(mkRes, p, TN)
% FDRT  Central difference of the residual in thrust at fixed unknowns,
% h_rel = 1e-3 (the plateau found 2026-09-08).  INPUTS: mkRes; p; TN.
% OUTPUTS: Rt [n x 1] (UNSCALED d R / d T; the engine applies sq).
h = 1e-3*max(TN, 1e-3);
Rt = (feval(mkRes(TN + h), p) - feval(mkRes(TN - h), p))/(2*h);
end
