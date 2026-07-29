% DIAG_ARCCHECK  Component-level breakdown of the seed_from_duals arc-check
% mismatch (arcCheckErr=0.35 at M=24, factor 1.12) before committing to the
% long LM solves.
setup_paths;
M = 24;
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';

S = load(matFile);
X = S.out.X;  U = S.out.U;  lamDef = S.out.lamDef;
p = cr3bp_lt_params(0.025, 15, 2100);
[beta, bInfo] = beta_from_duals(X, U, lamDef, p.c);
lamNode = -beta*[lamDef(1:7, :), lamDef(1:7, end)];
tSrc = X(8, :);
yNode = [X(1:7, :); lamNode];
tJsrc = arc_boundaries_tau(tSrc, X(1:3, :), M, p.muStar);
[tU, iu] = unique(tSrc);
ySrc = interp1(tU.', yNode(:, iu).', tJsrc.', 'linear').';

kMid = floor(M/2);
fprintf('kMid=%d  tJsrc(kMid)=%.6f  tJsrc(kMid+1)=%.6f  arc dt=%.6f\n', ...
    kMid, tJsrc(kMid), tJsrc(kMid+1), tJsrc(kMid+1)-tJsrc(kMid));

prob0 = ms_problem(S.factor, 1e-3);
[t1, Yarc] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob0.Tmax, prob0.c, ...
            prob0.muStar, prob0.epsSmooth), ...
            [tJsrc(kMid) tJsrc(kMid+1)], ySrc(:, kMid), prob0.odeOpts);

err7 = abs(Yarc(end, 1:7).' - ySrc(1:7, kMid+1));
labels = {'rx','ry','rz','vx','vy','vz','m'};
for k = 1:7
    fprintf('  %-3s  end=%12.6g  ref=%12.6g  |err|=%10.4g\n', ...
        labels{k}, Yarc(end,k), ySrc(k, kMid+1), err7(k));
end
fprintf('max err (rows1:7) = %.4e\n', max(err7));

% costate rows too, for context (not part of arcCheckErr but diagnostic)
errLam = abs(Yarc(end, 8:14).' - ySrc(8:14, kMid+1));
lamLabels = {'lr1','lr2','lr3','lv1','lv2','lv3','lm'};
for k = 1:7
    fprintf('  %-4s end=%12.6g  ref=%12.6g  |err|=%10.4g\n', ...
        lamLabels{k}, Yarc(end,7+k), ySrc(7+k, kMid+1), errLam(k));
end

% how many source-grid switches fall within this arc?
idxArc = tSrc >= tJsrc(kMid) & tSrc <= tJsrc(kMid+1);
sArc = U(4, idxArc);
nSwitch = nnz(diff(double(sArc > 0.5)) ~= 0);
fprintf('N source nodes in arc = %d,  throttle switches within arc = %d\n', ...
    nnz(idxArc), nSwitch);
fprintf('throttle range within arc: [%.4g, %.4g], mean=%.4g\n', ...
    min(sArc), max(sArc), mean(sArc));

% sub-divide the arc into 4 pieces and check divergence growth (chaos vs bug)
fprintf('\n--- sub-arc growth check (quarter arcs) ---\n');
tQ = linspace(tJsrc(kMid), tJsrc(kMid+1), 5);
yQ = interp1(tU.', yNode(:, iu).', tQ.', 'linear').';
y0 = ySrc(:, kMid);
for q = 1:4
    [~, Yq] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob0.Tmax, prob0.c, ...
                prob0.muStar, prob0.epsSmooth), ...
                [tQ(q) tQ(q+1)], y0, prob0.odeOpts);
    y0 = Yq(end, :).';
    errQ = max(abs(y0(1:7) - yQ(1:7, q+1)));
    fprintf('  after quarter %d (t=%.4f): max|err rows1:7| (propagated-from-scratch vs interp) = %.4e\n', ...
        q, tQ(q+1), errQ);
    % also reset to interpolated node each quarter (no error accumulation)
end

fprintf('\n--- reset-each-quarter (no compounding) ---\n');
for q = 1:4
    [~, Yq] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob0.Tmax, prob0.c, ...
                prob0.muStar, prob0.epsSmooth), ...
                [tQ(q) tQ(q+1)], yQ(:, q), prob0.odeOpts);
    errQ = max(abs(Yq(end, 1:7).' - yQ(1:7, q+1)));
    fprintf('  quarter %d alone (t=%.4f->%.4f): max|err rows1:7| = %.4e\n', ...
        q, tQ(q), tQ(q+1), errQ);
end
