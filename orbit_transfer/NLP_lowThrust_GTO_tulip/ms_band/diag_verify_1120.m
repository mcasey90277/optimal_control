% DIAG_VERIFY_1120  Dig: the 12-vs-10 switch-structure row + the arc-40 blowup.
%
% For each direct throttle switch of legacy_ms_f1120: location (node, tau,
% t), direction, dual-S value at the switch, min |S| within +-20 nodes,
% and distance to the nearest dual-S zero crossing — distinguishes a
% GRAZING S (resolution/noise-limited crossing: |S|_min ~ dual noise) from
% a real inconsistency (S far from 0 where the direct control switches).
% Also inspects the final arc (arc 40 of M = 40): S and throttle along the
% PROPAGATED arc vs the direct solution, and where the state defect
% accumulates.
setup_paths;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[~, prob, info] = sms_seed_duals(matFile, 40, 1e-4, 'd');
tauN = info.tauN;  Y16 = info.Y16;  X = info.X;  U = info.U;
nN = size(X, 2);
s  = U(4, :);

Snode = 1 - sqrt(sum(Y16(12:14, :).^2, 1))*prob.c./X(7, :) - Y16(15, :);
crossI = find(diff(sign(Snode)) ~= 0);
swI    = find(diff(double(s > 0.5)) ~= 0);

fprintf('=== switch-by-switch (12 direct vs %d dual-S crossings) ===\n', numel(crossI));
fprintf('%-3s %-6s %-9s %-8s %-6s %-10s %-10s %-9s\n', ...
        '#', 'node', 'tau', 't', 'dir', 'S@switch', 'min|S|+-20', 'distNode');
for q = 1:numel(swI)
    k = swI(q);
    dirStr = 'on->off';  if s(k+1) > 0.5, dirStr = 'off->on'; end
    win = max(1, k-20):min(nN, k+20);
    dN  = min(abs(crossI - k));
    fprintf('%-3d %-6d %-9.3f %-8.4f %-6s %-10.3e %-10.3e %-9d\n', ...
            q, k, tauN(k), X(8, k), dirStr, Snode(k), min(abs(Snode(win))), dN);
end
fprintf('crossing taus: %s\n', sprintf('%.3f ', ...
        interp1(1:nN, tauN, crossI + 0.5)));
fprintf('switch taus  : %s\n', sprintf('%.3f ', (tauN(swI) + tauN(swI+1))/2));

% burn/coast segment durations around each switch (tau units)
segB = diff([1, swI, nN]);
fprintf('segment lengths (nodes): %s\n', sprintf('%d ', segB));

% ---- arc 40 ----------------------------------------------------------------
yJ = interp1(tauN.', Y16.', prob.sJ.', 'linear').';
k  = 40;
[sk, Yk] = ode113(@(ss, y) sms_eom(ss, y, prob.Tmax, prob.c, prob.muStar, ...
           1e-4, prob.pSund), [prob.sJ(k) prob.sJ(k+1)], yJ(:, k), prob.odeOpts);
Sarc = 1 - sqrt(sum(Yk(:, 12:14).^2, 2)).'*prob.c./Yk(:, 7).' - Yk(:, 15).';
uArc = (1 - tanh(Sarc/(2*1e-4)))/2;
sDir = interp1(tauN, s, sk.', 'nearest');
r1a  = sqrt(sum((Yk(:, 1:3).' - [-prob.muStar; 0; 0]).^2, 1));
fprintf('\n=== arc 40 [%0.2f, %0.2f] ===\n', prob.sJ(40), prob.sJ(41));
fprintf('min r1 along arc: %.4f;  S range [%.3e, %.3e]\n', min(r1a), ...
        min(Sarc), max(Sarc));
% where do propagated u and direct s disagree?
dis = abs(uArc - sDir) > 0.5;
if any(dis)
    fprintf('u vs direct s disagree on %.1f%% of arc samples; first at tau %.3f (t %.4f)\n', ...
            100*mean(dis), sk(find(dis, 1)), Yk(find(dis, 1), 8));
else
    fprintf('u matches direct s along the whole arc\n');
end
% defect growth: compare propagated state to direct at several taus
tauChk = linspace(prob.sJ(40), prob.sJ(41), 9);
Xchk   = interp1(tauN.', Y16(1:8, :).', tauChk.').';
Ychk   = interp1(sk, Yk(:, 1:8), tauChk.').';
dChk   = max(abs(Ychk - Xchk), [], 1);
fprintf('arc-40 state defect growth at 9 stations: %s\n', sprintf('%.2e ', dChk));
fprintf('S at last 5 direct nodes: %s  (lamM(end) = %.3e)\n', ...
        sprintf('%.3e ', Snode(end-4:end)), Y16(15, end));
