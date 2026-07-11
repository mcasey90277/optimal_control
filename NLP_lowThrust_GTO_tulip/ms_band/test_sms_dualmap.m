% TEST_SMS_DUALMAP  Adjudicate the dual->node-costate conversion candidates.
%
% For each SMS_SEED_DUALS mode ('a' baseline, 'b' h-weighted, 'c'
% adjacent-h averaged, 'd' midpoint-principled), evaluates along the FULL
% stored 1.12x direct trajectory (no MS solves — cheap):
%   (i)   one-arc propagation error at five arcs of the M = 40 joint grid
%         (early / perigee-crossing / switch-containing / mid / late),
%         split into state rows 1:8 and costate rows 9:16;
%   (ii)  max & rms |Ht + lamT| over all nodes, Ht computed from the
%         direct solution's own states/controls with the candidate
%         costates (the true costates satisfy Ht + lamT = 0);
%   (iii) adjoint-defect residual: nonuniform central-FD d(lam)/d(tau)
%         from the candidate node costates vs the SMS_EOM adjoint RHS
%         evaluated on the direct trajectory (rows lamR/lamV/lamM/lamT;
%         switch-layer nodes +-3 and the ends excluded), RMS + max.
% Prints a comparison table; error() only on harness failure. The winner
% call (orders-of-magnitude separation expected if a map is right) is for
% the campaign record.
%
% REFERENCES:
%   [1] .superpowers/sdd/gpt56_review_S1.md (validation prescription).
setup_paths;
matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
M = 40;
epsEval = 1e-3;
modes = {'a', 'b', 'c', 'd'};
T = struct('mode', {}, 'beta', {}, 'arcState', {}, 'arcCost', {}, ...
           'HtMax', {}, 'HtRms', {}, 'adjRmsR', {}, 'adjRmsV', {}, ...
           'adjRmsM', {}, 'adjRmsT', {}, 'adjMax', {});

for mi = 1:numel(modes)
    md = modes{mi};
    [Zseed, prob, info] = sms_seed_duals(matFile, M, epsEval, md);
    [~, yJ] = sms_unpack(Zseed, prob);
    tauN = info.tauN;  Y16 = info.Y16;  X = info.X;  U = info.U;
    nN = size(X, 2);
    s  = U(4, :);

    % ---- select arcs: early, perigee-crossing, switch, mid, late ----------
    r1 = sqrt(sum((X(1:3, :) - [-prob.muStar; 0; 0]).^2, 1));
    inner = tauN > 0.1*prob.sigf & tauN < 0.9*prob.sigf;
    rInner = r1;  rInner(~inner) = Inf;
    [~, kPeriNode] = min(rInner);
    swNodes = find(diff(double(s > 0.5)) ~= 0);
    [~, swPick] = min(abs(tauN(swNodes) - prob.sigf/2));
    arcOf  = @(nodeIdx) max(1, min(M, ...
             find(prob.sJ <= tauN(nodeIdx), 1, 'last')));
    arcSel = unique([2, arcOf(kPeriNode), arcOf(swNodes(swPick)), 20, M-1]);

    % ---- (i) one-arc propagation on the seed joints ------------------------
    aS = 0;  aC = 0;
    fprintf('mode %s: beta = %.5f  arcs [%s]\n', md, info.beta, ...
            sprintf('%d ', arcSel));
    for k = arcSel
        [~, Yarc] = ode113(@(ss, y) sms_eom(ss, y, prob.Tmax, prob.c, ...
                    prob.muStar, epsEval, prob.pSund), ...
                    [prob.sJ(k) prob.sJ(k+1)], yJ(:, k), prob.odeOpts);
        eS = max(abs(Yarc(end, 1:8).'  - yJ(1:8,  k+1)));
        eC = max(abs(Yarc(end, 9:16).' - yJ(9:16, k+1)));
        fprintf('  arc %2d: state err %.3e   costate err %.3e\n', k, eS, eC);
        aS = max(aS, eS);  aC = max(aC, eC);
    end

    % ---- (ii) Ht + lamT along the direct trajectory ------------------------
    HtErr = zeros(1, nN);
    for k = 1:nN
        HtErr(k) = ht_node(X(:, k), U(:, k), Y16(9:16, k), prob);
    end
    HtMax = max(abs(HtErr));  HtRms = rms(HtErr);

    % ---- (iii) adjoint-defect residual --------------------------------------
    keep = true(1, nN);
    keep([1:3, nN-2:nN]) = false;
    for w = -3:3, keep(max(1, min(nN, swNodes + w))) = false; end
    kIdx = find(keep);  kIdx = kIdx(kIdx >= 2 & kIdx <= nN-1);
    dFD  = (Y16(9:16, kIdx+1) - Y16(9:16, kIdx-1)) ...
           ./ (tauN(kIdx+1) - tauN(kIdx-1));
    dRHS = zeros(8, numel(kIdx));
    for q = 1:numel(kIdx)
        dY = sms_eom(0, Y16(:, kIdx(q)), prob.Tmax, prob.c, prob.muStar, ...
                     epsEval, prob.pSund);
        dRHS(:, q) = dY(9:16);
    end
    D = dFD - dRHS;
    T(end+1) = struct('mode', md, 'beta', info.beta, 'arcState', aS, ...
        'arcCost', aC, 'HtMax', HtMax, 'HtRms', HtRms, ...
        'adjRmsR', rms(D(1:3, :), 'all'), 'adjRmsV', rms(D(4:6, :), 'all'), ...
        'adjRmsM', rms(D(7, :)), 'adjRmsT', rms(D(8, :)), ...
        'adjMax', max(abs(D(:)))); %#ok<SAGROW>
end

fprintf('\n%-4s %-9s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
        'mode', 'beta', 'arcState', 'arcCost', 'Ht+lamT mx', 'Ht+lamT rms', ...
        'adjR rms', 'adjV rms', 'adjM rms', 'adjT rms', 'adj max');
for mi = 1:numel(T)
    fprintf('%-4s %-9.5f %-10.3e %-10.3e %-10.3e %-11.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n', ...
            T(mi).mode, T(mi).beta, T(mi).arcState, T(mi).arcCost, ...
            T(mi).HtMax, T(mi).HtRms, T(mi).adjRmsR, T(mi).adjRmsV, ...
            T(mi).adjRmsM, T(mi).adjRmsT, T(mi).adjMax);
end
score = [T.HtRms] + [T.adjRmsR] + [T.adjRmsV];   % dominant blocks
[~, win] = min(score);
fprintf('WINNER: mode %s (score %.3e vs baseline %.3e)\n', ...
        T(win).mode, score(win), score(1));
save('dualmap_table.mat', 'T', 'score', 'win');
fprintf('PASS test_sms_dualmap (table computed)\n');

% -------------------------------------------------------------------------
function e = ht_node(x, u, lam, prob)
% HT_NODE  |Ht + lamT| at one direct-solution node with candidate costates.
%
% Ht uses the DIRECT solution's own throttle/direction (hard running cost
% (Tmax/c)*u; entropy term ~0 on a converged near-bang solution).
%
% INPUTS:
%   x    - direct state [8x1] ([r;v;m;t])
%   u    - direct control [4x1] ([alpha;s])
%   lam  - candidate costates [8x1] ([lamR;lamV;lamM;lamT])
%   prob - problem struct (Tmax, c, muStar)
%
% OUTPUTS:
%   e - Ht + lamT [scalar]
r = x(1:3);  v = x(4:6);  m = x(7);
al = u(1:3);  s = u(4);
dd = [r(1) + prob.muStar;     r(2); r(3)];
rr = [r(1) - 1 + prob.muStar; r(2); r(3)];
gr = [r(1); r(2); 0] - (1 - prob.muStar)*dd./sqrt(sum(dd.^2))^3 ...
     - prob.muStar*rr./sqrt(sum(rr.^2))^3;
hv = [2*v(2); -2*v(1); 0];
Ht = (prob.Tmax/prob.c)*s + lam(1:3).'*v ...
     + lam(4:6).'*(gr + hv + s*prob.Tmax/m.*al) ...
     + lam(7)*(-s*prob.Tmax/prob.c);
e  = Ht + lam(8);
end
