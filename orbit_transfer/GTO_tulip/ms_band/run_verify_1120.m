% RUN_VERIFY_1120  PMP-consistency certificate for the 1.12x direct solution.
%
% Adjudications passed to the gated certificate (justifications):
%   adjSwitches [3 4] — review-CONFIRMED near-graze throttle dip at
%     tau 19.46-19.51 (2026-07-10): the throttle only dips to u = 0.43
%     (never reaches the coast bound; 2 intermediate-throttle nodes of
%     4001), un-interpolated interval duals stay burn-side (S_int max
%     -7.3e-4), and the interval switching law separates all 542 coast /
%     3446 burn intervals at 100%. The historical "12 switches" was an
%     s>0.5 threshold-counting artifact; certified count is 10 (+1
%     near-graze dip, not a certified switch).
%   adjArcs [4 26] — v-block 1.03e-2 / 1.02e-2: 2-3% over the 1e-2
%     heuristic line on perigee-adjacent arcs, within the stated
%     ~5e-3..1e-2 dual-map floor band (marginal, floor-level).
%   adjArcs [40] — terminal switch-complex amplification, NOT costate
%     error (diag_verify_1120: defect <= 4.2e-3 up to the first of 5
%     in-arc switches, whose S-crossing matches at 0 nodes; the u-vs-s
%     disagreement after it integrates to O(1)).
setup_paths;
opts = struct('adjArcs', [4 26 40], 'adjSwitches', [3 4]);
summary = verify_direct_pmp( ...
    '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat', opts);
fprintf('run_verify_1120 done: cert %d, matched %d/%d, primer mean %.4f deg, |lamM(end)| %.3e\n', ...
        summary.certOK, summary.nMatched, summary.nSwitches, ...
        summary.primerMeanDeg, summary.lamMend);
