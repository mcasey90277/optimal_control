% TEST_PSR_MEE  No-solve-where-possible coverage for the MEE PSR port
% (psr_mee_refine.m + its two pure-function helpers, psr_switch_score_mee.m
% and psr_refine_sigma_mee.m). Mirrors test_warmstart_mee.m's structure:
% Parts 1-3 exercise the pure mesh-insertion logic on synthetic fixtures (no
% solve, no CasADi needed); Part 4 is ONE cheap budgeted live smoke of the
% full psr_mee_refine.m round loop against the certified 5 N fuel anchor
% (fast to re-solve, unlike the 1 N problem this port is validated against
% in the task-8 report).
%
% REFERENCES: [1] psr_switch_score_mee.m, psr_refine_sigma_mee.m (functions
%   under test in Parts 1-2). [2] interp_warmstart.m (handoff shape checked
%   in Part 3). [3] psr_mee_refine.m (Part 4). [4] test_warmstart_mee.m (the
%   no-certification-required live-smoke pattern Part 4 mirrors).

% =============================================================================
% Part 1: psr_switch_score_mee.m -- pure function, no solve
% =============================================================================

% Synthetic 20-interval (21-node) uniform sigma grid, throttle with two
% well-separated bang-bang switches at interval 5 (thr 1->0) and interval 14
% (thr 0->1), so the two switch windows (nbr=2 -> +-2 intervals) do not
% overlap (5+2=7 < 14-2=12).
N = 20;
sigma1 = linspace(0, 1, N + 1).';
thr = ones(1, N + 1);
thr(6:14) = 0;          % thr(1:5)=1, thr(6:14)=0, thr(15:21)=1
% diff(double(thr>0.5)) transitions: index 5 (1->0) and index 14 (0->1)
U1 = [zeros(3, N + 1); thr];   % beta rows unused by the scorer

optsSc = struct('nbr', 2);
[swIdx, score] = psr_switch_score_mee(sigma1, U1, optsSc);

assert(isequal(swIdx, [5 14]), 'test_psr_mee: expected switches at intervals [5 14], got %s', ...
    mat2str(swIdx));

% --- (a) score is exactly zero outside both +-2 windows -------------------
farIdx = [1 2 8 9 10 11 17 18 19 20];   % nowhere near interval 5 or 14 (windows: 3-7, 12-16)
assert(all(score(farIdx) == 0), ['test_psr_mee: score should be exactly zero far from any ' ...
    'switch, nonzero entries at %s'], mat2str(find(score(farIdx) ~= 0)));

% --- (b) score is strictly positive inside each window, peak at the switch
% interval itself, and decays monotonically with distance (taper) ----------
win1 = 3:7;   % interval 5 +- 2
assert(all(score(win1) > 0), 'test_psr_mee: score should be positive inside the switch-1 window');
[~, peakLoc] = max(score(win1));
assert(win1(peakLoc) == 5, 'test_psr_mee: score peak in window 1 should be AT interval 5');
% monotone non-increasing moving away from the peak on each side
left  = score(3:5);   right = score(5:7);
assert(all(diff(left) > 0),  'test_psr_mee: score should increase moving toward switch 5 from the left');
assert(all(diff(right) < 0), 'test_psr_mee: score should decrease moving away from switch 5 to the right');

fprintf(['test_psr_mee: Part 1 (psr_switch_score_mee.m, no solve) ALL PASS -- switches ' ...
    'detected at the right intervals, score zero outside local windows, peak at the switch, ' ...
    'tapered decay confirmed\n']);

% =============================================================================
% Part 2: psr_refine_sigma_mee.m -- pure function, no solve
% =============================================================================

% --- (a) basic bisection: sigma strictly increasing, no duplicates, all
% original nodes preserved, new nodes are interval midpoints -----------------
optsRef = struct('nbr', 2, 'K', Inf, 'hFloor', 1e-9, 'maxAdd', 2000);
[sigmaNew, isNew, nDropped] = psr_refine_sigma_mee(sigma1, score, optsRef);

assert(all(diff(sigmaNew) > 0), 'test_psr_mee: sigmaNew must be strictly increasing');
assert(numel(unique(sigmaNew)) == numel(sigmaNew), 'test_psr_mee: sigmaNew must have no duplicates');
for k = 1:numel(sigma1)
    assert(any(abs(sigmaNew - sigma1(k)) < 1e-14), ...
        'test_psr_mee: original node sigma1(%d)=%.6f missing from sigmaNew', k, sigma1(k));
end
assert(nDropped == 0, 'test_psr_mee: no interval should be dropped (all wide, all within maxAdd)');

% --- (b) insertion windows are LOCAL to the switches: every inserted node
% must be the midpoint of an interval within +-nbr of switch 5 or switch 14,
% i.e. sigmaNew(isNew) must all fall in [sigma1(3), sigma1(8)] U [sigma1(12), sigma1(17)]
% (interval k's midpoint lies in (sigma1(k), sigma1(k+1)); intervals 3-7 and 12-16 selected) --
insSig = sigmaNew(isNew);
inWin1 = insSig > sigma1(3) & insSig < sigma1(8);
inWin2 = insSig > sigma1(12) & insSig < sigma1(17);
assert(all(inWin1 | inWin2), ['test_psr_mee: all inserted nodes must be local to one of the ' ...
    'two switch windows -- found insertions outside both: %s'], mat2str(insSig(~(inWin1 | inWin2))));
assert(numel(insSig) == 10, 'test_psr_mee: expected 10 inserted nodes (5 intervals x 2 switches), got %d', ...
    numel(insSig));

fprintf(['test_psr_mee: Part 2a (psr_refine_sigma_mee.m basic bisection) ALL PASS -- sigma ' ...
    'strictly increasing, no duplicates, %d original nodes preserved, %d inserted nodes all ' ...
    'local to a switch window\n'], numel(sigma1), numel(insSig));

% --- (c) hFloor guard: an already-thin interval at a switch must NOT be
% bisected, and must be counted as dropped -----------------------------------
sigma2 = sigma1;
sigma2(6) = sigma2(5) + 1e-12;   % interval 5 now far below hFloor=1e-9
[~, score2] = psr_switch_score_mee(sigma2, U1, optsSc);
[sigmaNew2, isNew2, nDropped2] = psr_refine_sigma_mee(sigma2, score2, optsRef);
% interval 5 (now width 1e-12) must not have been bisected: no new node in (sigma2(5), sigma2(6))
midThin = (sigma2(5) < sigmaNew2) & (sigmaNew2 < sigma2(6));
assert(~any(midThin & isNew2), 'test_psr_mee: hFloor guard failed -- a sub-hFloor interval was bisected');
assert(nDropped2 >= 1, 'test_psr_mee: hFloor-thin interval should be counted as dropped (got nDropped=%d)', ...
    nDropped2);

% --- (d) maxAdd cap: capping to fewer insertions than the viable count must
% still produce a strictly increasing, duplicate-free grid, with the excess
% viable intervals counted as dropped ----------------------------------------
optsCap = struct('nbr', 2, 'K', Inf, 'hFloor', 1e-9, 'maxAdd', 3);
[sigmaNew3, isNew3, nDropped3] = psr_refine_sigma_mee(sigma1, score, optsCap);
assert(nnz(isNew3) == 3, 'test_psr_mee: maxAdd=3 should insert exactly 3 nodes, got %d', nnz(isNew3));
assert(nDropped3 == 10 - 3, 'test_psr_mee: maxAdd=3 should drop %d of the 10 viable intervals, got %d', ...
    10 - 3, nDropped3);
assert(all(diff(sigmaNew3) > 0), 'test_psr_mee: capped sigmaNew must still be strictly increasing');
assert(numel(unique(sigmaNew3)) == numel(sigmaNew3), 'test_psr_mee: capped sigmaNew must have no duplicates');
% the 3 kept insertions must be the 3 HIGHEST-scored intervals (closest to a switch, i.e. exactly
% intervals 5 and 14, plus one immediate neighbor -- score(5)==score(14) by symmetry, so the
% cap keeps both peaks plus the next-highest of the tied neighbors)
keptMid = sigmaNew3(isNew3);
assert(any(abs(keptMid - 0.5*(sigma1(5)+sigma1(6))) < 1e-14), ...
    'test_psr_mee: maxAdd cap should always keep the top-scored switch-5 interval');
assert(any(abs(keptMid - 0.5*(sigma1(14)+sigma1(15))) < 1e-14), ...
    'test_psr_mee: maxAdd cap should always keep the top-scored switch-14 interval');

fprintf(['test_psr_mee: Part 2b-c (hFloor guard, maxAdd cap) ALL PASS -- thin interval not ' ...
    'bisected and counted dropped; capped insertion keeps the top-scored intervals, grid stays ' ...
    'strictly increasing with no duplicates\n']);

% =============================================================================
% Part 3: handoff shapes into interp_warmstart.m (no solve)
% =============================================================================

% Synthetic source trajectory on sigma1 (21 nodes): distinct per-row state
% values, unit-norm rotating beta, bang-bang throttle matching U1 above.
Xsrc = zeros(7, N + 1);
for r = 1:7
    Xsrc(r, :) = 10*sigma1.' + r;
end
theta = linspace(0, 2*pi, N + 1);
betaSrc = [cos(theta); sin(theta); zeros(1, N + 1)];
Usrc = [betaSrc; thr];
dLsrc = 12.3;

W = interp_warmstart(Xsrc, Usrc, dLsrc, sigma1, sigmaNew);
assert(isequal(size(W.X), [7, numel(sigmaNew)]), ...
    'test_psr_mee: handoff W.X size mismatch, got %s, expected [7 %d]', ...
    mat2str(size(W.X)), numel(sigmaNew));
assert(isequal(size(W.U), [4, numel(sigmaNew)]), ...
    'test_psr_mee: handoff W.U size mismatch, got %s, expected [4 %d]', ...
    mat2str(size(W.U)), numel(sigmaNew));
betaNorms = sqrt(sum(W.U(1:3, :).^2, 1));
assert(max(abs(betaNorms - 1)) < 1e-10, ...
    'test_psr_mee: handoff beta rows must be unit-norm at every node (max dev %.3e)', ...
    max(abs(betaNorms - 1)));
srcThrSet = unique(thr);
for k = 1:numel(W.U(4, :))
    dmin = min(abs(W.U(4, k) - srcThrSet));
    assert(dmin < 1e-12, 'test_psr_mee: handoff throttle at node %d not a member of the bang-bang set', k);
end
assert(W.dL == dLsrc, 'test_psr_mee: handoff dL must pass through unchanged');

fprintf(['test_psr_mee: Part 3 (psr_refine_sigma_mee.m -> interp_warmstart.m handoff) ALL PASS ' ...
    '-- output shapes match the refined grid, beta renormalized, throttle bang-bang preserved, ' ...
    'dL passthrough correct\n']);

% =============================================================================
% Part 4: ONE cheap live smoke of psr_mee_refine.m's full round loop, budgeted
% (mirrors test_warmstart_mee.m Part 2's non-convergence-gated pattern, but
% here maxIter is generous enough that certification is plausible, since
% psr_mee_refine's round-accept/cache logic is only exercised by an actual
% certified resolve) -- runs against the certified 5 N fuel anchor (small,
% fast to re-solve), NOT the 1 N problem this port is validated against in
% the task-8 report (that lives in a separate, longer-running validation run).
% =============================================================================
here    = fileparts(mfilename('fullpath'));
resDir  = fullfile(here, 'results');
srcFile = fullfile(resDir, 'MEE_M2_5N.mat');
assert(isfile(srcFile), ['test_psr_mee: prerequisite %s not found -- this is the certified ' ...
    '25/rev 5 N fuel anchor Part 4 smoke-tests PSR against; re-run ' ...
    'run_transfer_mee(struct(''thrustN'',10)) to regenerate it first'], srcFile);
S = load(srcFile);
baseResult = S.res;

smokeTag = 'TEST_psr_smoke_T5N';
% clear any stale cache from a prior interrupted run of this exact smoke test
for f = {sprintf('%s_psr_history.mat', smokeTag), sprintf('%s_psr_r1.mat', smokeTag), ...
         sprintf('%s_psr_final.mat', smokeTag)}
    ff = fullfile(resDir, f{1});
    if isfile(ff), delete(ff); end
end

optsPsr = struct('maxRounds', 1, 'tag', smokeTag, 'maxIter', 1500, 'nbr', 2);

tStart = tic;
out = psr_mee_refine(baseResult, optsPsr);
wallSec = toc(tStart);

budget_s = 180;
assert(wallSec < budget_s, 'test_psr_mee: Part 4 live smoke took %.1fs, budget %.0fs', wallSec, budget_s);

assert(isfield(out, 'history') && numel(out.history) >= 1, ...
    'test_psr_mee: psr_mee_refine output missing a populated .history');
assert(isfield(out, 'finalOut') && isfield(out, 'finalSigma'), ...
    'test_psr_mee: psr_mee_refine output missing .finalOut/.finalSigma');
assert(out.history(1).nNodes == numel(baseResult.sigma) - 1, ...
    'test_psr_mee: round-0 history nNodes should match the seed''s own node count');
assert(out.history(1).switches == baseResult.fuel.switches, ...
    ['test_psr_mee: round-0 history switch count should match the seed''s own out.switches ' ...
     '(psr_switch_score_mee''s bracketing must agree with casadi_lt_mee''s)']);

fprintf(['test_psr_mee: Part 4 (psr_mee_refine.m, live, 5 N anchor) ALL PASS -- wallSec=%.2f, ' ...
    'rounds measured=%d, stopReason=%s, final certified=%d, final mf=%.4f kg, final sw=%d, ' ...
    'no exception\n'], wallSec, numel(out.history), out.stopReason, out.certified, ...
    out.finalOut.m_f_kg, out.finalOut.switches);

fprintf('test_psr_mee: ALL PASS (Parts 1-4)\n');
