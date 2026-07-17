function out = run_mintime_mee(thrustN, nodesPerRev, cfg)
% RUN_MINTIME_MEE  MEE/L-domain min-time anchor at one thrust level (Task 6).
%
% MEE analog of run_mintime.m: same two-stage-plus-continuation SPIRIT, but
% there is no rendezvous/manifold distinction here (the MEE terminal set is
% ALREADY the 5 fixed element constraints baked into casadi_lt_mee.m; total
% longitude span DeltaL is the free DOF, playing the role tau_f/L_end played
% on the Cartesian side). The two MEE-specific stages are:
%
%   Stage A - warm-start 'mintime' directly from the CERTIFIED fuel solution
%             at this thrust (results/<fuelTag>.mat, e.g. run_transfer_mee's
%             MEE_M2_10N.mat), throttle row forced to 1 (thr==1 is a hard
%             equality in 'mintime' mode regardless -- this just removes an
%             avoidable primal infeasibility at iteration 0, same lesson as
%             casadi_lt_mee.m's own U0w priming). The fuel solution's shape
%             is already a converged, defect-free trajectory satisfying the
%             SAME terminal set at a nearby (longer) transfer time, so it is
%             a far better warm start than a raw cold spiral -- only
%             available when a certified fuel anchor already exists for this
%             exact thrustN (true today only at 10 N; every lower thrust
%             falls through to Stage B until run_ladder.m produces one).
%   Stage B - a FRESH cold tangential seed (mee_seed, thr=cfg.seedThrB
%             [default 0.4 -- see STAGE-B SEED THROTTLE FIX below],
%             betaMode='tangential', nRev=cfg.nRevSeed [default 3, an
%             arbitrary but reasonable starting shape; DeltaL is free, so
%             the solve is not committed to nRevSeed revs, only seeded near
%             it]) at the same node density, then the same continuation
%             recipe.
%
% STAGE-B SEED THROTTLE FIX (2026-07-17, Task 7's first live ladder descent):
% the raw seed's throttle was originally hardcoded thr=1 (matching mintime's
% own all-burn physics). This is harmless at nRevSeed=3 (10 N, Task 6) but
% CATASTROPHIC once run_ladder.m's C-law warm hint pushes nRevSeed higher at
% lower thrust (nRevSeed=9 at the 5 N rung): a raw thr=1 constant-throttle
% ode113 integration is EXACTLY the seed shape run_transfer_mee.m's own
% DEVIATION note already warns about ("crosses GEO (P=1) at rev~3.06 and
% goes coordinate-singular past ~rev 4") -- confirmed live here too (offline
% probe: thr=1/nRev=9 at 5 N hits an ode113 tolerance failure with P blowing
% up to 1.3e6 by the end of the achieved partial integration), and feeding
% that garbage tail into casadi_lt_mee as a warm start produced instant
% NaN-in-constraint IPOPT failures and a hard stall (defect ~2e4, unmoving
% across 2 continuation rounds). FIX: Stage B's seed throttle now defaults
% to cfg.seedThrB=0.4 -- THE SAME constant-throttle value already proven
% robust for the FUEL seed (run_transfer_mee.m's seedThr) -- which stays
% cleanly conditioned across the WHOLE nRevSeed range the ladder needs
% (offline probe: thr=0.4/nRev=9 at 5 N: P range [0.276, 0.477], ex range
% [0.668, 0.750], no ode113 warning). This is a SEED-SHAPE-ONLY change: the
% actual mintime physics is unaffected, since casadi_lt_mee.m already forces
% the warm-start throttle row to 1 for 'mintime' mode regardless of what the
% seed's own U0 throttle contains (the "U0w(4,:)=1" priming line) -- thr=0.4
% only slows the RAW ballistic seed integration enough to stay well short of
% the coordinate singularity while still spanning the needed rev range.
% BASIN-MULTIPLICITY FINDING (2026-07-17, this task's own Step-2 validation
% run) -- CONTROLLER-AUTHORIZED DEVIATION from the brief's literal "Stage B
% only if Stage A fails" wording: at 10 N, Stage A (warm-started from the
% c_tf=1.5 fuel solution's ~7.3-rev shape, thr forced to 1) CONVERGES
% CLEANLY (Solve_Succeeded, defect 3.4e-15, 0 continuation rounds needed)
% but lands in a DIFFERENT local optimum -- 7.156 revs, t_f=25.350 ND, 14.1%
% off the Cartesian cross-formulation anchor (22.2248 ND). Stage B (cold
% nRev=3 seed, same node density) converges independently to 4.503 revs,
% t_f=22.2206 ND -- a 0.019% match. Min-time is genuinely multi-basin here
% (min-time revolution count is itself a discrete-ish degree of freedom
% the continuous NLP does not search globally): a Stage-A "success" per
% isGood (defect/termErr gate) is NOT sufficient evidence of having found
% the globally-relevant anchor, only of having found *a* locally-optimal
% one. FIX: Stage B is now ALWAYS attempted whenever Stage A is (not merely
% on Stage-A failure), and BOTH candidates are run through the same
% keep-if-improved selector (lower t_f wins, both must be certified) used
% for cross-invocation caching -- repurposed here as a same-run best-of-
% basins selector. This costs one extra ~1-min solve at 10 N (Stage B here
% is a small N=75 problem) and is skipped automatically at every other
% thrust in the eventual ladder anyway, since Stage A requires a cached
% fuel anchor for the EXACT thrustN being anchored, and the ladder's own
% per-rung order (anchor before fuel solve) means no such file exists
% except at 10 N, where Task 4 happened to produce one ahead of time. A
% cfg.alwaysTryStageB=false override is provided for a future caller
% confident enough in Stage A's basin to skip the cross-check (NOT used by
% default; not exercised by this task).
%
% FEASIBILITY-SELECTED BARRIER POLICY (review finding, binding -- supersedes
% the Cartesian file's round-number-indexed warmTight rule): every
% continuation round chooses warmTight from the MEASURED defect of the
% INCOMING warm-start point, not from which round index it is --
%   warmTight = incoming.success && incoming.maxDefect < 1e-6
% -- so a round warm-started from a point that did not exit cleanly
% (incoming.success==false: restoration-phase abort, Infeasible_Problem_
% Detected, Maximum_Iterations_Exceeded, ...) ALWAYS gets loose/adaptive-mu
% treatment, regardless of what its self-reported maxDefect happens to read.
% This is how "never reuse restoration-phase multipliers" is enforced here:
% casadi_lt_mee.m builds a fresh casadi.Opti() every call (no dual state
% persists across calls at all -- only the primal X/U/dL warm start does),
% so the only lever available to "not trust a restoration-phase exit's
% duals" is to keep IPOPT on adaptive-mu (loose) for the NEXT round rather
% than asking it to warm-start tight from an untrustworthy point.
%
% PER-ROUND RETAIN-IF-IMPROVED / RETRY-LOOSE (review finding, binding --
% DESIGN_thrust_ladder.md Phase 0 item 1): a continuation round's result is
% NO LONGER accepted unconditionally. After each round, the new iterate is
% retained (resetting the retry budget) only if it cleared the SAME
% decadeMin floor the stall guard has always used -- decadeImprove =
% log10(prevDefect) - log10(newDefect) >= decadeMin -- NOT merely
% newDefect<prevDefect (re-review finding, 2026-07-17: a bare-inequality
% "improved" let a sub-floor crawl, positive but below decadeMin, run to
% roundsMax without ever reaching the stall guard, silently defeating it).
% Anything below the floor -- sub-floor crawl OR true regression -- is
% retried EXACTLY ONCE, forced to the loose (adaptive-mu) regime regardless
% of what warmTight would otherwise have been -- mintime_guard_constants.m's
% maxLooseRetries (currently 1); the retained warm-start point for that
% retry is the BETTER of the two by defect (a sub-floor-but-positive
% candidate is kept as the new baseline even though it didn't reset the
% retry budget). Only if that loose retry ALSO fails to clear the floor
% does the existing stall guard (below) get to fire -- which it always
% does at that point, since not-clearing-the-floor is exactly the stall
% guard's own trigger condition. The retain/retry decision itself is a
% small pure function, round_advance_decision.m, factored out specifically
% so it is unit-testable with synthetic defect numbers without a solve
% (test_mintime_mee_guard.m) -- see that file for the improved /
% regressed-first-time / sub-floor-crawl / regressed-after-retry cases.
%
% Stall/converge guard: identical arithmetic to run_mintime.m's continuation
% loop, reusing mintime_guard_constants.m directly (NOT a hardcoded local
% copy -- pinned by test_mintime_mee_guard.m).
%
% CACHING: converged-only (never cache an uncertified result, campaign
% rule), config-fingerprinted (mirrors run_transfer_mee.m's check_cache_fp),
% and keep-if-improved: the final anchor file is only overwritten by a new
% solve that is BOTH certified and has a lower (better, for a minimization)
% t_f,min than whatever is already cached under the same tag. Per-round
% solver outputs are ALSO config-fingerprinted (same WARN-only backward-
% compat pattern for pre-fix round files with no stored fingerprint) and
% cached (resDir/<tag>_A_round%02d.mat / <tag>_B_round%02d.mat) so a killed
% process (sporadic MEX init crash) can be resumed by simply re-invoking
% this function -- already-solved rounds are loaded, not re-solved, and the
% retain/retry state machine is reconstructed by replaying the cached
% rounds in order (round_advance_decision is pure, so replay reproduces the
% same keep/retry decisions deterministically from the same defect trace).
%
% INPUTS:  thrustN     - max thrust [N]
%          nodesPerRev - node density [scalar, default 25 -- Task 5's
%                        production density]
%          cfg         - optional struct: .m0kg (1500), .ispS (2000),
%                        .maxIter (3000), .printLevel (3), .fuelTag (auto,
%                        see mee_fuel_tag.m), .nRevSeed (3, Stage B cold
%                        seed's target revs), .seedThrB (0.4, Stage B cold
%                        seed's constant throttle -- see STAGE-B SEED
%                        THROTTLE FIX above; NOT the mintime physics, only
%                        the raw-seed conditioning), .warmStartAnchor
%                        (optional struct .X [7x(Np+1)] .U [4x(Np+1)] .dL
%                        [scalar, ALREADY C-law rescaled by the caller]
%                        .N [scalar, Np] -- the PREVIOUS ladder rung's own
%                        converged min-time anchor result; see ANCHOR
%                        WARM-START FIX below), .tag (auto, see below),
%                        .alwaysTryStageB (true -- see BASIN-MULTIPLICITY
%                        FINDING above; set false to trust a converged
%                        Stage A outright and skip the Stage-B cross-check)
% OUTPUTS: out - struct: .tfmin (ND) .tfmin_h .dL_mt .revs .thrustN
%          .nodesPerRev .N .stage ('A'|'B') .continuationRounds .certified
%          .solverOut (full casadi_lt_mee output, for downstream
%          warm-starting e.g. run_ladder.m) .fp
%
% REFERENCES: [1] earth_elliptic_to_geo/run_mintime.m (Cartesian template;
%   two-stage/continuation recipe, stall guard, seedAnchor precedent).
%   [2] earth_elliptic_to_geo/mintime_guard_constants.m (shared guard pair).
%   [3] earth_elliptic_to_geo/casadi_lt_mee.m (MEE solver core, Task 3).
%   [4] DESIGN_thrust_ladder.md sec "Phase 2" (feasibility-selected barrier
%       policy, Phase-0 review finding, carried into the MEE driver).
%   [5] .superpowers/sdd/task-6-brief.md (this task's spec).
%   [6] earth_elliptic_to_geo/mee_fuel_tag.m (shared fuel-tag convention,
%       also used by run_ladder.m).
%   [7] earth_elliptic_to_geo/round_advance_decision.m (retain-if-improved /
%       retry-loose decision, Phase 0 item 1, factored out for unit testing).
if nargin < 2 || isempty(nodesPerRev), nodesPerRev = 25; end
if nargin < 3, cfg = struct(); end
d = @(f,v) getdef7(cfg, f, v);

m0kg     = d('m0kg', 1500);
ispS     = d('ispS', 2000);
maxIter  = d('maxIter', 3000);
printLvl = d('printLevel', 3);
fuelTag  = d('fuelTag', mee_fuel_tag(thrustN));
nRevSeed = d('nRevSeed', 3);
seedThrB = d('seedThrB', 0.4);   % see STAGE-B SEED THROTTLE FIX below
warmStartAnchor = d('warmStartAnchor', []);   % see ANCHOR WARM-START FIX below
tag      = d('tag', sprintf('MEE_mintime_T%d', round(10*thrustN)));
alwaysTryB = d('alwaysTryStageB', true);   % see BASIN-MULTIPLICITY FINDING above

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[roundsMax, decadeMin, maxLooseRetries] = mintime_guard_constants();   % shared, NOT duplicated

fp = struct('thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'm0kg', m0kg, ...
    'ispS', ispS, 'maxIter', maxIter, 'fuelTag', fuelTag, 'nRevSeed', nRevSeed, ...
    'seedThrB', seedThrB, 'roundsMax', roundsMax, 'decadeMin', decadeMin, ...
    'alwaysTryStageB', alwaysTryB);
if ~isempty(warmStartAnchor)
    fp.warmStartAnchorDL = warmStartAnchor.dL;   % fingerprint the warm-start input itself
    fp.warmStartAnchorN  = warmStartAnchor.N;
end

finalFile = fullfile(resDir, [tag '.mat']);
if isfile(finalFile)
    S = load(finalFile);
    check_cache_fp_mt(S, fp, finalFile, tag);
    out = S.out;
    fprintf(['cached %s: stage=%s tfmin=%.4f ND (%.1f h) dL=%.4f revs=%.3f ' ...
             '(%d continuation round(s))\n'], finalFile, out.stage, out.tfmin, ...
             out.tfmin_h, out.dL_mt, out.revs, out.continuationRounds);
    return;
end

par = kepler_lt_params(thrustN, m0kg, ispS);
isGood = @is_certified;   % shared with pack_candidate's .certified (Fix 3, no divergent copies)

candidates = {};

% --- Stage A: warm-start from the certified fuel solution at this thrust ---
fuelFile = fullfile(resDir, [fuelTag '.mat']);
stageAOK = false;
if isfile(fuelFile)
    Sf = load(fuelFile);
    fres = Sf.res;
    sigmaA = fres.sigma;
    X0A = fres.fuel.X;  U0A = fres.fuel.U;  U0A(4,:) = 1;  dL0A = fres.fuel.dL;
    x0state = X0A(:,1);
    fprintf('MINTIME_MEE Stage A: T=%g N, warm-started from fuel anchor %s (N=%d, dL0=%.4f)\n', ...
            thrustN, fuelTag, size(X0A,2)-1, dL0A);
    try
        [outA, roundA] = mintime_mee_continue(sigmaA, X0A, U0A, dL0A, x0state, par, ...
            isGood, roundsMax, decadeMin, maxLooseRetries, maxIter, printLvl, resDir, tag, 'A', fp);
        candidates{end+1} = pack_candidate(outA, thrustN, nodesPerRev, size(X0A,2)-1, 'A', roundA, par); %#ok<AGROW>
        stageAOK = true;
    catch ME_A
        fprintf('MINTIME_MEE Stage A FAILED (%s)\n', ME_A.message);
    end
else
    fprintf('MINTIME_MEE Stage A SKIPPED: no cached fuel anchor %s\n', fuelFile);
end

% --- Stage B: EITHER a fresh cold tangential seed OR (Task 7 addition, see
% ANCHOR WARM-START FIX below) a mesh-refined/rescaled warm start from the
% PREVIOUS ladder rung's own converged min-time anchor. Per the BASIN-
% MULTIPLICITY FINDING above, this is run whenever Stage A is unavailable/
% failed (the brief's literal fallback case) AND, by default
% (cfg.alwaysTryStageB), even when Stage A already converged -- a converged
% Stage A result is only a LOCAL optimum candidate, not proof of having
% found the globally-relevant one, so a second independent starting shape
% is cheap insurance.
%
% ANCHOR WARM-START FIX (2026-07-17, Task 7's first live ladder descent):
% the raw cold tangential seed -- even with the STAGE-B SEED THROTTLE FIX
% above (thr=0.4, well-conditioned) -- turned out to be a GENUINELY hard
% NLP start at 5 N, independent of seed conditioning: offline probing
% showed dual infeasibility exploding past 1e6 within 30 iterations and the
% defect stuck at ~1.3e-2 (essentially zero net progress) across multiple
% full-maxIter continuation rounds -- a real stall, not a fast-converging
% problem merely needing more compute. Root cause: the raw tangential-
% steering seed, however well-conditioned numerically, is nothing like the
% true (bang-bang, costate-driven) min-time control law's shape, and its
% large terminal-manifold miss (P ending well short of 1) gives Newton a
% hard multi-constraint target to hit from a poor starting direction. FIX:
% when cfg.warmStartAnchor is supplied (struct .X/.U/.dL/.N -- the PREVIOUS
% rung's own CONVERGED min-time anchor result, already thr=1 throughout,
% already satisfying the terminal manifold EXACTLY), Stage B mesh-refines
% THAT trajectory (same interp1 pattern as run_transfer_mee.m's fuel
% cfg.warmStart: linear for X, linear for the RTN thrust-direction rows of
% U -- thr itself is irrelevant, forced to 1 regardless) onto a grid sized
% by warmStartAnchor.dL/(2*pi) (already C-law rescaled by the caller,
% run_ladder.m: dL_guess = dL_mt(T_prev)*(T_prev/T_new)) at cfg.nRevSeed
% node density, INSTEAD OF calling mee_seed. This is a self-similar-shape,
% already-optimal-adjacent starting point (same rationale as the fuel warm
% start) rather than an arbitrary constant-throttle spiral, and is
% expected to be dramatically better conditioned. When warmStartAnchor is
% NOT supplied (e.g. the very first/10 N rung, which has no predecessor),
% Stage B falls back to the original raw cold tangential seed unchanged.
if ~stageAOK || alwaysTryB
    if ~isempty(warmStartAnchor)
        revsGuessB = warmStartAnchor.dL / (2*pi);
        N = round(nodesPerRev * revsGuessB);
        sigmaPrevB = linspace(0, 1, warmStartAnchor.N + 1).';
        sigmaB = linspace(0, 1, N + 1).';
        X0B    = interp1(sigmaPrevB, warmStartAnchor.X.',       sigmaB, 'linear').';
        Ubeta  = interp1(sigmaPrevB, warmStartAnchor.U(1:3,:).', sigmaB, 'linear').';
        Uthr   = interp1(sigmaPrevB, warmStartAnchor.U(4,:).',   sigmaB, 'nearest').';
        U0B    = [Ubeta; Uthr];
        dL0B   = warmStartAnchor.dL;
        x0state = X0B(:,1);
        fprintf(['MINTIME_MEE Stage B: T=%g N, WARM-STARTED from prior rung''s converged ' ...
                 'anchor (dL_guess=%.4f rad -> revs_guess=%.4f, N=%d nodes, %d nodes/rev)\n'], ...
                thrustN, dL0B, revsGuessB, N, nodesPerRev);
    else
        N = round(nodesPerRev * nRevSeed);
        [sigmaB, X0B, U0B, dL0B, seedInfoB] = mee_seed(par, struct('thr', seedThrB, ...
            'betaMode', 'tangential', 'nRev', nRevSeed, 'N', N));
        x0state = X0B(:,1);
        fprintf(['MINTIME_MEE Stage B: T=%g N, cold tangential thr=%.2g seed, nRev=%g target ' ...
                 '(N=%d, dL0=%.4f, achieved seed nRev=%.4f)\n'], ...
                thrustN, seedThrB, nRevSeed, N, dL0B, seedInfoB.nRev);
    end
    try
        [outB, roundB] = mintime_mee_continue(sigmaB, X0B, U0B, dL0B, x0state, par, ...
            isGood, roundsMax, decadeMin, maxLooseRetries, maxIter, printLvl, resDir, tag, 'B', fp);
        candidates{end+1} = pack_candidate(outB, thrustN, nodesPerRev, N, 'B', roundB, par); %#ok<AGROW>
    catch ME_B
        fprintf('MINTIME_MEE Stage B FAILED (%s)\n', ME_B.message);
    end
end

assert(~isempty(candidates), 'run_mintime_mee:allStagesFailed', ...
    'T=%g N: both Stage A and Stage B failed to converge -- no certified candidate (BLOCKED)', thrustN);

% keep-if-improved best-of-basins selection: lower t_f wins among certified
% candidates (repurposing the same selector used for cross-invocation
% caching -- see save_anchor_keep_if_improved).
tfs = cellfun(@(c) c.tfmin, candidates);
[~, ibest] = min(tfs);
for kc = 1:numel(candidates)
    tag_kc = ifelse(kc == ibest, 'BEST', 'discarded');
    fprintf('  candidate stage=%s tf=%.4f ND revs=%.3f -> %s\n', candidates{kc}.stage, ...
            candidates{kc}.tfmin, candidates{kc}.revs, tag_kc);
end
candidate = candidates{ibest};

candidate.fp = fp;
out = save_anchor_keep_if_improved(finalFile, candidate, tag);
fprintf(['MINTIME_MEE T=%g N: stage=%s tf=%.4f ND = %.1f h, dL=%.4f rad ' ...
         '(%.3f revs), defect=%.2e, %d continuation round(s), %s\n'], ...
        thrustN, out.stage, out.tfmin, out.tfmin_h, out.dL_mt, out.revs, ...
        out.solverOut.maxDefect, out.continuationRounds, out.solverOut.ipoptStatus);
end

% ---------------------------------------------------------------------------
function s = ifelse(cond, a, b)
% IFELSE  Tiny inline conditional for a printf label (MATLAB has no ternary).
if cond, s = a; else, s = b; end
end

% ---------------------------------------------------------------------------
function [out, round_] = mintime_mee_continue(sigma, X0, U0, dL0, x0state, par, ...
    isGood, roundsMax, decadeMin, maxLooseRetries, maxIter, printLvl, resDir, tag, label, fp)
% MINTIME_MEE_CONTINUE  Shared initial-solve + warm-continuation loop with
% the feasibility-selected barrier policy, the retain-if-improved / retry-
% loose mechanism (DESIGN_thrust_ladder.md Phase 0 item 1), and the shared
% stall/converge guard. Per-round outputs are config-fingerprinted and
% cached to resDir/<tag>_<label>_round%02d.mat (round 0 = the initial loose
% solve) so a killed process resumes without re-solving completed rounds.
%
% RETAIN-IF-IMPROVED / RETRY-LOOSE: each round's candidate (outNew) is
% retained (out = outNew, retry budget reset) only if round_advance_decision.m
% says decadeImprove = log10(out.maxDefect)-log10(outNew.maxDefect) cleared
% decadeMin -- the SAME floor the stall guard itself uses, not a bare
% newDefect<prevDefect test (a bare inequality would let a sub-floor crawl
% run to roundsMax without ever reaching the stall guard). On a sub-floor
% round (positive-but-below-floor progress, OR true regression), the
% warm-start point for the retry is the BETTER of {out, outNew} by defect,
% and the NEXT round is a forced-loose (warmTight=false) retry from that
% point, consuming one of mintime_guard_constants's maxLooseRetries. Only
% once that budget is exhausted (and the retry still didn't clear the
% floor) does the pre-existing decadeImprove/decadeMin stall guard get
% evaluated -- which it always fires at that point (not-clearing-the-floor
% is exactly decadeImprove<decadeMin), so exhausting the retry budget on a
% non-improving round is what actually throws run_mintime_mee:stall.
%
% INPUTS:  sigma [(N+1)x1]; X0/U0 [7/4 x (N+1)] warm-start state/control;
%          dL0 [scalar]; x0state [7x1] initial MEE state (opts.x0); par -
%          kepler_lt_params struct; isGood [function handle]; roundsMax/
%          decadeMin/maxLooseRetries [scalar, from mintime_guard_constants];
%          maxIter/printLvl [scalar]; resDir [char]; tag/label [char];
%          fp [struct, the caller's config fingerprint -- extended here
%          with .label/.N for the per-round cache fingerprint check]
% OUTPUTS: out - converged (or best-effort) casadi_lt_mee output struct;
%          round_ - number of continuation rounds beyond round 0 (includes
%          any loose-retry rounds -- they consume the same round budget)
%
% REFERENCES: [1] run_mintime.m>mintime_continue_only (pattern this mirrors,
%   adapted for casadi_lt_mee's opts contract and the feasibility-selected
%   barrier policy in place of the round-number-indexed warmTight rule).
%   [2] round_advance_decision.m (the retain/retry decision, unit-tested
%   separately). [3] DESIGN_thrust_ladder.md Phase 0 item 1 (the mandate).
fpRound = fp;  fpRound.label = label;  fpRound.N = size(X0,2) - 1;
fp = fpRound;   % shadow the input fp with the round-level fingerprint for the
                % remainder of this function -- the caller's copy is untouched
                % (MATLAB passes by value), and every save() below stores this
                % as field name 'fp', matching check_cache_fp_round's read and
                % the run_transfer_mee.m/check_cache_fp convention.

round0File = fullfile(resDir, sprintf('%s_%s_round%02d.mat', tag, label, 0));
if isfile(round0File)
    S = load(round0File);  out = S.out;
    check_cache_fp_round(S, fpRound, round0File, tag, label, 0);
    fprintf('  [cached] %s round 0\n', label);
else
    out = casadi_lt_mee(sigma, X0, U0, dL0, struct('par', par, 'mode', 'mintime', ...
        'x0', x0state, 'maxIter', maxIter, 'warmTight', false, 'printLevel', printLvl));
    save(round0File, 'out', 'fp');
end
fprintf('  %s round 0 (loose, adaptive mu): status=%s defect=%.3e termErr=%.3e\n', ...
        label, out.ipoptStatus, out.maxDefect, out.termErr);

round_ = 0;
defect0 = out.maxDefect;
gateDefect = 1e-8;
guardC = struct('decadeMin', decadeMin, 'maxLooseRetries', maxLooseRetries);
retriedThisStall = false;   % true iff the CURRENTLY retained `out` has
                             % already consumed its one loose-regime retry
                             % (Phase 0 item 1) without improving on it.
while ~isGood(out) && round_ < roundsMax
    round_ = round_ + 1;
    roundFile = fullfile(resDir, sprintf('%s_%s_round%02d.mat', tag, label, round_));
    if isfile(roundFile)
        S = load(roundFile);  outNew = S.outNew;
        check_cache_fp_round(S, fpRound, roundFile, tag, label, round_);
        fprintf('  [cached] %s round %d\n', label, round_);
    else
        % FEASIBILITY-SELECTED BARRIER POLICY, with the retry-loose override
        % (Phase 0 item 1): if this round IS the one-shot loose retry for a
        % regressed prior round, warmTight is forced false regardless of the
        % measured feasibility of the (retained) incoming point.
        if retriedThisStall
            warmTight = false;
        else
            warmTight = out.success && out.maxDefect < 1e-6;
        end
        outNew = casadi_lt_mee(sigma, out.X, out.U, out.dL, struct('par', par, ...
            'mode', 'mintime', 'x0', x0state, 'maxIter', maxIter, ...
            'warmTight', warmTight, 'printLevel', printLvl));
        save(roundFile, 'outNew', 'fp');
    end
    warmTightUsed = ~retriedThisStall && out.success && out.maxDefect < 1e-6;
    fprintf(['  %s round %d (warmTight=%d, looseRetry=%d): defect %.3e -> %.3e, ' ...
             'termErr=%.3e, status=%s\n'], label, round_, warmTightUsed, retriedThisStall, ...
            out.maxDefect, outNew.maxDefect, outNew.termErr, outNew.ipoptStatus);
    decadeImprove = log10(max(out.maxDefect, realmin)) - log10(max(outNew.maxDefect, realmin));
    cumDecades = log10(max(defect0, realmin)) - log10(max(outNew.maxDefect, realmin));
    decadesRemaining = log10(max(outNew.maxDefect, realmin)) - log10(gateDefect);
    fprintf(['    cumulative: %.2f decades closed since round 0 (defect0=%.3e), ' ...
             '~%.2f decades remaining to gate (%d/%d rounds used)\n'], ...
            cumDecades, defect0, max(decadesRemaining, 0), round_, roundsMax);

    [keepNew, retryLoose] = round_advance_decision(out.maxDefect, outNew.maxDefect, ...
        retriedThisStall, guardC);
    if keepNew
        out = outNew;
        retriedThisStall = false;
    elseif retryLoose
        % Sub-floor crawl OR true regression -- either way this round does
        % NOT clear decadeMin, so it does not reset the retry budget. But
        % the RETAINED point for the retry is the better of the two by
        % defect (review finding, re-review): a sub-floor-but-positive
        % candidate is still numerically better than the incoming point and
        % should not be thrown away, it just doesn't count as "progress"
        % for stall-guard purposes.
        if outNew.maxDefect < out.maxDefect
            fprintf(['    round %d sub-floor progress (defect %.3e -> %.3e, %.2f decades ' ...
                     '< floor %.2f) -- retaining the improved point as the new baseline, ' ...
                     'retrying once with the loose regime before declaring a stall\n'], ...
                    round_, out.maxDefect, outNew.maxDefect, decadeImprove, decadeMin);
            out = outNew;
        else
            fprintf(['    round %d regressed (defect %.3e -> %.3e) -- keeping the retained ' ...
                     'prior iterate, retrying once with the loose regime before declaring a stall\n'], ...
                    round_, out.maxDefect, outNew.maxDefect);
        end
        retriedThisStall = true;
    else
        if ~isGood(outNew) && decadeImprove < decadeMin
            error('run_mintime_mee:stall', ['continuation stalled at %s round %d ' ...
                '(loose retry exhausted, maxLooseRetries=%d): defect %.3e -> %.3e ' ...
                '(%.2f decades, need >=%.2f), termErr=%.3e, status=%s'], label, round_, ...
                maxLooseRetries, out.maxDefect, outNew.maxDefect, decadeImprove, decadeMin, ...
                outNew.termErr, outNew.ipoptStatus);
        end
        out = outNew;
        retriedThisStall = false;
    end
end
if ~isGood(out)
    error('run_mintime_mee:noConverge', ...
        ['%s failed to converge after %d continuation round(s): defect=%.3e ' ...
         'termErr=%.3e status=%s'], label, round_, out.maxDefect, out.termErr, out.ipoptStatus);
end
end

% ---------------------------------------------------------------------------
function tf = is_certified(o)
% IS_CERTIFIED  Single shared certification predicate (review finding, Fix
% 3): the feasibility + terminal-accuracy gate used BOTH as the
% continuation loop's isGood test and pack_candidate's .certified flag, so
% the two definitions cannot silently diverge.
%
% INPUTS:  o - casadi_lt_mee output struct (.success/.maxDefect/.termErr)
% OUTPUTS: tf - true iff o is a certified (converged, tight) result [logical]
tf = o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;
end

% ---------------------------------------------------------------------------
function cand = pack_candidate(o, thrustN, nodesPerRev, N, stage, roundsUsed, par)
% PACK_CANDIDATE  Build the anchor summary struct from a converged
% casadi_lt_mee output (mirrors run_mintime.m's revs/dL_mt reporting).
cand = struct('tfmin', o.tf, 'tfmin_h', o.tf*par.TU_s/3600, 'dL_mt', o.dL, ...
    'revs', o.dL/(2*pi), 'thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'N', N, ...
    'stage', stage, 'continuationRounds', roundsUsed, 'certified', is_certified(o), ...
    'solverOut', o);
end

% ---------------------------------------------------------------------------
function out = save_anchor_keep_if_improved(finalFile, candidate, tag)
% SAVE_ANCHOR_KEEP_IF_IMPROVED  Converged-only, keep-if-improved cache write:
% a candidate that is not certified is never written (campaign rule -- an
% uncertified anchor must never sit in a path a downstream task could
% silently load). If a cache already exists under this tag (a race or a
% resumed run with a stale finalFile), overwrite it only when the new
% candidate is BOTH certified and has a lower (better, for a min-time
% objective) t_f,min than the cached value; otherwise keep the existing
% file and return it.
%
% INPUTS:  finalFile - path [char]; candidate - struct from pack_candidate,
%          plus .fp; tag - [char, for messages]
% OUTPUTS: out - the struct actually left on disk (candidate or the
%          pre-existing cache, whichever is authoritative)
if ~candidate.certified
    error('run_mintime_mee:uncertified', ['%s: candidate NOT certified ' ...
        '(defect=%.2e termErr=%.2e status=%s) -- campaign rule: never cache ' ...
        'an uncertified anchor'], tag, candidate.solverOut.maxDefect, ...
        candidate.solverOut.termErr, candidate.solverOut.ipoptStatus);
end
if isfile(finalFile)
    S = load(finalFile);
    if S.out.certified && S.out.tfmin <= candidate.tfmin
        fprintf(['%s: existing cached anchor (tfmin=%.4f) is already as good ' ...
                 'or better than the new candidate (tfmin=%.4f) -- keeping it\n'], ...
                tag, S.out.tfmin, candidate.tfmin);
        out = S.out;
        return;
    end
    fprintf('%s: new candidate (tfmin=%.4f) improves on cached (tfmin=%.4f) -- overwriting\n', ...
            tag, candidate.tfmin, S.out.tfmin);
end
out = candidate;
save(finalFile, 'out');
end

% ---------------------------------------------------------------------------
function check_cache_fp_mt(S, fp, file, tag)
% CHECK_CACHE_FP_MT  Fail-loud cache-fingerprint guard (mirrors
% run_transfer_mee.m's check_cache_fp). BACKWARD COMPAT: a pre-fix cache
% with no .fp field only WARNs and is trusted as-is.
if ~isfield(S, 'out') || ~isfield(S.out, 'fp')
    warning('run_mintime_mee:noCachedFingerprint', ['%s has no cached config ' ...
        'fingerprint -- trusting it because tag=''%s'' matches; use a new ' ...
        'tag to regain fingerprint protection'], file, tag);
    return;
end
sfp = S.out.fp;
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(sfp, f) || ~isequal(sfp.(f), fp.(f))
        error('run_mintime_mee:fingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s: field ''%s'' differs between the ' ...
            'cache and the current config -- stale cache under tag=''%s''; ' ...
            'delete the file or use a new tag'], file, f, tag);
    end
end
end

% ---------------------------------------------------------------------------
function check_cache_fp_round(S, fpRound, file, tag, label, roundIdx)
% CHECK_CACHE_FP_ROUND  Fail-loud cache-fingerprint guard for the PER-ROUND
% cache files (<tag>_<label>_round%02d.mat) -- review finding (Fix 2): these
% were previously loaded by round number with no config check at all,
% unlike the final anchor file (check_cache_fp_mt above) and
% run_transfer_mee.m's check_cache_fp, which this mirrors exactly.
% BACKWARD COMPAT: a pre-fix round file with no stored .fp field (e.g. the
% pre-existing MEE_mintime_T100_A/B_round00.mat) only WARNs and is trusted
% as-is -- tag/label/round index are already an exact match (that's how the
% file was found by name), and no per-field comparison is possible without
% a stored fingerprint.
%
% INPUTS:  S - loaded round-cache struct; fpRound - current round-level
%          config fingerprint [struct]; file - path [char]; tag/label
%          [char]; roundIdx - round number [scalar, for messages]
if ~isfield(S, 'fp')
    warning('run_mintime_mee:noCachedRoundFingerprint', ['%s (round %d, %s) has ' ...
        'no cached config fingerprint (pre-fix round cache) -- trusting it ' ...
        'because tag/label/round match; use a new tag to regain fingerprint ' ...
        'protection'], file, roundIdx, label);
    return;
end
flds = fieldnames(fpRound);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f) || ~isequal(S.fp.(f), fpRound.(f))
        error('run_mintime_mee:roundFingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s (round %d, %s): field ''%s'' differs ' ...
            'between the cache and the current config -- stale round cache ' ...
            'under tag=''%s''; delete the file or use a new tag'], file, roundIdx, ...
            label, f, tag);
    end
end
end

% ---------------------------------------------------------------------------
function v = getdef7(s, f, dflt)
% GETDEF7  Optional-field default (mirrors casadi_lt_2body's helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
