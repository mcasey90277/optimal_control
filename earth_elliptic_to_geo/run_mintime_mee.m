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
%   Stage B - a FRESH cold tangential seed (mee_seed, thr=1,
%             betaMode='tangential', nRev=3 -- an arbitrary but reasonable
%             all-burn starting shape; DeltaL is free, so the solve is not
%             committed to 3 revs, only seeded near it) at the same node
%             density, then the same continuation recipe.
%
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
% Stall/converge guard: identical arithmetic to run_mintime.m's continuation
% loop, reusing mintime_guard_constants.m directly (NOT a hardcoded local
% copy -- pinned by test_mintime_mee_guard.m).
%
% CACHING: converged-only (never cache an uncertified result, campaign
% rule), config-fingerprinted (mirrors run_transfer_mee.m's check_cache_fp),
% and keep-if-improved: the final anchor file is only overwritten by a new
% solve that is BOTH certified and has a lower (better, for a minimization)
% t_f,min than whatever is already cached under the same tag. Per-round
% solver outputs are also cached (resDir/<tag>_A_round%02d.mat /
% <tag>_B_round%02d.mat) so a killed process (sporadic MEX init crash) can
% be resumed by simply re-invoking this function -- already-solved rounds
% are loaded, not re-solved.
%
% INPUTS:  thrustN     - max thrust [N]
%          nodesPerRev - node density [scalar, default 25 -- Task 5's
%                        production density]
%          cfg         - optional struct: .m0kg (1500), .ispS (2000),
%                        .maxIter (3000), .printLevel (3), .fuelTag (auto,
%                        see mee_fuel_tag.m), .nRevSeed (3, Stage B cold
%                        seed's target revs), .tag (auto, see below),
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
if nargin < 2 || isempty(nodesPerRev), nodesPerRev = 25; end
if nargin < 3, cfg = struct(); end
d = @(f,v) getdef7(cfg, f, v);

m0kg     = d('m0kg', 1500);
ispS     = d('ispS', 2000);
maxIter  = d('maxIter', 3000);
printLvl = d('printLevel', 3);
fuelTag  = d('fuelTag', mee_fuel_tag(thrustN));
nRevSeed = d('nRevSeed', 3);
tag      = d('tag', sprintf('MEE_mintime_T%d', round(10*thrustN)));
alwaysTryB = d('alwaysTryStageB', true);   % see BASIN-MULTIPLICITY FINDING above

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[roundsMax, decadeMin] = mintime_guard_constants();   % shared, NOT duplicated

fp = struct('thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'm0kg', m0kg, ...
    'ispS', ispS, 'maxIter', maxIter, 'fuelTag', fuelTag, 'nRevSeed', nRevSeed, ...
    'roundsMax', roundsMax, 'decadeMin', decadeMin, 'alwaysTryStageB', alwaysTryB);

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
isGood = @(o) o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;

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
            isGood, roundsMax, decadeMin, maxIter, printLvl, resDir, tag, 'A');
        candidates{end+1} = pack_candidate(outA, thrustN, nodesPerRev, size(X0A,2)-1, 'A', roundA, par); %#ok<AGROW>
        stageAOK = true;
    catch ME_A
        fprintf('MINTIME_MEE Stage A FAILED (%s)\n', ME_A.message);
    end
else
    fprintf('MINTIME_MEE Stage A SKIPPED: no cached fuel anchor %s\n', fuelFile);
end

% --- Stage B: fresh cold tangential thr=1 seed. Per the BASIN-MULTIPLICITY
% FINDING above, this is run whenever Stage A is unavailable/failed (the
% brief's literal fallback case) AND, by default (cfg.alwaysTryStageB),
% even when Stage A already converged -- a converged Stage A result is only
% a LOCAL optimum candidate, not proof of having found the globally-
% relevant one, so a second independent starting shape is cheap insurance.
if ~stageAOK || alwaysTryB
    N = round(nodesPerRev * nRevSeed);
    [sigmaB, X0B, U0B, dL0B, seedInfoB] = mee_seed(par, struct('thr', 1, ...
        'betaMode', 'tangential', 'nRev', nRevSeed, 'N', N));
    x0state = X0B(:,1);
    fprintf(['MINTIME_MEE Stage B: T=%g N, cold tangential thr=1 seed, nRev=%g target ' ...
             '(N=%d, dL0=%.4f, achieved seed nRev=%.4f)\n'], ...
            thrustN, nRevSeed, N, dL0B, seedInfoB.nRev);
    try
        [outB, roundB] = mintime_mee_continue(sigmaB, X0B, U0B, dL0B, x0state, par, ...
            isGood, roundsMax, decadeMin, maxIter, printLvl, resDir, tag, 'B');
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
    isGood, roundsMax, decadeMin, maxIter, printLvl, resDir, tag, label)
% MINTIME_MEE_CONTINUE  Shared initial-solve + warm-continuation loop with
% the feasibility-selected barrier policy and the shared stall/converge
% guard. Per-round outputs are cached to resDir/<tag>_<label>_round%02d.mat
% (round 0 = the initial loose solve) so a killed process resumes without
% re-solving completed rounds.
%
% INPUTS:  sigma [(N+1)x1]; X0/U0 [7/4 x (N+1)] warm-start state/control;
%          dL0 [scalar]; x0state [7x1] initial MEE state (opts.x0); par -
%          kepler_lt_params struct; isGood [function handle]; roundsMax/
%          decadeMin [scalar, from mintime_guard_constants]; maxIter/
%          printLvl [scalar]; resDir [char]; tag/label [char]
% OUTPUTS: out - converged (or best-effort) casadi_lt_mee output struct;
%          round_ - number of continuation rounds beyond round 0
%
% REFERENCES: [1] run_mintime.m>mintime_continue_only (pattern this mirrors,
%   adapted for casadi_lt_mee's opts contract and the feasibility-selected
%   barrier policy in place of the round-number-indexed warmTight rule).
round0File = fullfile(resDir, sprintf('%s_%s_round%02d.mat', tag, label, 0));
if isfile(round0File)
    S = load(round0File);  out = S.out;
    fprintf('  [cached] %s round 0\n', label);
else
    out = casadi_lt_mee(sigma, X0, U0, dL0, struct('par', par, 'mode', 'mintime', ...
        'x0', x0state, 'maxIter', maxIter, 'warmTight', false, 'printLevel', printLvl));
    save(round0File, 'out');
end
fprintf('  %s round 0 (loose, adaptive mu): status=%s defect=%.3e termErr=%.3e\n', ...
        label, out.ipoptStatus, out.maxDefect, out.termErr);

round_ = 0;
defect0 = out.maxDefect;
gateDefect = 1e-8;
while ~isGood(out) && round_ < roundsMax
    round_ = round_ + 1;
    roundFile = fullfile(resDir, sprintf('%s_%s_round%02d.mat', tag, label, round_));
    if isfile(roundFile)
        S = load(roundFile);  outNew = S.out;
        fprintf('  [cached] %s round %d\n', label, round_);
    else
        % FEASIBILITY-SELECTED BARRIER POLICY: warmTight only if the
        % INCOMING point exited cleanly (out.success) AND its measured
        % defect is already tight (<1e-6) -- never from the round index,
        % never on a restoration/failed exit (out.success==false forces
        % loose/adaptive-mu regardless of the reported defect number).
        warmTight = out.success && out.maxDefect < 1e-6;
        outNew = casadi_lt_mee(sigma, out.X, out.U, out.dL, struct('par', par, ...
            'mode', 'mintime', 'x0', x0state, 'maxIter', maxIter, ...
            'warmTight', warmTight, 'printLevel', printLvl));
        save(roundFile, 'outNew');
    end
    fprintf(['  %s round %d (warmTight=%d): defect %.3e -> %.3e, termErr=%.3e, status=%s\n'], ...
            label, round_, out.success && out.maxDefect < 1e-6, out.maxDefect, ...
            outNew.maxDefect, outNew.termErr, outNew.ipoptStatus);
    decadeImprove = log10(max(out.maxDefect, realmin)) - log10(max(outNew.maxDefect, realmin));
    cumDecades = log10(max(defect0, realmin)) - log10(max(outNew.maxDefect, realmin));
    decadesRemaining = log10(max(outNew.maxDefect, realmin)) - log10(gateDefect);
    fprintf(['    cumulative: %.2f decades closed since round 0 (defect0=%.3e), ' ...
             '~%.2f decades remaining to gate (%d/%d rounds used)\n'], ...
            cumDecades, defect0, max(decadesRemaining, 0), round_, roundsMax);
    if ~isGood(outNew) && decadeImprove < decadeMin
        error('run_mintime_mee:stall', ['continuation stalled at %s round %d: ' ...
            'defect %.3e -> %.3e (%.2f decades, need >=%.2f), termErr=%.3e, status=%s'], ...
            label, round_, out.maxDefect, outNew.maxDefect, decadeImprove, decadeMin, ...
            outNew.termErr, outNew.ipoptStatus);
    end
    out = outNew;
end
if ~isGood(out)
    error('run_mintime_mee:noConverge', ...
        ['%s failed to converge after %d continuation round(s): defect=%.3e ' ...
         'termErr=%.3e status=%s'], label, round_, out.maxDefect, out.termErr, out.ipoptStatus);
end
end

% ---------------------------------------------------------------------------
function cand = pack_candidate(o, thrustN, nodesPerRev, N, stage, roundsUsed, par)
% PACK_CANDIDATE  Build the anchor summary struct from a converged
% casadi_lt_mee output (mirrors run_mintime.m's revs/dL_mt reporting).
cand = struct('tfmin', o.tf, 'tfmin_h', o.tf*par.TU_s/3600, 'dL_mt', o.dL, ...
    'revs', o.dL/(2*pi), 'thrustN', thrustN, 'nodesPerRev', nodesPerRev, 'N', N, ...
    'stage', stage, 'continuationRounds', roundsUsed, 'certified', ...
    o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8, 'solverOut', o);
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
function v = getdef7(s, f, dflt)
% GETDEF7  Optional-field default (mirrors casadi_lt_2body's helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
