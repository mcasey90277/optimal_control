function res = run_mintime(thrustN, hx0, N, seedAnchor, routeB)
% RUN_MINTIME  Free-L (manifold) min-time anchor at one thrust level.
%
% t_f,min sets every c_tf scale downstream (paper's TfMin is free-longitude).
%
% METHODOLOGICAL FINDING (Task 8, 2026-07-16): a direct one-shot 'mintime'
% solve of casadi_lt_2body against the free-longitude insertion MANIFOLD,
% warm-started straight from the cold tangential seed_2body seed, does NOT
% converge (tested coplanar case: 3000 iters, Maximum_Iterations_Exceeded,
% maxDefect ~2.6e-3, termErr ~0.14 -- gates missed by 5+ orders of magnitude).
% The manifold terminal set is too weakly constraining relative to the
% free-time mintime objective for the cold seed to find the right basin.
%
% What DOES converge (used here as the real algorithm, not a fallback):
%   Stage 1 - 'mintime' vs a FIXED terminal rendezvous at the seed's own
%             arrival longitude (geo_terminal('fixed', p, sinfo.Larr)) --
%             zero terminal gap by construction, so this is an easy warm-up
%             even if IPOPT does not fully converge it (max-iter here is
%             fine; only the trajectory SHAPE is used downstream).
%   Stage 2 - 'mintime' vs the free-longitude MANIFOLD, warm-started from
%             stage 1's (possibly non-converged) trajectory, warmTight=true.
%   Continuation - if stage 2 itself does not reach gate tolerance, repeat
%             the manifold solve warm-started from the previous iterate
%             (warmTight=true) up to CONTIN_MAX_ROUNDS times. Empirically
%             this needed 0 extra rounds for the 3D case (stage 2 converged
%             directly) and 1 extra round for the coplanar case (stage 2
%             stalled at maxDefect 2.9e-4, one more warm continuation reached
%             2.6e-12). A round that improves maxDefect by less than one
%             decade without reaching the gate tolerance is treated as a
%             stall and raises an error (never grind indefinitely).
%
% Caching is CONVERGED-ONLY: a result is written to the canonical cache path
% only if success && maxDefect<1e-8 && termErr<1e-8, so a non-converged
% iterate can never sit in a path that downstream tasks silently load.
%
% THRUST CONTINUATION (Task 14 controller triage, 2026-07-16, TWO rounds):
% the brute cold-seed stage1+stage2+continuation recipe above stalls outright
% at lower thrust (5 N, N=1200: continuation round 1 improved defect only
% 0.24 decades, far under the 1-decade/round stall guard -- 2.5 N was never
% reached). Round-1 fix attempt (seedAnchor warm-starting stage 2 DIRECTLY,
% skipping stage 1) was itself wrong and re-stalled (5 N: round 1 went
% Infeasible_Problem_Detected, -0.02 decades) -- the fixed-rendezvous stage-1
% warm-up turns out to be load-bearing even when the warm start is already a
% converged min-time SHAPE at a different thrust, because it resettles the
% dynamics defects the time-stretch introduces (only X(8,:) is rescaled;
% X(9,:)=cScale is left as-is, so raw stretched state is dynamically
% inconsistent until stage 1 relaxes it). Round-2 fix (current):
%   1. seedAnchor path REINSTATES stage 1: 'mintime' vs a FIXED terminal at
%      the STRETCHED seed's OWN arrival longitude (positions are unchanged by
%      the stretch, so Lend = atan2(X0(2,end),X0(1,end)) already sits exactly
%      on the anchor's converged endpoint -- ​a zero-gap warm-up, same spirit
%      as the cold-seed recipe), THEN stage 2 (manifold, warmTight), then
%      continuation rounds.
%   2. The stall guard is RECALIBRATED for anchor-seeded continuation (this
%      was tightened further in round 3/4, see below): a round is only
%      flagged a stall if it improves less than ~0.15 decades (observed rate
%      at 5 N was 0.24 decades/round -- slower convergence at lower thrust is
%      expected and not itself pathological). Still converged-only caching.
%   3. If the anchor-seeded path (mesh = anchor's own, decoupled from the
%      fuel-stage N) still fails after (1)+(2), ONE fallback attempt is made:
%      a FRESH cold tangential seed at N=1200 (nodes-per-rev parity -- 8.4
%      revs at this thrust needs ~1200 for the same mesh density as the 10 N
%      anchor's 600-node/~4.2-rev mesh), still with stage 1 + the same
%      relaxed round/decade guard. If that also fails, the error propagates
%      out of run_mintime uncaught (BLOCKED -- report with both stall
%      trajectories; do not grind further).
%   4. Anchors are meant to be chained thrust-by-thrust (10 N -> 5 N -> 2.5 N
%      in run_ctf_sweep) -- a lower thrust is only attempted once its
%      immediate neighbor's anchor has actually converged and cached.
%
% seedAnchor PATH IS DEPRECATED (Task 14 controller triage round 3,
% 2026-07-16): it has a TOPOLOGY FLAW. Time-stretching a converged anchor
% only rescales X(8,:) (physical time) -- it does NOT add revolutions, so a
% 4.5-rev 10 N shape was being forced toward the ~8.4-rev 5 N min-time
% topology (same class of error as a wrong-rev fuel seed). The round-2 fix
% above (reinstating stage 1) improved the SYMPTOM (stage 1 no longer
% rejects the warm start outright) but not the root cause. CURRENT RECIPE:
% every thrust gets its OWN cold tangential seed at nodes-per-rev parity (N
% scaled with the thrust's own min-time rev count: 600 @ 10 N/~4.2 revs,
% 1200 @ 5 N/~8.4 revs, 2400 @ 2.5 N/~17 revs) -- the seed_2body event-stop
% spiral winds the CORRECT rev count by construction, by definition. The
% code below is kept (not deleted) for reference/future use but is no longer
% called by run_ctf_sweep.
%
% ROUND-3 RECALIBRATION: the N=1200/5 N cold-seed continuation genuinely
% converges (round 1: 0.24 decades, termErr already 3.2e-4) but at ~5
% min/round and ~0.24 decades/round, closing 5.7 decades needs ~24 rounds
% (~2 h) -- affordable, and empirically convergence rate tends to accelerate
% near the basin. Continuation cap raised 6 -> 24 (same 0.15-decade stall
% floor, same converged-only caching); a cumulative decades-closed figure is
% now printed every round so a long run's progress is legible in the log.
%
% ROUND-4 BUG FIX: the round-3 recalibration (CONTIN_MAX_ROUNDS=24,
% CONTIN_DECADE_MIN=0.15) was plumbed into the seedAnchor branch (and its
% N=1200 fallback) but NOT into the plain cold-seed branch below -- which,
% since seedAnchor is deprecated, is the branch every non-cached thrust
% actually runs through. That branch was still calling
% mintime_stage12_continue with the old tight 3-round/1.0-decade pair, so a
% 5 N round that legitimately cleared 0.244 decades (>= 0.15) was rejected
% against a floor of 1.0 instead. The guard arithmetic itself was correct
% throughout (pinned in test_stall_guard.m); only the wrong constants were
% reaching it. Fixed: one shared CONTIN_MAX_ROUNDS/CONTIN_DECADE_MIN pair now
% feeds every call site below, including the plain cold-seed path.
%
% ROUND-5: ROUTE-B (energy-warm-start min-time). With the guard bug fixed,
% the plain cold-seed continuation at 5 N genuinely ran (round 1: 0.24
% decades, matches test_stall_guard.m) but round 2 hit defect 5.336e-03 ->
% 4.877e-03 (0.04 decades, status=Infeasible_Problem_Detected, termErr~4e-16)
% -- a LEGITIMATE stall, not a guard bug: IPOPT's restoration phase declares
% local infeasibility around a ~5e-3 defect floor that the tangential cold
% seed's basin cannot climb out of. This exact signature (defect floor ~5e-3
% + false local-infeasibility) is documented in the parent tulip campaign for
% the min-time cell, cured there by the ELFO ROUTE-B pattern (see
% NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m,
% docs/superpowers/specs/2026-07-15-elfo-mintime-route-b-design.md): warm-
% start min-time from a CONVERGED SMOOTH ENERGY solution (mode='fixedtf',
% eps=1) instead of from a raw thrust-propagation seed. The energy problem is
% convex-ish and (per every energy solve run so far, e.g. test_energy_stage.m)
% converges cleanly; its trajectory SHAPE -- already satisfying the manifold
% terminal (term='manifold' used for the energy stage too) -- is then a much
% better min-time warm start than the cold tangential spiral. Implementation:
% the routeB argument (see INPUTS) triggers a dedicated branch that (1) solves
% the energy stage at the caller-given (sbar, tDur) with N nodes, term=
% 'manifold', tfTarget=tDur, warmTight=false; (2) warm-starts 'mintime' from
% that solution directly (U(4,:) forced to 1 -- 'mintime' mode pins s==1 as a
% hard constraint regardless, so this is a warm-start convenience, not a
% relaxation), warmTight=false for the first pass; (3) hands off to the SAME
% mintime_continue_only rounds/guard as every other path. This is the FINAL
% solver strategy for this task -- if it also stalls under the (already
% validated-correct) guard, the anchor is BLOCKED and reported as such rather
% than attempting a further solver variant.
%
% INPUTS:  thrustN - max thrust [N];  hx0 - initial hx (0 coplanar | 0.0612);
%          N - mesh segments for the cold-seed/Route-B recipe (default 600;
%          pass the thrust's own nodes-per-rev-parity mesh, e.g. 1200 @ 5 N,
%          2400 @ 2.5 N; ignored when seedAnchor succeeds -- the anchor's own
%          mesh is used; overridden to 1200 internally if the anchor path
%          falls back);  seedAnchor - DEPRECATED (see above), optional path
%          to a neighboring-thrust's cached mintime_*.mat; when given, both
%          stages warm-start from its time-stretched trajectory;  routeB -
%          optional struct .sbar .tDur [ND] to trigger the energy-warm-start
%          recipe (round 5, see above) instead of the plain cold-seed path;
%          tDur is also used as the energy stage's tfTarget (pick ~2x the
%          thrust's own min-time estimate from the C-law, T*tf~const, so the
%          energy problem sits comfortably interior)
% OUTPUTS: res - .out .tfmin .tfmin_h .dL_mt .revs (also saved/cached in
%          results/); .continuationRounds records how many extra warm
%          continuation rounds beyond stage 1+2 (or beyond the energy-warm-
%          started mintime solve, on the Route-B path) were needed;
%          .seedInfo is the seed_2body info struct on the cold-seed path, a
%          provenance struct (.continuationFrom .anchorThrustN .stretch
%          .fallback [.fallbackN .anchorAttemptError]) on the (deprecated)
%          anchor-seeded path, or (.routeB .sbar .tDur .energyDefect) on the
%          Route-B path.
%
% REFERENCES: [1] DESIGN.md sec 4 step 1.  [2] PLAN.md Task 8.
%   [3] task-8-brief.md contingency comment (origin of the stage-1/stage-2
%       idea; promoted here from documented fallback to the default path
%       after the one-shot manifold solve was shown not to converge).
%   [4] Task 14 controller triage rounds 1-5 (thrust-continuation fix,
%       Route-B). [5] NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m (Route-B
%       origin, parent tulip campaign).
if nargin < 3 || isempty(N), N = 600; end
if nargin < 4, seedAnchor = ''; end
if nargin < 5, routeB = []; end
% ROUND-4 BUG FIX (Task 14 controller triage): round 3 introduced a patient
% ANCHOR_MAX_ROUNDS/ANCHOR_DECADE_MIN pair but only plumbed it into the
% (now-deprecated) seedAnchor branch; the plain cold-seed branch below --
% which is the ACTUAL recipe every non-cached thrust now uses, since
% seedAnchor is deprecated -- still called mintime_stage12_continue with the
% old tight CONTINGENCY_MAX_ROUNDS=3/decadeMin=1 pair (calibrated for the 10 N
% single-shot case). That is why the 5 N round-1 result (9.353e-03 ->
% 5.336e-03 = 0.244 decades) was rejected against a floor of 1.0 instead of
% the intended 0.15 -- the guard arithmetic itself was correct throughout
% (see test_stall_guard.m); only the WRONG constants were reaching it. FIX:
% one shared patient recipe (CONTIN_MAX_ROUNDS/CONTIN_DECADE_MIN) used by
% every call site below -- the fast 10 N case still exits the continuation
% loop immediately via isGood (0 extra rounds needed, per the M0-M2 logs),
% so raising its ceiling and lowering its floor costs nothing there.
[CONTIN_MAX_ROUNDS, CONTIN_DECADE_MIN] = mintime_guard_constants();
% Documented in mintime_guard_constants.m: ~0.24 decades/round at 5 N, N=1200,
% ~5 min/round; recalibrated from the original 1-decade/round 10 N floor.

resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
tag = sprintf('mintime_T%d_i%d', round(10*thrustN), round(hx0 > 0)*7);
fn  = fullfile(resDir, [tag '.mat']);
if isfile(fn), S = load(fn); res = S.res; fprintf('cached %s\n', fn); return; end

p  = kepler_lt_params(thrustN, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, hx0, 0, pi, p.mu);
rv0 = [r0; v0];
term_manifold = geo_terminal('manifold', p, []);
warmOpts = struct('par',p, 'mode','mintime', 'rv0',rv0, 'maxIter',3000, ...
                   'printLevel',3, 'warmTight',true);
isGood = @(o) o.success && o.maxDefect < 1e-8 && o.termErr < 1e-8;

if ~isempty(routeB)
    fprintf('MINTIME ROUTE-B: energy warm-start @ T=%g N (sbar=%.3f, tDur=%.2f ND, N=%d)\n', ...
            thrustN, routeB.sbar, routeB.tDur, N);
    [sgE, X0E, U0E, tauf0E, ~] = seed_2body(p, rv0, struct('sbar',routeB.sbar,'tDur',routeB.tDur,'N',N));
    energyOpts = struct('par',p, 'mode','fixedtf', 'eps',1, 'tfTarget',routeB.tDur, ...
                         'rv0',rv0, 'maxIter',3000, 'warmTight',false, 'printLevel',3);
    outE = casadi_lt_2body(sgE, X0E, U0E, tauf0E, term_manifold, energyOpts);
    fprintf('  energy stage: status=%s defect=%.3e tf=%.4f (target %.4f)\n', ...
            outE.ipoptStatus, outE.maxDefect, outE.tf, routeB.tDur);
    if ~(outE.success && outE.maxDefect < 1e-8)
        error('run_mintime:routeBEnergyFail', ...
            'Route-B energy stage failed for thrustN=%g N: status=%s defect=%.3e', ...
            thrustN, outE.ipoptStatus, outE.maxDefect);
    end
    X0mt = outE.X;  U0mt = outE.U;  U0mt(4,:) = 1;   % s=1 warm-start row ('mintime' pins s==1 anyway)
    mtOpts = warmOpts;  mtOpts.warmTight = false;    % loose first pass from the energy shape
    outMT = casadi_lt_2body(sgE, X0mt, U0mt, tauf0E, term_manifold, mtOpts);
    fprintf('  mintime (from energy): status=%s defect=%.3e termErr=%.3e\n', ...
            outMT.ipoptStatus, outMT.maxDefect, outMT.termErr);
    [out, round_] = mintime_continue_only(sgE, tauf0E, term_manifold, warmOpts, outMT, isGood, ...
        CONTIN_MAX_ROUNDS, CONTIN_DECADE_MIN, thrustN, hx0, 'Route-B');
    sinfo = struct('routeB', true, 'sbar', routeB.sbar, 'tDur', routeB.tDur, ...
                    'energyDefect', outE.maxDefect);
elseif ~isempty(seedAnchor)
    A = load(seedAnchor);
    Aout = A.res.out;
    Nanchor = size(Aout.X, 2) - 1;                % anchor's own mesh (decoupled)
    sgA = linspace(0, 1, Nanchor+1).';
    stretch = A.res.thrustN / thrustN;            % C-law: T*tf ~ const
    X0A = Aout.X;  X0A(8,:) = Aout.X(8,:) * stretch;   % only the time row is rescaled
    U0A = Aout.U;  tauf0A = Aout.tauf0;
    LendA = atan2(X0A(2,end), X0A(1,end));        % positions unchanged by the stretch
    term_fixedA = geo_terminal('fixed', p, LendA);
    baseOptsA = warmOpts;  baseOptsA.warmTight = false;
    fprintf('MINTIME thrust-continuation seed from %s (T=%g N -> T=%g N, stretch=%.3f)\n', ...
            seedAnchor, A.res.thrustN, thrustN, stretch);
    try
        [out, round_] = mintime_stage12_continue(sgA, X0A, U0A, tauf0A, term_fixedA, ...
            term_manifold, baseOptsA, warmOpts, isGood, CONTIN_MAX_ROUNDS, CONTIN_DECADE_MIN, ...
            thrustN, hx0, sprintf('anchor-seeded @ L=%.4f', LendA));
        N = Nanchor;
        sinfo = struct('continuationFrom', seedAnchor, 'anchorThrustN', A.res.thrustN, ...
                        'stretch', stretch, 'fallback', false);
    catch ME_anchor
        fprintf(['ANCHOR-SEEDED continuation FAILED (%s)\n  falling back to a FRESH cold ' ...
                 'tangential seed at N=1200 (nodes-per-rev parity)...\n'], ME_anchor.message);
        N = 1200;
        [sg, X0, U0, tauf0, sinf] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',N));
        term_fixed = geo_terminal('fixed', p, sinf.Larr);
        baseOpts = warmOpts;  baseOpts.warmTight = false;
        [out, round_] = mintime_stage12_continue(sg, X0, U0, tauf0, term_fixed, term_manifold, ...
            baseOpts, warmOpts, isGood, CONTIN_MAX_ROUNDS, CONTIN_DECADE_MIN, thrustN, hx0, ...
            'fallback cold-seed N=1200');
        sinfo = struct('continuationFrom', seedAnchor, 'fallback', true, 'fallbackN', N, ...
                        'anchorAttemptError', ME_anchor.message);
    end
else
    [sg, X0, U0, tauf0, sinf] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',N));
    term_fixed = geo_terminal('fixed', p, sinf.Larr);
    baseOpts = warmOpts;  baseOpts.warmTight = false;
    [out, round_] = mintime_stage12_continue(sg, X0, U0, tauf0, term_fixed, term_manifold, ...
        baseOpts, warmOpts, isGood, CONTIN_MAX_ROUNDS, CONTIN_DECADE_MIN, thrustN, hx0, 'cold-seed');
    sinfo = sinf;
end

Lun = unwrap(atan2(out.X(2,:), out.X(1,:)));
res = struct('out', out, 'tfmin', out.tf, 'tfmin_h', out.tf*p.TU_s/3600, ...
             'dL_mt', Lun(end)-Lun(1), 'revs', (Lun(end)-Lun(1))/(2*pi), ...
             'thrustN', thrustN, 'hx0', hx0, 'N', N, 'seedInfo', sinfo, ...
             'continuationRounds', round_);
save(fn, 'res');
fprintf('MINTIME T=%g N: tf=%.4f ND = %.1f h, dL=%.1f rad (%.2f revs), defect %.2e, %s\n', ...
        thrustN, res.tfmin, res.tfmin_h, res.dL_mt, res.revs, out.maxDefect, out.ipoptStatus);
end

% ---------------------------------------------------------------------------
function [out, round_] = mintime_stage12_continue(sg, X0, U0, tauf0, term_fixed, ...
    term_manifold, baseOpts, warmOpts, isGood, roundsMax, decadeMin, thrustN, hx0, label)
% MINTIME_STAGE12_CONTINUE  Shared stage-1/stage-2/continuation recipe.
%
% Runs the fixed-rendezvous warm-up (stage 1, max-iter outcome accepted),
% then the free-longitude manifold solve (stage 2, warmTight), then up to
% roundsMax further warm continuation rounds against the manifold, stopping
% early (error) on a truly flat round (< decadeMin decades of defect
% improvement) and erroring if the gate (isGood) is still unmet after
% roundsMax rounds. Used identically by the cold-seed and anchor-seeded
% (and its N=1200 fallback) paths in RUN_MINTIME, with different
% roundsMax/decadeMin calibration per path.
%
% INPUTS:  sg [(N+1)x1] tau grid;  X0/U0 warm-start state/control;  tauf0
%          [scalar];  term_fixed/term_manifold [struct, see GEO_TERMINAL];
%          baseOpts/warmOpts [struct, see CASADI_LT_2BODY opts, stage-1/2];
%          isGood [function handle];  roundsMax/decadeMin [scalar];
%          thrustN/hx0 [scalar, for error messages];  label [char, for
%          fprintf/error messages]
% OUTPUTS: out - converged (or best-effort) CASADI_LT_2BODY output struct;
%          round_ - number of continuation rounds actually run
%
% REFERENCES: [1] Task 14 controller triage round 2 (factored out of the
%   former inline cold-seed-only recipe to share with anchor-seeded paths).
fprintf('MINTIME %s stage 1 (fixed rendezvous)...\n', label);
out1 = casadi_lt_2body(sg, X0, U0, tauf0, term_fixed, baseOpts);
fprintf('  stage1: status=%s defect=%.3e termErr=%.3e (max-iter here is fine -- warm-up only)\n', ...
        out1.ipoptStatus, out1.maxDefect, out1.termErr);

fprintf('MINTIME %s stage 2 (manifold, warm-started from stage 1)...\n', label);
out = casadi_lt_2body(sg, out1.X, out1.U, tauf0, term_manifold, warmOpts);
fprintf('  stage2: status=%s defect=%.3e termErr=%.3e\n', out.ipoptStatus, out.maxDefect, out.termErr);

[out, round_] = mintime_continue_only(sg, tauf0, term_manifold, warmOpts, out, isGood, ...
    roundsMax, decadeMin, thrustN, hx0, label);
end

% ---------------------------------------------------------------------------
function [out, round_] = mintime_continue_only(sg, tauf0, term_manifold, warmOpts, out, ...
    isGood, roundsMax, decadeMin, thrustN, hx0, label)
% MINTIME_CONTINUE_ONLY  Shared warm-continuation loop + stall/converge guard.
%
% Given an ALREADY-COMPUTED 'out' (from stage 2, or from a Route-B energy-
% warm-started mintime solve), repeats the manifold solve warm-started from
% the previous iterate (warmTight=true, per warmOpts) up to roundsMax times,
% stopping early (error) on a truly flat round (< decadeMin decades of defect
% improvement) and erroring if the gate (isGood) is still unmet after
% roundsMax rounds. Factored out of the former MINTIME_STAGE12_CONTINUE
% (round 5) so the Route-B path can reuse the identical, already-validated
% guard (see test_stall_guard.m) without going through a stage-1/stage-2
% warm-up it does not need.
%
% INPUTS:  sg [(N+1)x1] tau grid;  tauf0 [scalar];  term_manifold [struct,
%          see GEO_TERMINAL];  warmOpts [struct, see CASADI_LT_2BODY opts];
%          out - CASADI_LT_2BODY output to continue from;  isGood [function
%          handle];  roundsMax/decadeMin [scalar];  thrustN/hx0 [scalar, for
%          error messages];  label [char, for fprintf/error messages]
% OUTPUTS: out - converged (or best-effort) CASADI_LT_2BODY output struct;
%          round_ - number of continuation rounds actually run
%
% REFERENCES: [1] Task 14 controller triage round 2 (origin, as part of
%   mintime_stage12_continue). [2] Round 5 (factored out for Route-B reuse).
round_ = 0;
defect0 = out.maxDefect;                          % for the cumulative progress print
gateDefect = 1e-8;
while ~isGood(out) && round_ < roundsMax
    round_ = round_ + 1;
    prevDefect = out.maxDefect;
    outNew = casadi_lt_2body(sg, out.X, out.U, tauf0, term_manifold, warmOpts);
    fprintf('  continuation round %d: defect %.3e -> %.3e, termErr=%.3e, status=%s\n', ...
            round_, prevDefect, outNew.maxDefect, outNew.termErr, outNew.ipoptStatus);
    decadeImprove   = log10(max(prevDefect, realmin)) - log10(max(outNew.maxDefect, realmin));
    cumDecades      = log10(max(defect0, realmin)) - log10(max(outNew.maxDefect, realmin));
    decadesRemaining = log10(max(outNew.maxDefect, realmin)) - log10(gateDefect);
    fprintf(['    cumulative: %.2f decades closed since start (defect0=%.3e), ' ...
             '~%.2f decades remaining to gate (%d/%d rounds used)\n'], ...
            cumDecades, defect0, max(decadesRemaining, 0), round_, roundsMax);
    if ~isGood(outNew) && decadeImprove < decadeMin
        error('run_mintime:stall', ['continuation stalled at round %d for thrustN=%g hx0=%g ' ...
            '(%s): defect %.3e -> %.3e (%.2f decades, need >=%.2f), termErr=%.3e, status=%s'], ...
            round_, thrustN, hx0, label, prevDefect, outNew.maxDefect, decadeImprove, ...
            decadeMin, outNew.termErr, outNew.ipoptStatus);
    end
    out = outNew;
end
if ~isGood(out)
    error('run_mintime:noConverge', ...
        ['thrustN=%g hx0=%g (%s) failed to converge after %d continuation round(s): ' ...
         'defect=%.3e termErr=%.3e status=%s'], ...
        thrustN, hx0, label, round_, out.maxDefect, out.termErr, out.ipoptStatus);
end
end
