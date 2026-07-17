function results = run_ladder(thrustList, cfg)
% RUN_LADDER  Thrust-continuation orchestrator (Task 6 skeleton) -- descends
% a strictly-decreasing thrust list, at each rung producing (1) a min-time
% anchor via run_mintime_mee.m and (2) a fixed-tf fuel (eps:1->0) solution
% via run_transfer_mee.m at c_tf*t_f,min, recording per-rung structure
% counts (revs, switches, m_f_kg, edge, tfmin) into a resume-safe cache.
%
% GEOMETRIC WARM CHAIN (per rung k>1): the min-time anchor's Stage B is
% given the PREVIOUS rung's converged anchor as a FULL-TRAJECTORY warm
% start (X/U/dL, mesh-refined by run_mintime_mee.m -- see its ANCHOR
% WARM-START FIX), DeltaL rescaled by the C-law (T_max*t_f,min ~ const,
% hence total winding DeltaL ~ 1/T_max at fixed trajectory shape --
% consistent with paper Table 3's revs roughly doubling every thrust
% halving: 7.5 -> 15 -> 30 -> 74.5 -> ...):
%     dL_guess(T_new) = dL_mt(T_prev) * (T_prev / T_new)
% passed through as cfg.warmStartAnchor.dL (plus .X/.U/.N from the prior
% anchor) to run_mintime_mee. A live 5 N probe (Task 7) found the ORIGINAL
% design -- a raw cold tangential seed, only hinted via nRevSeed (an
% integer target passed to mee_seed) -- to be a genuinely hard NLP start
% once nRevSeed exceeds ~3-4 (dual infeasibility diverges, hard stall,
% independent of seed throttle conditioning); the full-trajectory warm
% start closes that gap the same way the fuel solve's warm start does.
% cfg.nRevSeed is still passed too (used only as Stage B's cold-seed
% FALLBACK inside run_mintime_mee.m when no warmStartAnchor is available,
% e.g. a future direct/standalone call).
%
% FUEL-SOLVE WARM START (Task 7 -- CLOSES the extension point Task 6 left
% documented-but-unimplemented): run_transfer_mee.m now accepts cfg.warmStart
% (struct .sigma/.X/.U/.dL, see its header). Every rung k>1 that needs a
% FRESH fuel solve (i.e. no existing certified cache under mee_fuel_tag) is
% given the PREVIOUS rung's converged fuel trajectory (X/U/dL/sigma, kept
% in `prevFuelFull` across loop iterations -- sourced from `fuelRes` whether
% that rung was freshly solved OR loaded from its own <fuelTag>.mat cache,
% so this warm chain survives a killed-and-resumed run_ladder invocation
% with no extra persistence: each rung's <fuelTag>.mat already carries the
% full trajectory), with DeltaL C-law rescaled exactly like the anchor's
% own warm hint above:
%     dL_guess(T_new) = dL_prevFuel * (T_prev / T_new)
% passed as cfg.warmStart.dL; run_transfer_mee.m sizes its own node count
% and interpolation from this. The 10 N leg (k=1, no prevFuelFull) is
% unaffected -- reuses its existing certified artifact exactly as in Task 6.
%
% RESUME-SAFE, PER-RUNG CACHING: each rung's result is cached to
% resDir/MEE_ladder_T<10*thrustN>.mat (keyed by THRUST, not list position,
% so re-invoking with a different-length/reordered thrustList never
% mis-attributes a cached rung); a rung whose cache file already exists is
% loaded, not re-solved -- this is what makes the 10 N leg a genuine
% NO-RE-SOLVE reuse of Task 4/6's existing certified artifacts (both the
% mintime anchor MEE_mintime_T100.mat and the fuel solution MEE_M2_10N.mat
% predate this file and are loaded as-is).
%
% INPUTS:  thrustList - thrust levels [N], STRICTLY DESCENDING (each rung
%          warm-hints the next lower one), e.g. [10 5 2.5 1]
%          cfg - optional struct: .ctf (1.5), .nodesPerRev (25), .m0kg
%          (1500), .ispS (2000), .maxIter (1500, fuel homotopy), .seedThr
%          (0.4, fuel cold-seed throttle), .betaMode ('tangential')
% OUTPUTS: results - struct array [1 x numel(thrustList)], fields
%          .thrustN .anchor (run_mintime_mee out struct) .fuelTag .fuel
%          (run_transfer_mee report struct) .tf .certified .reused (true
%          if BOTH the anchor and the fuel solve were loaded from a
%          pre-existing cache -- no solve at all for that rung)
%
% REFERENCES: [1] DESIGN_thrust_ladder.md sec 2 "Phase 2" (thrust-
%   continuation backbone, C-law rev-growth rationale).
%   [2] PLAN_thrust_ladder.md Task 6/7 (this file's scope: orchestration +
%       10 N leg here, the real 5/2.5/1 N descent in Task 7).
%   [3] earth_elliptic_to_geo/run_mintime_mee.m (per-rung anchor).
%   [4] earth_elliptic_to_geo/run_transfer_mee.m (per-rung fuel solve).
%   [5] earth_elliptic_to_geo/mee_fuel_tag.m (shared fuel-tag convention).
if nargin < 2, cfg = struct(); end
d = @(f,v) getdef_ladder(cfg, f, v);

ctf         = d('ctf', 1.5);
nodesPerRev = d('nodesPerRev', 25);
m0kg        = d('m0kg', 1500);
ispS        = d('ispS', 2000);
maxIter     = d('maxIter', 1500);
seedThr     = d('seedThr', 0.4);
betaMode    = d('betaMode', 'tangential');
% mtMaxIter -- OPERATIONAL checkpoint-granularity knob (Task 7, execution-
% environment constraint, NOT a solver-quality change): run_mintime_mee's
% own per-round cache only saves AFTER a full casadi_lt_mee call returns; at
% large N (2.5/1 N rungs) a single MUMPS-dominated Newton solve can run
% longer than one execution window with zero intermediate checkpoint. A
% capped mtCfg.maxIter turns "one long uninterruptible solve" into several
% cheaper, independently-cached continuation rounds via the EXISTING
% retain-if-improved machinery -- same final converged answer (still gated
% at defect<1e-8/termErr<1e-8), just reached in more, individually-
% resumable steps. Default 3000 leaves single-invocation callers unaffected.
%
% mtMaxIter is DELIBERATELY EXCLUDED from the per-rung fingerprint `fp`
% built below (review finding, Fix 4): it only changes HOW MANY per-round
% checkpoints the continuation loop takes to reach the same converged
% answer (checkpoint granularity only), not the physics/config the rung
% actually solves -- so two invocations differing only in mtMaxIter must be
% treated as the SAME rung (cache hit, no re-solve), not flagged as a
% fingerprint mismatch. Every other field in `fp` IS solver-relevant and
% must stay fingerprinted.
mtMaxIter   = d('mtMaxIter', 3000);

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

thrustList = thrustList(:).';
assert(numel(thrustList) >= 1, 'run_ladder:empty', 'thrustList must have at least one thrust');
assert(all(diff(thrustList) < 0), 'run_ladder:notDescending', ...
    'thrustList must be STRICTLY DESCENDING (each rung warm-hints the next lower one): got %s', ...
    mat2str(thrustList));

nRungs  = numel(thrustList);
results = repmat(empty_rung(), 1, nRungs);
prevAnchor = [];  prevThrust = [];  prevFuelFull = [];

for k = 1:nRungs
    thrustN = thrustList(k);
    fprintf('\n=== LADDER RUNG %d/%d: T=%g N ===\n', k, nRungs, thrustN);
    rungFile = fullfile(resDir, sprintf('MEE_ladder_T%d.mat', round(10*thrustN)));
    % NOTE: mtMaxIter is intentionally NOT a field here -- see its
    % declaration above (checkpoint granularity only, not a config change).
    fp = struct('thrustN', thrustN, 'ctf', ctf, 'nodesPerRev', nodesPerRev, ...
        'm0kg', m0kg, 'ispS', ispS, 'maxIter', maxIter, 'seedThr', seedThr, ...
        'betaMode', betaMode);

    if isfile(rungFile)
        S = load(rungFile);
        check_cache_fp_ladder(S, fp, rungFile, thrustN);
        rung = S.rung;
        fprintf('  [cached rung] %s: tfmin=%.4f ND revs_mt=%.3f | fuel mf=%.2f kg sw=%d revs=%.3f\n', ...
                rungFile, rung.anchor.tfmin, rung.anchor.revs, rung.fuel.m_f_kg, ...
                rung.fuel.switches, rung.fuel.revs);
        % Reload the full fuel trajectory (X/U/dL/sigma) from its own
        % <fuelTag>.mat cache -- the rung cache only stores the SUMMARY
        % report (rung.fuel), not the full state, but the next rung's
        % warm-start chain (below) needs the full trajectory regardless of
        % whether THIS rung was freshly solved or loaded from cache.
        Sfc = load(fullfile(resDir, [rung.fuelTag '.mat']));
        fuelRes = Sfc.res;
    else
        % --- min-time anchor, warm-hinted by the C-law rev-growth rescale --
        mtCfg = struct('m0kg', m0kg, 'ispS', ispS, 'maxIter', mtMaxIter);
        if ~isempty(prevAnchor)
            dLGuess = prevAnchor.dL_mt * (prevThrust / thrustN);
            nRevGuess = max(1, round(dLGuess / (2*pi)));
            mtCfg.nRevSeed = nRevGuess;
            % FULL-TRAJECTORY anchor warm start (Task 7 addition -- see
            % run_mintime_mee.m's ANCHOR WARM-START FIX): a raw cold
            % tangential seed turned out to be a genuinely hard NLP start
            % once the C-law hint pushes nRevSeed beyond ~3-4 (offline
            % probe at the 5 N rung: dual infeasibility exploding past 1e6,
            % defect stuck ~1.3e-2 across multiple full continuation
            % rounds -- a real stall, not a slow-but-convergent solve).
            % The PREVIOUS rung's own converged min-time anchor (already
            % thr=1 throughout, already exactly satisfying the terminal
            % manifold) is a far better-conditioned starting shape than a
            % constant-throttle spiral -- same self-similar-shape argument
            % as the fuel warm start, applied here to close the SAME class
            % of gap for the anchor.
            mtCfg.warmStartAnchor = struct('X', prevAnchor.solverOut.X, ...
                'U', prevAnchor.solverOut.U, 'dL', dLGuess, 'N', prevAnchor.N);
            fprintf(['  anchor warm hint: dL_guess=%.4f rad -> nRevSeed=%d ' ...
                     '(rescaled %.4f rad * %g/%g N, C-law T*tf~const; full-trajectory ' ...
                     'warm start from prior rung''s converged anchor)\n'], ...
                    dLGuess, nRevGuess, prevAnchor.dL_mt, prevThrust, thrustN);
        end
        anchorOut = run_mintime_mee(thrustN, nodesPerRev, mtCfg);

        % --- fuel solve: reuse an existing certified cache verbatim (no
        % re-solve, no fingerprint re-derivation against it -- this IS the
        % "10 N leg reuses existing artifacts" requirement), else produce a
        % fresh cold-seed fuel solve at this rung's own anchor ------------
        fuelTag  = mee_fuel_tag(thrustN);
        fuelFile = fullfile(resDir, [fuelTag '.mat']);
        if isfile(fuelFile)
            Sf = load(fuelFile);
            fuelRes = Sf.res;
            reused = true;
            fprintf('  fuel solve: REUSED existing certified artifact %s (no re-solve)\n', fuelFile);
        else
            fuelCfg = struct('thrustN', thrustN, 'ctf', ctf, 'tfMinAnchor', anchorOut.tfmin, ...
                'tag', fuelTag, 'seedThr', seedThr, 'betaMode', betaMode, ...
                'nodesPerRev', nodesPerRev, 'maxIter', maxIter, 'm0kg', m0kg, 'ispS', ispS);
            if ~isempty(prevFuelFull)
                dLGuessFuel = prevFuelFull.fuel.dL * (prevThrust / thrustN);
                fuelCfg.warmStart = struct('sigma', prevFuelFull.sigma, ...
                    'X', prevFuelFull.fuel.X, 'U', prevFuelFull.fuel.U, 'dL', dLGuessFuel);
                fprintf(['  fuel solve warm hint: dL_guess=%.4f rad -> revs_guess=%.3f ' ...
                         '(rescaled %.4f rad * %g/%g N, prior rung''s converged fuel dL)\n'], ...
                        dLGuessFuel, dLGuessFuel/(2*pi), prevFuelFull.fuel.dL, prevThrust, thrustN);
            end
            fuelRes = run_transfer_mee(fuelCfg);
            reused = false;
            assert(fuelRes.report.certified, 'run_ladder:fuelUncertified', ...
                'T=%g N: fuel solve did not certify (defect=%.2e) -- rung BLOCKED', ...
                thrustN, fuelRes.report.defect);
        end

        rung = struct('thrustN', thrustN, 'anchor', anchorOut, 'fuelTag', fuelTag, ...
            'fuel', fuelRes.report, 'tf', fuelRes.tf, 'certified', fuelRes.report.certified, ...
            'reused', reused, 'fp', fp);
        save(rungFile, 'rung');
        fprintf('  rung DONE: anchor tfmin=%.4f ND, fuel mf=%.2f kg sw=%d revs=%.3f edge=%.1f%%\n', ...
                anchorOut.tfmin, rung.fuel.m_f_kg, rung.fuel.switches, rung.fuel.revs, ...
                100*rung.fuel.edge);
    end

    results(k) = rung;
    prevAnchor = rung.anchor;  prevThrust = thrustN;  prevFuelFull = fuelRes;
end

save(fullfile(resDir, 'MEE_ladder.mat'), 'results');
print_ladder_table(results);
end

% ---------------------------------------------------------------------------
function print_ladder_table(results)
% PRINT_LADDER_TABLE  Fixed-width console summary + law-R0 spread (T*tfmin).
fprintf('\n--- LADDER SUMMARY ---\n');
fprintf('%-8s %-10s %-10s %-6s %-10s %-8s %-8s %-8s\n', ...
    'T [N]', 'tfmin[ND]', 'tfmin[h]', 'revsMT', 'mf[kg]', 'sw', 'revs', 'reused');
C = zeros(1, numel(results));
for k = 1:numel(results)
    r = results(k);
    fprintf('%-8g %-10.4f %-10.2f %-6.3f %-10.2f %-8d %-8.3f %-8d\n', ...
        r.thrustN, r.anchor.tfmin, r.anchor.tfmin_h, r.anchor.revs, ...
        r.fuel.m_f_kg, r.fuel.switches, r.fuel.revs, r.reused);
    C(k) = r.thrustN * r.anchor.tfmin_h;
end
if numel(C) > 1
    fprintf('law R0: T*tfmin = %s N.h (spread %.1f%%, paper ~850, Cartesian 846.6)\n', ...
            mat2str(round(C, 1)), 100*(max(C)-min(C))/mean(C));
else
    fprintf('law R0: T*tfmin = %.2f N.h (single rung, no spread; paper ~850, Cartesian 846.6)\n', C(1));
end
end

% ---------------------------------------------------------------------------
function r = empty_rung()
% EMPTY_RUNG  Placeholder struct matching the real rung's field set (used
% only for repmat preallocation; every entry is overwritten in the loop).
r = struct('thrustN', 0, 'anchor', struct('tfmin', 0, 'tfmin_h', 0, 'dL_mt', 0, 'revs', 0), ...
    'fuelTag', '', 'fuel', struct('m_f_kg', 0, 'switches', 0, 'revs', 0, 'edge', 0, 'certified', false), ...
    'tf', 0, 'certified', false, 'reused', false, 'fp', struct());
end

% ---------------------------------------------------------------------------
function check_cache_fp_ladder(S, fp, file, thrustN)
% CHECK_CACHE_FP_LADDER  Fail-loud cache-fingerprint guard (mirrors
% run_transfer_mee.m's check_cache_fp / run_mintime_mee.m's check_cache_fp_mt).
% BACKWARD COMPAT: a pre-fix cache with no .fp field only WARNs.
if ~isfield(S, 'rung') || ~isfield(S.rung, 'fp')
    warning('run_ladder:noCachedFingerprint', ['rung cache for T=%g N has no ' ...
        'stored config fingerprint -- trusting it as-is'], thrustN);
    return;
end
sfp = S.rung.fp;
flds = fieldnames(fp);
for kf = 1:numel(flds)
    f = flds{kf};
    if ~isfield(sfp, f) || ~isequal(sfp.(f), fp.(f))
        error('run_ladder:fingerprintMismatch', ['cached config fingerprint ' ...
            'mismatch in %s: field ''%s'' differs between the cache and the ' ...
            'current config -- stale rung cache for T=%g N; delete the file ' ...
            'or change cfg'], file, f, thrustN);
    end
end
end

% ---------------------------------------------------------------------------
function v = getdef_ladder(s, f, dflt)
% GETDEF_LADDER  Optional-field default (mirrors casadi_lt_2body's helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
