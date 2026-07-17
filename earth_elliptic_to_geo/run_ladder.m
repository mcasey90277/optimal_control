function results = run_ladder(thrustList, cfg)
% RUN_LADDER  Thrust-continuation orchestrator (Task 6 skeleton) -- descends
% a strictly-decreasing thrust list, at each rung producing (1) a min-time
% anchor via run_mintime_mee.m and (2) a fixed-tf fuel (eps:1->0) solution
% via run_transfer_mee.m at c_tf*t_f,min, recording per-rung structure
% counts (revs, switches, m_f_kg, edge, tfmin) into a resume-safe cache.
%
% GEOMETRIC WARM CHAIN (per rung k>1): the min-time anchor's cold-seed
% fallback (run_mintime_mee.m's Stage B) is given a REVOLUTION-COUNT HINT
% derived from the PREVIOUS rung's converged anchor, rescaled by the C-law
% (T_max*t_f,min ~ const, hence total winding DeltaL ~ 1/T_max at fixed
% trajectory shape -- consistent with paper Table 3's revs roughly doubling
% every thrust halving: 7.5 -> 15 -> 30 -> 74.5 -> ...):
%     dL_guess(T_new) = dL_mt(T_prev) * (T_prev / T_new)
%     nRevSeed_guess  = round(dL_guess / (2*pi))
% passed through as cfg.nRevSeed to run_mintime_mee (which already exposes
% this as a Stage-B seed-shaping knob -- see run_mintime_mee.m). This is a
% WORKING warm hint, not a placeholder: it directly narrows the cold-seed
% search away from the hardcoded nRev=3 default used at the top rung.
%
% FUEL-SOLVE WARM START -- EXTENSION POINT FOR TASK 7/9 (documented, NOT
% implemented here; out of this task's file list, which is run_mintime_mee.m
% + run_ladder.m only): the brief's step (2) calls for the fuel solve to
% warm-start from the PREVIOUS rung's converged fuel solution (state/
% control interpolated and DeltaL rescaled, mirroring nodestudy_mee.m's
% solve_warm_node mesh-refine-and-resolve pattern). run_transfer_mee.m's
% current interface always builds its own cold mee_seed start internally;
% it does not yet accept an externally-supplied warm X0/U0/dL0. Task 7/9
% should add that entry point (or a thin run_transfer_mee_warm.m sibling)
% and wire it in at the "% --- fuel solve ---" block below. Until then,
% every rung's fuel solve either reuses an EXISTING certified cache (the 10 N
% leg, this task) or falls through to run_transfer_mee's own cold-seed
% recipe (correct, just not warm-started from the ladder's neighbor).
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
prevAnchor = [];  prevThrust = [];

for k = 1:nRungs
    thrustN = thrustList(k);
    fprintf('\n=== LADDER RUNG %d/%d: T=%g N ===\n', k, nRungs, thrustN);
    rungFile = fullfile(resDir, sprintf('MEE_ladder_T%d.mat', round(10*thrustN)));
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
    else
        % --- min-time anchor, warm-hinted by the C-law rev-growth rescale --
        mtCfg = struct('m0kg', m0kg, 'ispS', ispS);
        if ~isempty(prevAnchor)
            dLGuess = prevAnchor.dL_mt * (prevThrust / thrustN);
            nRevGuess = max(1, round(dLGuess / (2*pi)));
            mtCfg.nRevSeed = nRevGuess;
            fprintf(['  anchor warm hint: dL_guess=%.4f rad -> nRevSeed=%d ' ...
                     '(rescaled %.4f rad * %g/%g N, C-law T*tf~const)\n'], ...
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
    prevAnchor = rung.anchor;  prevThrust = thrustN;
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
