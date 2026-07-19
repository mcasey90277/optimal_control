function row = reproduce_row(T, opts)
% REPRODUCE_ROW  Table-3 reproducer ENGINE: re-solve one thrust rung's
% min-fuel Earth-elliptic-to-GEO transfer FROM SCRATCH, verify it against
% the campaign's certified numbers, and cache the result.
%
% Composes the existing certified drivers (run_mintime_mee.m,
% run_transfer_mee.m, psr_mee_refine.m, anchor_smallN_first.m) under a
% strategy dispatched from table3_recipes.m (the anchor stage's strategy:
% coldB/chain/R0law/smallN_first -- the last, the 1 N low-node-density-first
% strategy, is anchor_smallN_first.m, Task 3). Every solve uses REPRO_-prefixed
% tags so a from-scratch run can
% NEVER load or clobber the campaign's own production caches under the
% plain MEE_mintime_T.../MEE_M2_...N tags -- this is the whole point of
% the "engine": prove the certified numbers are reproducible by an
% independent cold build, not merely re-read a cached answer.
%
% FROM-SCRATCH ISOLATION: cfg.tag/.fuelTag for every driver call below are
% REPRO_-prefixed (see ttag/local helper), landing in the SAME
% earth_elliptic_to_geo/results/ directory the campaign uses (run_mintime_mee/
% run_transfer_mee fix their own resDir internally; psr_mee_refine does accept
% an outDir, which the PSR stage points at results/repro/) but under names that
% can never collide with a campaign tag. The anchor stage's
% cfg.fuelTag is set to a tag that is GUARANTEED never to exist
% ('REPRO_none_<ttag>'), forcing run_mintime_mee.m's Stage A (warm-start
% from an existing certified fuel anchor) to be skipped unconditionally,
% so Stage B is a genuinely cold solve, not a reuse of Task 6/7's own
% campaign anchor. This engine's OWN row summary (the .mat this function
% itself writes) is kept in a separate results/repro/ subdirectory, so it
% is visually and physically distinct from both the campaign caches and
% the driver-internal REPRO_ cache files.
%
% opts.reuseCampaignCache (default false) is plumbed through but NOT
% exercised by default: when true, the driver tags fall back to the
% campaign's own canonical tags (mee_fuel_tag.m / run_mintime_mee.m's own
% default 'MEE_mintime_T<10*thrustN>') instead of the REPRO_ tags, so a
% caller who explicitly opts in can let the drivers' own file-based
% caching short-circuit an already-solved rung (useful for fast iteration
% on THIS function without re-paying a multi-minute solve) -- this
% deliberately breaks the from-scratch guarantee and must never be the
% default.
%
% THROW, NEVER FAKE: if the fuel stage (run_transfer_mee) or the optional
% PSR refinement pass (psr_mee_refine) comes back uncertified, this
% function errors out immediately -- a row is never assembled from an
% uncertified solve. The anchor stage (run_mintime_mee) already enforces
% the same discipline internally (it either returns a certified anchor or
% throws), so no additional anchor-side guard is required here beyond a
% defensive assert.
%
% INPUTS:
%   T    - max thrust level [N]; one of 10, 5, 2.5 (coldB/chain), 1
%          (smallN_first -- anchor_smallN_first.m, Task 3), 0.5 (R0law;
%          its fuel-stage chain predecessor is the 1 N smallN_first rung) [scalar]
%   opts - optional struct:
%     .reuseCampaignCache - if true, use the campaign's own canonical
%                    driver tags instead of REPRO_ tags (see above);
%                    default false (from-scratch isolation preserved)  [logical]
%     .m0kg        - initial mass [kg]; default 1500                  [scalar]
%     .ispS        - specific impulse [s]; default 2000                [scalar]
%     .certifySosc - Task 9, OPT-IN NLP-level SOSC certificate on the
%                    keep-best-mass WINNER only (never per multi-start
%                    candidate -- too costly); default false, so the
%                    reproducer stays fast by default, matching
%                    run_transfer_mee.m's cfg.certifySosc opt-in default.
%                    When true, row.sosc is populated (verify_sosc_mee.m
%                    output) and a FAIL verdict emits a loud warning AND
%                    demotes row.certified to false (via apply_sosc_gate,
%                    the same gate run_transfer_mee.m uses) -- a proven
%                    saddle is never adopted as a clean certified
%                    reproduction, per process/DESIGN_sosc.md sec 11.6.
%                    This does NOT change which candidate won -- mass
%                    alone still decides the keep-best-mass selection;
%                    only the CERTIFIED status of that winner is gated.
%                    PASS/WEAK_MIN/INCONCLUSIVE/ERROR all leave
%                    row.certified true                                [logical]
%
% OUTPUTS:
%   row - gergaud_row.m row struct for this rung (+ .sosc, [] unless
%         opts.certifySosc=true), VERIFIED against
%         table3_certified(T) within defaultTol(T) (verify_row throws on
%         any breach, so a returned row is, by construction, a passing
%         row). Also written to
%         results/repro/REPRO_row_T<round(10*T)>.mat (variables row,
%         anchor, sol, rep).
%
% REFERENCES:
%   [1] .superpowers/sdd/task-2-brief.md (this task's spec + composition
%       pseudocode).
%   [2] earth_elliptic_to_geo/table3_recipes.m / table3_certified.m /
%       verify_row.m (Task 1 foundation this engine consumes).
%   [3] earth_elliptic_to_geo/run_mintime_mee.m / run_transfer_mee.m /
%       psr_mee_refine.m (the existing certified drivers this engine
%       composes; NOT edited by this task).
%   [4] earth_elliptic_to_geo/gergaud_row.m / gergaud_row_str.m (row
%       assembly + Table-3-style formatting, Task 6 of the prior build).
%   [5] earth_elliptic_to_geo/mee_fuel_tag.m (the numeric-tag convention
%       ttag, below, mirrors).

if nargin < 1
    error('reproduce_row:badInput', 'T (thrustN) is required');
end
if nargin < 2, opts = struct(); end
d = @(f,v) optdef(opts, f, v);

reuseCampaignCache = d('reuseCampaignCache', false);
m0kg = d('m0kg', 1500);
ispS = d('ispS', 2000);
certifySosc = d('certifySosc', false);   % Task 9, OPT-IN: see SOSC block below

recipe = table3_recipes(T);
cert   = table3_certified(T);
par    = kepler_lt_params(T, m0kg, ispS);   % dimensional context (TU_s, Tmax); logged below

reproDir = fullfile(module_root(), 'results', 'repro');
if ~exist(reproDir, 'dir'), mkdir(reproDir); end

if reuseCampaignCache
    % Opt-in escape hatch: reuse the campaign's own canonical tags (and
    % therefore its own cached solves) instead of forcing a from-scratch
    % run. NOT the default -- see header note.
    tagMt     = sprintf('MEE_mintime_T%d', round(10*T));   % run_mintime_mee.m's own default
    tagNoFuel = mee_fuel_tag(T);                            % real fuel anchor -> Stage A allowed
    tagFuel   = mee_fuel_tag(T);
else
    tagMt     = sprintf('REPRO_MEE_mintime_%s', ttag(T));
    tagNoFuel = sprintf('REPRO_none_%s', ttag(T));           % never exists -> Stage A forced-skip
    tagFuel   = sprintf('REPRO_MEE_M2_%sN', ttag(T));
end

fprintf(['REPRODUCE_ROW T=%g N: strategy=%s (TU_s=%.2f s, Tmax_ND=%.4e, ' ...
         'reuseCampaignCache=%d)\n'], T, recipe.anchor.strategy, par.TU_s, par.Tmax, ...
        reuseCampaignCache);

% ---------------------------------------------------------------------------
% ANCHOR (strategy dispatch)
% ---------------------------------------------------------------------------
switch recipe.anchor.strategy
    case 'coldB'
        anchor = run_mintime_mee(T, recipe.anchor.npr, struct('tag', tagMt, ...
            'fuelTag', tagNoFuel, 'maxIter', 3000));
        anchor.anchorSource = 'solved';

    case 'chain'
        prev = load_prev(T, recipe.anchor.warmFrom, reproDir);
        dLGuess = prev.anchor.dL_mt * (recipe.anchor.warmFrom / T);
        anchor = run_mintime_mee(T, recipe.anchor.npr, struct('tag', tagMt, ...
            'fuelTag', tagNoFuel, 'maxIter', recipe.anchor.mtMaxIter, ...
            'nRevSeed', max(1, round(dLGuess/(2*pi))), ...
            'warmStartAnchor', struct('X', prev.anchor.solverOut.X, ...
                                       'U', prev.anchor.solverOut.U, ...
                                       'dL', dLGuess, ...
                                       'N', prev.anchor.N)));
        anchor.anchorSource = 'solved';

    case 'R0law'
        anchor = struct('tfmin', recipe.tfmin_or_R0, 'anchorSource', 'R0law', ...
            'thrustN', T, 'dL_mt', [], 'N', [], 'solverOut', []);

    case 'smallN_first'
        prevSN = load_prev(T, recipe.anchor.warmFrom, reproDir);
        aopts  = struct('nprLo', recipe.anchor.nprLo, 'nprHi', recipe.anchor.nprHi, ...
            'tag', tagMt);
        anchor = anchor_smallN_first(T, par, prevSN.anchor, aopts);

    otherwise
        error('reproduce_row:unknownStrategy', ...
            'unrecognized anchor strategy ''%s'' for T=%g N', recipe.anchor.strategy, T);
end

assert(isfield(anchor, 'tfmin') && ~isempty(anchor.tfmin) && ...
       (~isfield(anchor, 'certified') || anchor.certified), ...
       'reproduce_row:anchorNotCertified', ...
       'T=%g N: anchor stage did not produce a certified tfmin -- refusing to proceed', T);

tfMinAnchor = anchor.tfmin;

% ---------------------------------------------------------------------------
% FUEL (keep-best-mass multi-start at fixed tf = ctf*tfMinAnchor; A->B hybrid)
% ---------------------------------------------------------------------------
% The fuel bang-bang basin can be RAZOR-SENSITIVE to tf (see memory
% tenN-minfuel-razor-basin): a single solve at the EXACT min-time anchor may
% land in a slightly-worse local optimum (10 N: 24 sw / 8.118 rev / 1377.086
% kg) while a different seed at the same tf reaches the better one (19 sw /
% ~7.3 rev / up to 1377.19 kg). Since minimum-fuel means MAXIMIZE final mass,
% the fuel stage does a keep-best-mass multi-start rather than trusting one
% solve: (A) a seed set at the exact tf; (B, fallback) a tiny tf-bracket, run
% only if A does not reach the campaign mass reference. See fuel_multistart.
warmStartFuel = [];
if ~isempty(recipe.fuel.warmFrom)
    prevF = load_prev(T, recipe.fuel.warmFrom, reproDir);
    dLGuessFuel = prevF.sol.dL * (recipe.fuel.warmFrom / T);
    warmStartFuel = struct('sigma', prevF.sol.sigma, 'X', prevF.sol.X, ...
        'U', prevF.sol.U, 'dL', dLGuessFuel);
end

res = fuel_multistart(T, tfMinAnchor, recipe, tagFuel, m0kg, ispS, ...
    warmStartFuel, cert, defaultTol(T));

sol = struct('X', res.fuel.X, 'U', res.fuel.U, 'dL', res.fuel.dL, 'sigma', res.sigma);
rep = res.report;

% ---------------------------------------------------------------------------
% PSR (optional post-solution refinement)
% ---------------------------------------------------------------------------
if ~isempty(recipe.psr)
    psrOpts = struct('tag', [tagFuel '_PSR'], 'outDir', reproDir, ...
        'maxRounds', recipe.psr.maxRounds, 'nbr', recipe.psr.nbr, ...
        'maxIter', recipe.fuel.maxIter);
    % Pass globalEvery/globalFactor ONLY when the recipe actually specifies them,
    % so a recipe that omits them keeps psr_mee_refine's OWN default -- do NOT
    % inject 3/1.3 into a rung that intends no global-sweep round (reviewer
    % Important: the 1 N recipe omits these; forcing 3 silently overrode the
    % driver's disabled default).
    if isfield(recipe.psr, 'globalEvery') && ~isempty(recipe.psr.globalEvery)
        psrOpts.globalEvery = recipe.psr.globalEvery;
    end
    if isfield(recipe.psr, 'globalFactor') && ~isempty(recipe.psr.globalFactor)
        psrOpts.globalFactor = recipe.psr.globalFactor;
    end
    psrOut = psr_mee_refine(res, psrOpts);
    if ~psrOut.certified
        error('reproduce_row:psrNotCertified', ...
            ['T=%g N: PSR refinement (psr_mee_refine) did NOT certify its final ' ...
             'measured round (stopReason=%s) -- refusing to assemble a row from an ' ...
             'uncertified solve'], T, psrOut.stopReason);
    end
    sol = struct('X', psrOut.finalOut.X, 'U', psrOut.finalOut.U, ...
        'dL', psrOut.finalOut.dL, 'sigma', psrOut.finalSigma);
    rep = psr_report(psrOut);
end

% ---------------------------------------------------------------------------
% SOSC (Task 9, OPT-IN): certify the WINNING candidate ONCE, not per keep-
% best-mass multi-start candidate -- each verify_sosc_mee call is a ~1-2 min
% warm re-solve + eig inertia, far too costly to run on every fuel_multistart
% seed (mass alone still decides the winner; see fuel_multistart above).
% Gated behind opts.certifySosc (default false) so the reproducer stays fast
% by default -- OPT-IN CHOICE, matching the driver's cfg.certifySosc default.
% When enabled, runs exactly ONCE on the FINAL winner (post-PSR if the
% recipe ran one, since `sol`/`rep` are reassigned to the PSR output above)
% and does NOT change which candidate won: SOSC only ANNOTATES the row via
% row.sosc (process/DESIGN_sosc.md sec 11.6 tiered gate: PASS/WEAK_MIN/
% INCONCLUSIVE/ERROR keep the row as-is; a FAIL verdict is a genuine finding
% -- a best-mass candidate that is a PROVEN SADDLE -- surfaced via a loud
% warning, never by silently swapping in a different candidate).
% Winner-tag provenance: res.cfg.tag is the ACTUAL winning multi-start
% candidate's tag (e.g. '..._A3'/'..._B2', set per-candidate in
% build_fuel_cfg/fuel_multistart above), not the un-suffixed tagFuel stem --
% using it here lets sosc.meta.tag trace to the candidate that actually won
% (reviewer Minor).
winnerTag = tagFuel;
if isfield(res, 'cfg') && isfield(res.cfg, 'tag') && ~isempty(res.cfg.tag)
    winnerTag = res.cfg.tag;
end

sosc = [];
if certifySosc
    savedT = struct('sigma', sol.sigma(:), 'X', sol.X, 'U', sol.U, 'dL', sol.dL, ...
        'tfTarget', 1.5*tfMinAnchor, 'xf', [1;0;0;0;0], 'thrustN', T, ...
        'm0kg', m0kg, 'ispS', ispS, 'maxIter', recipe.fuel.maxIter, ...
        'tag', winnerTag, 'kind', 'MEE_M2');
    sosc = verify_sosc_mee(savedT);
    fprintf('  [sosc] T=%g N: verdict=%s (%s)\n', T, sosc.verdict, sosc.status);
    if strcmp(sosc.verdict, 'FAIL')
        warning('reproduce_row:soscFail', ['T=%g N: the BEST-MASS candidate is a ' ...
            'PROVEN SADDLE (SOSC verdict FAIL) -- %s; keep-best-by-mass selection ' ...
            'is unchanged by design, but the row is demoted to uncertified ' ...
            '(row.certified=false) per DESIGN_sosc.md sec 11.6 -- a proven saddle ' ...
            'is never adopted as a clean certified reproduction'], T, sosc.reason);
    end
end

% ---------------------------------------------------------------------------
% ROW + VERIFY + SAVE
% ---------------------------------------------------------------------------
note = '';
if strcmp(anchor.anchorSource, 'R0law')
    note = '0.5 N: R0-law tfmin estimate (anchor-free)';
end

row = gergaud_row(struct('thrustN', T, 'tfmin_ND', tfMinAnchor, 'ctf', 1.5, ...
    'tf_ND', 1.5*tfMinAnchor, 'm_f_kg', rep.m_f_kg, 'switches', rep.switches, ...
    'revs', rep.revs, 'edge', rep.edge, 'incl_deg', rep.incDeg, 'defect', rep.defect, ...
    'certified', true, 'note', note, 'm0kg', m0kg, 'ispS', ispS, 'sosc', sosc));

% SOSC gate (Task 9 review fix, Important): a FAIL verdict on the winner must
% NEVER be adopted as a clean certified reproduction (DESIGN_sosc.md sec 11.6:
% "the reproducer adopts a candidate whose verdict in {PASS, WEAK_MIN,
% INCONCLUSIVE} (never FAIL)"). Reuse the SAME gate the driver
% (run_transfer_mee.m) uses -- row already carries a .certified field, so
% apply_sosc_gate demotes it to false on FAIL (and re-attaches row.sosc,
% a no-op here since gergaud_row already passed it through). PASS/WEAK_MIN/
% INCONCLUSIVE/ERROR (or certifySosc=false, sosc=[]) leave the row untouched.
if certifySosc
    row = apply_sosc_gate(row, sosc);
end

tol = defaultTol(T);
[~, vinfo] = verify_row(row, cert, tol);   % ONE-SIDED: throws only if WORSE than the campaign floor
if vinfo.improved
    fprintf(['  [verify] T=%g N: m_f=%.4f kg BEATS campaign %.4f by %.4f kg ' ...
             '(structure sw %d vs %d, revs %.3f vs %.3f -- reported, not gated)\n'], ...
        T, row.m_f_kg, cert.m_f_kg, vinfo.improvedKg, vinfo.switches, ...
        vinfo.campaignSwitches, vinfo.revs, vinfo.campaignRevs);
else
    fprintf(['  [verify] T=%g N: m_f=%.4f kg >= campaign floor %.4f (matches; ' ...
             'sw %d/%d revs %.3f/%.3f)\n'], ...
        T, row.m_f_kg, vinfo.massFloor, vinfo.switches, vinfo.campaignSwitches, ...
        vinfo.revs, vinfo.campaignRevs);
end

save(fullfile(reproDir, sprintf('REPRO_row_T%d.mat', round(10*T))), 'row', 'anchor', 'sol', 'rep');

fprintf('%s', gergaud_row_str(row));

end

% ---------------------------------------------------------------------------
function bestRes = fuel_multistart(T, tfMinAnchor, recipe, tagFuel, m0kg, ispS, warmStart, cert, tol)
% FUEL_MULTISTART  Keep-best-mass fuel solve (A->B hybrid) at fixed
% tf = 1.5*tfMinAnchor. Minimum-fuel = maximize final mass, so rather than
% trusting a single solve (which can land in a slightly-worse local optimum
% on the razor-thin fuel basin, memory tenN-minfuel-razor-basin), this tries a
% set of candidates and returns the highest-mass CERTIFIED run_transfer_mee
% output.
%   A: a seed set at the EXACT anchor tf.
%   B (fallback, only if A does not reach the campaign mass reference
%      cert.m_f_kg - tol.m_f_kg): a tiny tf-bracket (a few x1e-5, within
%      c_tf=1.5's numerical precision), keeping the best over A union B.
% A window-rejected seed or transient failure is SKIPPED (never fatal). A
% process-fatal MEX/MUMPS hang cannot be caught in-process; it is absorbed by
% the OUTER per-rung watchdog + run_transfer_mee's per-tag caching (a relaunch
% skips the cached completed seeds and retries the hung one).
%
% INPUTS:  T [scalar] thrust; tfMinAnchor [scalar]; recipe [struct]; tagFuel
%          [char] REPRO fuel tag stem; m0kg/ispS [scalar]; warmStart
%          [struct|[]] prior-rung fuel warm start; cert [struct] table3_certified;
%          tol [struct] verify tolerance (uses .m_f_kg for the A->B gate)
% OUTPUTS: bestRes [struct] the best-mass certified run_transfer_mee output
%
% REFERENCES: [1] memory tenN-minfuel-razor-basin. [2] the reproducer design
%   spec, sec 9b (keep-best-mass multi-start).
ctf = 1.5;
seedSet = fuel_seed_set(recipe, warmStart);
bestRes = [];  bestMf = -inf;

% --- A: seeds at the EXACT anchor tf ---
for s = 1:numel(seedSet)
    cfg = build_fuel_cfg(T, ctf, tfMinAnchor, recipe, m0kg, ispS, warmStart, ...
        seedSet(s), sprintf('%s_A%d', tagFuel, s));
    [bestRes, bestMf] = fuel_try_keep(cfg, bestRes, bestMf);
end

% --- B: tiny tf-bracket, ONLY if A did not reach the campaign mass reference ---
if isempty(bestRes) || bestMf < cert.m_f_kg - tol.m_f_kg
    fracs  = [-5e-5, 5e-5, 2e-4];   % within c_tf=1.5's numerical precision
    bcount = 0;
    for f = fracs
        tb = tfMinAnchor * (1 + f);
        for s = 1:numel(seedSet)
            bcount = bcount + 1;
            cfg = build_fuel_cfg(T, ctf, tb, recipe, m0kg, ispS, warmStart, ...
                seedSet(s), sprintf('%s_B%d', tagFuel, bcount));
            [bestRes, bestMf] = fuel_try_keep(cfg, bestRes, bestMf);
        end
    end
end

if isempty(bestRes)
    error('reproduce_row:fuelAllFailed', ...
        ['T=%g N: no certified fuel solution across the keep-best-mass ' ...
         'multi-start -- refusing to assemble a row'], T);
end
fprintf('  [fuel multi-start] best certified m_f = %.4f kg (%d sw, %.4f rev)\n', ...
    bestRes.report.m_f_kg, bestRes.report.switches, bestRes.report.revs);
end

% ---------------------------------------------------------------------------
function cfg = build_fuel_cfg(T, ctf, tfMinAnchor, recipe, m0kg, ispS, warmStart, seed, tag)
% BUILD_FUEL_CFG  One run_transfer_mee cfg for a given tf + seed + tag.
%
% INPUTS:  T/ctf/tfMinAnchor/m0kg/ispS [scalar]; recipe [struct]; warmStart
%          [struct|[]]; seed [struct .seedThr .betaMode]; tag [char]
% OUTPUTS: cfg [struct] run_transfer_mee config
cfg = struct('thrustN', T, 'ctf', ctf, 'tfMinAnchor', tfMinAnchor, 'tag', tag, ...
    'seedThr', seed.seedThr, 'betaMode', seed.betaMode, 'nodesPerRev', recipe.fuel.npr, ...
    'maxIter', recipe.fuel.maxIter, 'm0kg', m0kg, 'ispS', ispS);
if ~isempty(warmStart), cfg.warmStart = warmStart; end
end

% ---------------------------------------------------------------------------
function [bestRes, bestMf] = fuel_try_keep(cfg, bestRes, bestMf)
% FUEL_TRY_KEEP  Solve one fuel cfg; retain it iff CERTIFIED and higher-mass
% than the incumbent. A window-rejected seed (run_transfer_mee's seed-revs
% guard) or a transient error is caught and SKIPPED so the multi-start
% tolerates it. (An uncatchable MEX/MUMPS hang is handled one level up -- see
% fuel_multistart's header.)
%
% INPUTS:  cfg [struct] run_transfer_mee config; bestRes [struct|[]] incumbent;
%          bestMf [scalar] incumbent final mass [kg]
% OUTPUTS: bestRes/bestMf updated if this candidate certifies and improves
try
    res = run_transfer_mee(cfg);
    if res.report.certified && res.report.m_f_kg > bestMf
        bestRes = res;  bestMf = res.report.m_f_kg;
    end
catch ME
    % out-of-seed-window / transient failure -- skip this candidate, but LOG
    % which one and why: this stage is failure-prone by design (seed-window
    % rejects, transient IPOPT failures), and a silent drop leaves the terminal
    % fuelAllFailed error with no breadcrumbs (reviewer Important).
    fprintf('  [fuel multi-start] SKIP %s: %s\n', cfg.tag, ME.message);
end
end

% ---------------------------------------------------------------------------
function set = fuel_seed_set(recipe, warmStart)
% FUEL_SEED_SET  Candidate (seedThr, betaMode) seeds for the multi-start.
%   Cold-seed rung (no warmStart): span the recipe default plus higher
%   throttles that escape a worse basin, x both beta modes; run_transfer_mee's
%   own seed-revs window guard rejects out-of-window throttles (caught +
%   skipped in fuel_try_keep), so this list can be generous.
%   Warm rung (warmStart given): the warm trajectory IS the seed, so seedThr/
%   betaMode are inert on the warm path -- a SINGLE candidate; the razor-basin
%   variation (if any) then comes from fuel_multistart's tf-bracket (B).
%
% INPUTS:  recipe [struct]; warmStart [struct|[]]
% OUTPUTS: set [1xK struct] fields .seedThr .betaMode
if ~isempty(warmStart)
    set = struct('seedThr', recipe.fuel.seedThr, 'betaMode', 'tangential');
    return;
end
thrs  = unique([recipe.fuel.seedThr, 0.45, 0.50], 'stable');
betas = {'tangential', 'transverse'};
set   = struct('seedThr', {}, 'betaMode', {});
for t = thrs
    for b = 1:numel(betas)
        set(end+1) = struct('seedThr', t, 'betaMode', betas{b}); %#ok<AGROW>
    end
end
end

% ---------------------------------------------------------------------------
function s = ttag(T)
% TTAG  Numeric-only tag stem for a thrust level, matching mee_fuel_tag.m's
% own convention (integer -> plain digits; non-integer -> decimal point
% replaced with 'p'), but WITHOUT that function's 'MEE_M2_'/'N' wrapping --
% this is the bare numeric fragment reproduce_row.m's own REPRO_ tags are
% built from.
%
% INPUTS:  T - max thrust [N]                                       [scalar]
% OUTPUTS: s - numeric tag fragment, e.g. 10->'10', 2.5->'2p5'       [char]
%
% REFERENCES: [1] earth_elliptic_to_geo/mee_fuel_tag.m (the convention mirrored).
if abs(T - round(T)) < 1e-9
    s = sprintf('%d', round(T));
else
    s = strrep(sprintf('%g', T), '.', 'p');
end
end

% ---------------------------------------------------------------------------
function prevS = load_prev(T, prevT, reproDir)
% LOAD_PREV  Load the previous rung's engine-written row summary
% (results/repro/REPRO_row_T<10*prevT>.mat) for a 'chain' anchor or
% warm-started fuel stage. Errors loudly (not silently falling back to a
% cold solve) if that rung has not been reproduced yet -- a chained rung's
% correctness depends on the ACTUAL previous-rung solution, not a
% substitute.
%
% INPUTS:  T      - the rung currently being reproduced [N], for messages [scalar]
%          prevT   - the previous rung's thrust [N] to load               [scalar]
%          reproDir - results/repro directory                              [char]
%
% OUTPUTS: prevS - loaded struct with fields .row/.anchor/.sol/.rep (exactly
%          the variables reproduce_row.m saves for that rung)             [struct]
fname = fullfile(reproDir, sprintf('REPRO_row_T%d.mat', round(10*prevT)));
if ~isfile(fname)
    error('reproduce_row:prevMissing', ...
        ['T=%g N needs its predecessor rung''s reproduced result to warm-start ' ...
         'from -- run rung %g first (expected %s, not found)'], T, prevT, fname);
end
prevS = load(fname);
end

% ---------------------------------------------------------------------------
function rep = psr_report(psrOut)
% PSR_REPORT  Build a report struct (.m_f_kg/.switches/.revs/.edge/.incDeg/
% .defect/.certified) from psr_mee_refine.m's output, mirroring the shape of
% run_transfer_mee.m's res.report so downstream row assembly does not need
% to branch on whether PSR ran.
%
% INPUTS:  psrOut - psr_mee_refine.m output struct (.finalOut, .certified)  [struct]
% OUTPUTS: rep    - struct .m_f_kg .switches .revs .edge .incDeg .defect
%          .certified                                                     [struct]
fo = psrOut.finalOut;
rep = struct('m_f_kg', fo.m_f_kg, 'switches', fo.switches, 'revs', fo.dL/(2*pi), ...
    'edge', fo.edge, 'incDeg', fo.incDeg, 'defect', fo.maxDefect, ...
    'certified', psrOut.certified);
end

% ---------------------------------------------------------------------------
function tol = defaultTol(T)
% DEFAULTTOL  verify_row.m tolerance policy. Under the keep-best-mass /
% one-sided-verify decision (memory tenN-minfuel-razor-basin), the ONLY gate
% is mass: the reproduced solution must be at least as good as the campaign
% floor, i.e. m_f_kg >= cert.m_f_kg - tol.m_f_kg. tol.m_f_kg is a small slack
% (0.5 kg) absorbing numerical noise below the floor; the reproducer is
% expected to EQUAL or BEAT the campaign. Structure (switches/revs) is NOT
% gated (verify_row reports it), because the best min-fuel optimum can have a
% different bang-bang structure than the campaign row -- so no per-rung
% switch-tolerance policy is needed any more (the former +-3 PSR laxity is
% obsolete under one-sided-verify).
%
% INPUTS:  T   - max thrust [N]; one of 10, 5, 2.5, 1, 0.5              [scalar]
% OUTPUTS: tol - struct .m_f_kg (0.5 kg mass slack below the campaign floor) [struct]
if ~any(abs(T - [10, 5, 2.5, 1, 0.5]) < 1e-9)
    error('reproduce_row:noTolForThrust', ...
        'no default tolerance policy registered for thrustN=%g', T);
end
tol = struct('m_f_kg', 0.5);
end
