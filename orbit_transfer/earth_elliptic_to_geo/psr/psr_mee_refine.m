function out = psr_mee_refine(baseResult, opts)
% PSR_MEE_REFINE  PMP-Steered Refinement (PSR), ported to the L-domain MEE
% solver: switch-aware mesh refinement that inserts collocation nodes local
% to throttle switches instead of uniformly densifying the whole mesh.
%
% Each round: (1) locate throttle switches on the CURRENT sigma mesh and
% score each collocation interval by proximity to a switch (a local, tapered
% window, NOT the whole mesh); (2) bisect the top-scored intervals (original
% nodes preserved -- the "no-resample" discipline); (3) warm-interpolate the
% converged trajectory onto the new mesh (interp_warmstart.m: linear state +
% renormalized beta, nearest throttle -- keeps bang-bang edges crisp); (4)
% re-solve at eps=0, warmTight=true (feasibility-selected: only a solve that
% is BOTH success and defect/tfErr-certified is accepted and cached -- a
% looser stop, e.g. the legacy tulip-campaign 1e-3 stall stop, biases m_f,
% see REFERENCES [2]). Iterates until switch count and m_f stabilize between
% rounds, or a re-solve fails to certify, or maxRounds is hit.
%
% PORTED FROM vs ADAPTED (GTO_tulip/direct/PSR/lib/refine_loop.m +
% refine_sigma.m + pmp_refine_indicator.m -- read, NOT modified, NOT
% path-shared; this file is self-contained, mirroring PSR's own
% "vendored, not addpath'd" discipline for exactly the reason PSR's README
% documents: a later bare setup_paths() sharing a folder name can shadow
% masters):
%   - round loop / measure-then-refine-then-resolve structure: PORTED
%     directly (refine_loop.m's row-per-round history, "row 1 = seed, never
%     converged on its own" semantics, break-without-caching on a failed
%     resolve).
%   - mesh bisection (psr_refine_sigma_mee.m, factored into its own file for
%     no-solve testability, mirroring interp_warmstart.m's rationale):
%     PORTED directly from refine_sigma.m (preserve original nodes, insert
%     interval midpoints, hFloor guard, K/maxAdd caps, count-not-silently-
%     drop viable-but-capped selections).
%   - switch localization: ADAPTED, not ported. The tulip PSR scores
%     intervals from a PMP switching function S(tau) built off the direct
%     solve's own KKT defect DUALS (pmp_refine_indicator.m ->
%     sms_seed_duals's mode-'d' dual-to-costate map), letting it localize a
%     switch to a sub-interval crossing point. That machinery is a
%     substantial dependency (a dual->costate recovery pipeline) not yet
%     built for the MEE/L-domain formulation, and porting it is explicitly
%     OUT of this task's scope per the brief ("the primal throttle-crossing
%     version is acceptable for this port"). psr_switch_score_mee.m (also
%     factored into its own file, same testability rationale) instead
%     scores intervals from the PRIMAL throttle U(4,:) crossing 0.5
%     (exactly the crossings casadi_lt_mee.m's own out.switches already
%     counts via sum(abs(diff(thr>0.5)))), spread over a symmetric tapered
%     window of half-width opts.nbr intervals on each side (no sub-interval
%     S=0 root localization, since a bang-bang thr in {0,1} carries no
%     information between adjacent nodes to interpolate against).
%   - stabilization criterion: ADAPTED. refine_loop.m requires the switch
%     COUNT unchanged AND max switch MOVE < local mesh width AND
%     |dProp| < 1e-4 kg (a legacy tight campaign stop). This port instead
%     uses the brief's explicit criterion: |dm_f| < opts.propTol (default
%     0.1 kg) AND |d(switch count)| <= opts.swTol (default 2), over
%     opts.maxRounds (default 4) -- ~171 MEE switches vs the tulip campaign's
%     ~25-40 make an EXACT switch-count match too strict a bar to expect
%     within a 4-round budget, and a coarser |d(count)|<=2 is what the brief
%     asks for. No switch-move/localH check (no sub-interval switch position
%     is computed here -- see above).
%   - K/maxAdd DEFAULTS differ from refine_sigma.m's own (K=min(8,...)):
%     with ~171 switches (vs tulip's ~25-40), capping insertion to the
%     top-8 scored intervals per round would need >20 rounds to visit every
%     switch once, incompatible with a 4-round stabilization budget. Default
%     here is K = "all positively-scored intervals" (effectively unlimited),
%     capped only by a generous maxAdd safety valve (default 2000) -- see
%     psr_refine_sigma_mee.m, same mechanics as refine_sigma.m, different
%     numbers.
%
% INPUTS:
%   baseResult - the seed/current state to refine from, in run_transfer_mee's
%                'res' struct layout: .sigma [(N+1)x1], .fuel (struct with
%                .X [7x(N+1)] .U [4x(N+1)] .dL .success .ipoptStatus
%                .maxDefect .tfErr .m_f_kg -- i.e. a casadi_lt_mee.m out
%                struct), .cfg (.thrustN .m0kg .ispS .ctf .tfMinAnchor .tag),
%                .tf (fixedtf target, ND). baseResult.fuel need NOT already
%                be certified (the coarse-base node-economy probe feeds this
%                function an UNDER-resolved starting point on purpose -- PSR
%                is the mechanism that is supposed to fix that).
%   opts       - struct, all optional:
%                .maxRounds [4], .tag [baseResult.cfg.tag + '_PSR'],
%                .outDir [earth_elliptic_to_geo/results], .maxIter [1500] a
%                FLOOR, not a fixed cap (Task 9 Step 0 fix): the EFFECTIVE
%                per-round cap is psr_iter_cap(opts.maxIter, N_round) =
%                max(opts.maxIter, ceil(N_round)), so later/bigger rounds get
%                proportionally more iterations instead of the same fixed
%                budget every round; if a round's first attempt at that cap
%                still fails to certify, ONE retry is made at 2x the cap
%                before the round is treated as resolveFailed (see
%                REFERENCES [9]),
%                .nbr [2] switch-window half-width in COLLOCATION INTERVALS,
%                .K [Inf -> all positively-scored intervals], .hFloor [1e-9]
%                min sigma-interval width to bisect, .maxAdd [2000] safety
%                cap on inserted nodes per round, .propTol [0.1] kg,
%                .swTol [2] switches, .globalEvery [0 -> disabled] period (in
%                MEASURED rounds) for a periodic GLOBAL densify round (Task 9
%                Step 0 hybrid-refinement fix, see psr_is_global_round.m):
%                every opts.globalEvery-th measured round, the refinement
%                step uses psr_global_score_mee.m (evenly-spaced, switch-
%                blind densify by opts.globalFactor) INSTEAD OF the default
%                switch-local psr_switch_score_mee.m score, so a switch pair
%                hiding in a region with no visible throttle crossing yet
%                gets a chance to surface on the next resolve -- windowed
%                PSR scoring alone is structurally blind to switches it
%                hasn't already seen (Task 8's uniform-1.5x-vs-PSR
%                comparison: sw=177 vs PSR's flat 171 at the same node
%                budget). Switch COUNTING/history (nSw, swSigma) always uses
%                the switch scorer regardless of which round type is
%                refining, .globalFactor [1.3] densify factor for a global
%                round (see psr_global_score_mee.m).
%
% OUTPUTS:
%   out - struct: .history [1xR] struct array (round [0-based], nNodes,
%         switches, m_f_kg, maxDefect, tfErr, certified, converged, success,
%         ipoptStatus, swSigma [switch locations in sigma]), .finalSigma
%         [(N'+1)x1], .finalOut (casadi_lt_mee out struct for the final
%         measured round -- on a 'resolveFailed' stop this is the LAST
%         CERTIFIED round, NOT the failed attempt: the loop never advances
%         curSigma/curOut past an uncertified resolve, so PSR always hands
%         back its best certified state rather than a broken one),
%         .certified (logical, final round only),
%         .stabilized (logical, true iff the loop stopped because the
%         stabilization criterion was met, not maxRounds/failure/no-refine),
%         .stopReason ('stabilized'|'maxRounds'|'noRefinable'|'resolveFailed'),
%         .tag, .baseTag, .opts (resolved, defaults filled in). Saved to
%         resDir/[tag '_psr_final.mat'] ONLY if out.certified (never cache
%         uncertified -- house rule). Per-round certified re-solves cached to
%         resDir/[tag '_psr_r<k>.mat'] (config-fingerprinted, resume-safe:
%         an uncertified round is NEVER cached, so a resumed run always
%         re-attempts it rather than silently trusting a partial result).
%         The round-by-round history is unconditionally persisted every
%         round to resDir/[tag '_psr_history.mat'] (a diagnostic log, not a
%         warm-startable solution -- same distinction refine_loop.m draws
%         between its history file and its solFile).
%
% REFERENCES:
%   [1] GTO_tulip/direct/PSR/lib/refine_loop.m (round-loop structure,
%       read-only source of this port).
%   [2] GTO_tulip/direct/PSR/lib/refine_sigma.m (bisection mechanics,
%       read-only source of this port).
%   [3] GTO_tulip/direct/PSR/lib/pmp_refine_indicator.m (the
%       dual-based switching-function scorer this file's primal-throttle
%       scorer stands in for; NOT ported -- see ADAPTED note above).
%   [4] earth_elliptic_to_geo/psr_switch_score_mee.m (the ADAPTED primal
%       scorer this file calls, factored out for no-solve testability).
%   [5] earth_elliptic_to_geo/psr_refine_sigma_mee.m (the PORTED bisection
%       helper this file calls, factored out for no-solve testability).
%   [6] earth_elliptic_to_geo/interp_warmstart.m (mesh-refine warm-start,
%       Task 7 machinery, reused as-is).
%   [7] earth_elliptic_to_geo/casadi_lt_mee.m (the L-domain solver every
%       round re-solves against).
%   [8] .superpowers/sdd/task-8-brief.md (this task's spec).
%   [9] .superpowers/sdd/task-9-brief.md Step 0 (the N-scaled maxIter +
%       2x-retry fix, and the periodic global-round hybrid fix, both added
%       here). psr_iter_cap.m / psr_is_global_round.m / psr_global_score_mee.m
%       (the three pure-function decision helpers those fixes are factored
%       into for no-solve testability).

resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
if nargin < 2, opts = struct(); end
d = @(f, v) optdef(opts, f, v);

cfg = baseResult.cfg;
baseTagStr = optdef(cfg, 'tag', 'MEE_M2_unknown');

opts.maxRounds = d('maxRounds', 4);
opts.tag       = d('tag', [baseTagStr '_PSR']);
opts.outDir    = d('outDir', resDir);
opts.maxIter   = d('maxIter', 1500);
opts.nbr       = d('nbr', 2);
opts.K         = d('K', Inf);
opts.hFloor    = d('hFloor', 1e-9);
opts.maxAdd    = d('maxAdd', 2000);
opts.propTol   = d('propTol', 0.1);
opts.swTol     = d('swTol', 2);
opts.globalEvery  = d('globalEvery', 0);
opts.globalFactor = d('globalFactor', 1.3);
tag    = opts.tag;
outDir = opts.outDir;
if ~exist(outDir, 'dir'), mkdir(outDir); end

assert(isfield(cfg, 'thrustN') && ~isempty(cfg.thrustN), 'psr_mee_refine:noThrustN', ...
    ['baseResult.cfg has no .thrustN -- some early run_transfer_mee.m callers saved the ' ...
     'ORIGINAL (possibly empty) cfg argument verbatim into res.cfg rather than the resolved ' ...
     'config (e.g. MEE_M2_10N.mat''s res.cfg is a no-field struct); pass a baseResult whose ' ...
     'cfg.thrustN/.m0kg/.ispS are populated (e.g. MEE_M2_1N.mat / MEE_M2_5N.mat / ' ...
     'MEE_M2_2p5N.mat, or a hand-built struct for the node-economy probe)']);
m0kg = optdef(cfg, 'm0kg', 1500);
ispS = optdef(cfg, 'ispS', 2000);
par  = kepler_lt_params(cfg.thrustN, m0kg, ispS);
tf   = optdef(baseResult, 'tf', optdef(cfg, 'ctf', 1.5) * optdef(cfg, 'tfMinAnchor', NaN));
assert(~isnan(tf), 'psr_mee_refine:noTf', ...
    'baseResult carries neither .tf nor cfg.ctf/cfg.tfMinAnchor -- cannot fix the fixedtf target');

curSigma = baseResult.sigma(:);
curOut   = baseResult.fuel;
x0       = curOut.X(:, 1);

historyFile = fullfile(outDir, [tag '_psr_history.mat']);
history = struct([]);
prevSwitch = NaN;  prevMf = NaN;
stopReason = 'maxRounds';   % default if the loop runs out without another cause

for row = 1:(opts.maxRounds + 1)
    [swIdx, score] = psr_switch_score_mee(curSigma, curOut.U, opts);
    nSw = numel(swIdx);
    mfk = curOut.m_f_kg;
    tfE = optdef(curOut, 'tfErr', NaN);
    certified = is_psr_certified(curOut);

    converged = false;
    if row > 1
        dmf = abs(mfk - prevMf);
        dsw = abs(nSw - prevSwitch);
        converged = dsw <= opts.swTol && dmf < opts.propTol;
    end

    history(row).round     = row - 1; %#ok<AGROW>
    history(row).nNodes     = numel(curSigma) - 1;
    history(row).switches   = nSw;
    history(row).swSigma    = curSigma(swIdx).';
    history(row).m_f_kg     = mfk;
    history(row).maxDefect  = curOut.maxDefect;
    history(row).tfErr      = tfE;
    history(row).certified  = certified;
    history(row).converged  = converged;
    history(row).success    = optdef(curOut, 'success', false);
    history(row).ipoptStatus = optdef(curOut, 'ipoptStatus', '');
    save(historyFile, 'history');
    fprintf(['[psr round %2d] N=%5d sw=%4d mf=%9.4f kg defect=%.2e tfErr=%.2e ' ...
             'cert=%d conv=%d\n'], row - 1, numel(curSigma) - 1, nSw, mfk, ...
            curOut.maxDefect, tfE, certified, converged);

    if converged
        stopReason = 'stabilized';  break;
    end
    if row > opts.maxRounds
        stopReason = 'maxRounds';  break;
    end
    prevSwitch = nSw;  prevMf = mfk;

    measuredRound = row - 1;    % 0-based, matches history(row).round
    globalRound = psr_is_global_round(measuredRound, opts.globalEvery);
    if globalRound
        refScore = psr_global_score_mee(numel(curSigma) - 1, opts.globalFactor);
        fprintf('  [round %d] GLOBAL densify round (every %d, factor %.2f)\n', ...
                measuredRound, opts.globalEvery, opts.globalFactor);
    else
        refScore = score;
    end

    [sigmaNew, isNew, nDropped] = psr_refine_sigma_mee(curSigma, refScore, opts);
    if nDropped > 0
        fprintf('  psr_refine_sigma_mee dropped %d viable interval(s) (hFloor/maxAdd)\n', nDropped);
    end
    if nnz(isNew) == 0
        fprintf('  no refinable intervals (all sub-hFloor or no switches); stopping.\n');
        stopReason = 'noRefinable';  break;
    end

    W = interp_warmstart(curOut.X, curOut.U, curOut.dL, curSigma, sigmaNew);
    roundIdx  = row;   % this resolve produces the state measured at row+1
    roundFile = fullfile(outDir, sprintf('%s_psr_r%d.mat', tag, roundIdx));
    fpRound = struct('baseTag', baseTagStr, 'thrustN', par.thrustN, 'm0kg', m0kg, ...
        'ispS', ispS, 'tf', tf, 'nbr', opts.nbr, 'K', opts.K, 'hFloor', opts.hFloor, ...
        'maxAdd', opts.maxAdd, 'maxIter', opts.maxIter, 'round', roundIdx, ...
        'Nsrc', numel(curSigma) - 1, 'Ndst', numel(sigmaNew) - 1);

    if isfile(roundFile)
        S = load(roundFile);
        o = S.o;
        check_cache_fp_psr(S, fpRound, roundFile, tag);
        fprintf('  [round %d] loaded certified cache %s\n', roundIdx, roundFile);
        okCert = true;   % never cached unless certified at save time
    else
        [o, iterCapUsed, retried] = solve_psr_round(sigmaNew, W, par, tf, x0, opts.maxIter);
        okCert = is_psr_certified(o);
        if okCert
            save(roundFile, 'o', 'sigmaNew', 'isNew', 'fpRound');
        else
            retryMsg = '';
            if retried, retryMsg = sprintf(' (after retry at 2x cap=%d)', iterCapUsed); end
            fprintf(['  [round %d] resolve did NOT certify%s (iterCap=%d, defect=%.2e ' ...
                     'tfErr=%.2e status=%s); NOT cached, stopping refinement.\n'], ...
                    roundIdx, retryMsg, iterCapUsed, o.maxDefect, o.tfErr, o.ipoptStatus);
        end
    end

    if ~okCert
        stopReason = 'resolveFailed';  break;
    end
    curSigma = sigmaNew;  curOut = o;
end

out = struct('history', history, 'finalSigma', curSigma, 'finalOut', curOut, ...
    'certified', history(end).certified, 'stabilized', strcmp(stopReason, 'stabilized'), ...
    'stopReason', stopReason, 'tag', tag, 'baseTag', baseTagStr, 'opts', opts);

if out.certified
    finalFile = fullfile(outDir, [tag '_psr_final.mat']);
    fpFinal = struct('baseTag', baseTagStr, 'thrustN', par.thrustN, 'm0kg', m0kg, ...
        'ispS', ispS, 'tf', tf, 'nbr', opts.nbr, 'K', opts.K, 'hFloor', opts.hFloor, ...
        'maxAdd', opts.maxAdd, 'maxIter', opts.maxIter, 'maxRounds', opts.maxRounds, ...
        'propTol', opts.propTol, 'swTol', opts.swTol, 'nRoundsRun', row - 1, ...
        'Nfinal', numel(curSigma) - 1);
    save(finalFile, 'out', 'fpFinal'); %#ok<*NASGU>
    fprintf('psr_mee_refine: wrote certified final result %s\n', finalFile);
else
    warning('psr_mee_refine:uncertified', ['%s: final measured state is NOT certified ' ...
        '(defect=%.2e) -- not cached as final (house rule: never cache uncertified)'], ...
        tag, curOut.maxDefect);
end
end

% ---------------------------------------------------------------------------
function ok = is_psr_certified(o)
% IS_PSR_CERTIFIED  Single feasibility-selected certification predicate --
% success && maxDefect<1e-8 && (tfErr is NaN or tfErr<1e-8). Factored out of
% 3 separately-inlined copies of the same test (Task 9 Step 0 review finding:
% the round-0/seed measurement check, the post-resolve accept-and-cache
% check, and the loop's own resolveFailed check had each spelled out the
% same predicate independently, a drift risk if any one of them were ever
% edited alone). Semantics UNCHANGED from all three original inlined copies.
tfE = NaN;
if isfield(o, 'tfErr'), tfE = o.tfErr; end
ok = isfield(o, 'success') && o.success && isfield(o, 'maxDefect') ...
     && o.maxDefect < 1e-8 && (isnan(tfE) || tfE < 1e-8);
end

% ---------------------------------------------------------------------------
function [o, iterCapUsed, retried] = solve_psr_round(sigmaNew, W, par, tf, x0, floorIter)
% SOLVE_PSR_ROUND  One fresh casadi_lt_mee resolve at an N-scaled iteration
% cap (psr_iter_cap.m), with a SINGLE retry at 2x that cap if the first
% attempt fails to certify (Task 9 Step 0 fix: Task 8's validation runs both
% stalled via 'resolveFailed' at N~3000-3365 under the OLD flat
% opts.maxIter=1500 default, at defects only 1-2 orders of magnitude above
% the certify gate -- 5.7e-7/7.95e-7 vs <1e-8 -- suggesting the solves were
% close, not diverging). The retry reuses the SAME warm start W (the
% interpolated trajectory, not the maxed-out iterate) -- an IPOPT iterate
% that stopped short of the barrier-clearing point at maxIter is not
% obviously a better warm start than the original interpolation, and reusing
% it would complicate the never-cache-uncertified bookkeeping.
%
% INPUTS:  sigmaNew  - refined node grid for this round [(N'+1)x1]
%          W         - interp_warmstart.m output (.X .U .dL) for sigmaNew
%          par       - kepler_lt_params output for this rung's thrust
%          tf        - fixedtf target [ND]
%          x0        - initial state [7x1]
%          floorIter - opts.maxIter, the user-set FLOOR [scalar]
%
% OUTPUTS: o           - casadi_lt_mee out struct (the retry's result if a
%                        retry happened and improved on the first attempt,
%                        else the first attempt)
%          iterCapUsed - the iteration cap actually used to produce o
%                        [scalar]
%          retried     - true iff a second (2x-cap) attempt was made
%                        [logical]
%
% REFERENCES: [1] psr_iter_cap.m (the N-scaling formula). [2] task-9-brief.md
%   Step 0 item 1. [3] psr_mee_refine.m (sole caller).
Nround  = numel(sigmaNew) - 1;
iterCap = psr_iter_cap(floorIter, Nround);
o = casadi_lt_mee(sigmaNew, W.X, W.U, W.dL, struct('par', par, 'mode', 'fixedtf', ...
    'eps', 0, 'tfTarget', tf, 'x0', x0, 'maxIter', iterCap, ...
    'warmTight', true, 'printLevel', 0));
retried = false;
if ~is_psr_certified(o)
    iterCap2 = 2 * iterCap;
    fprintf(['  [psr retry] round did not certify at maxIter=%d (defect=%.2e, status=%s); ' ...
             'retrying once at maxIter=%d before giving up\n'], iterCap, o.maxDefect, o.ipoptStatus, iterCap2);
    o2 = casadi_lt_mee(sigmaNew, W.X, W.U, W.dL, struct('par', par, 'mode', 'fixedtf', ...
        'eps', 0, 'tfTarget', tf, 'x0', x0, 'maxIter', iterCap2, ...
        'warmTight', true, 'printLevel', 0));
    retried = true;
    if is_psr_certified(o2) || o2.maxDefect < o.maxDefect
        o = o2;  iterCap = iterCap2;
    end
end
iterCapUsed = iterCap;
end

% ---------------------------------------------------------------------------
function check_cache_fp_psr(S, fp, file, tag)
% CHECK_CACHE_FP_PSR  Fail-loud cache-fingerprint guard (mirrors
% run_transfer_mee.m's check_cache_fp): error out, naming the first
% mismatched field, if a loaded round cache's stored fingerprint disagrees
% with the current config under the same tag/round index.
if ~isfield(S, 'fpRound')
    warning('psr_mee_refine:noCachedFingerprint', ['%s has no cached config ' ...
        'fingerprint -- trusting it because tag=''%s'' matches'], file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fpRound, f) || ~isequal(S.fpRound.(f), fp.(f))
        error('psr_mee_refine:fingerprintMismatch', ['cached config fingerprint ' ...
            'mismatch in %s: field ''%s'' differs between the cache and the ' ...
            'current config -- stale cache from a different configuration under ' ...
            'tag=''%s''; delete the file or use a new opts.tag'], file, f, tag);
    end
end
end
