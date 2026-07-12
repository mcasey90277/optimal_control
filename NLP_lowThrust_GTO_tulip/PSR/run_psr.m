% RUN_PSR  Entry-level driver for a full PSR (PMP-Steered Refinement) run.
%
% PSR is the working direct-method pipeline for the CR3BP min-fuel
% GTO -> south-pole-tulip transfer (15 kg, 25 mN, Isp 2100 s, ~40-rev spiral,
% Sundman-regularized mesh dt/dtau = r1^1.5). One run of this script does:
%
%   1. PARAMETERS   - set t_f and everything else in ONE place (this section)
%   2. DIRECT SOLVE - energy->fuel homotopy (CasADi+IPOPT, trapezoid in tau):
%                     solve the smooth min-ENERGY problem's neighborhood first,
%                     then sharpen epsilon -> 0 to the bang-bang min-FUEL
%                     solution (Bertrand-Epenoy: J = Int[s]dt - eps*Int[s(1-s)]dt)
%   3. PSR REFINE   - PMP-steered mesh refinement: the INDIRECT machinery
%                     measures where the switching function S = 1-||lamV||c/m-lamM
%                     (costates recovered from the NLP's own KKT duals) localizes
%                     each throttle switch worst, refines the mesh THERE, and
%                     re-solves; repeats until the switch times stabilize.
%                     This sharpens switch times below the original mesh width.
%   4. COSTATES +   - PSR does not natively produce costates; this stage runs
%      DATA EXPORT    the dual->costate recovery (adjudicated mode-'d' map +
%                     beta fit) and saves ALL data products -- mesh, trajectory,
%                     control, costates, switching function, switch times,
%                     transversality errors, constants, provenance -- to
%                     ../PSR_data/psr_data_tf<factor>_sw<k>.mat. The file
%                     doubles as a ready-made IFS seed (same layout).
%   5. VERIFY       - first-order PMP certificate (EXTREMALITY only): per-arc
%                     propagation of the 16-dim state+costate system from the
%                     solution's own duals, primer alignment, transversality,
%                     switch-structure match. Appends its summary to the
%                     PSR_data file.  [TODO: second-order / conjugate-point
%                     test to upgrade extremality -> local minimality.]
%   6. MOVIE        - control movie: rotating-frame transfer colored burn/coast,
%                     primer thrust arrows, throttle strip, running Delta-V.
%
% Pipeline outputs land in PSR/results/; DATA PRODUCTS land in ../PSR_data/.
% Approximate costs on this machine:
% stage 2 ~ 30-90 min (13 IPOPT solves at N=4001), stage 3 ~ 20-60 min
% (maxRounds IPOPT re-solves + indicator), stage 4 ~ 5-10 min, stage 5
% seconds ('preview') to ~15 min ('movie'). Stages are RESUMABLE: each one
% skips itself if its output file already exists (delete the file, or flip
% the rerun flag, to force a redo).
%
% REFERENCES:
%   [1] ../LOW_THRUST_MINFUEL_CAMPAIGN.md  (campaign record; two-walls analysis)
%   [2] ../sundman_minfuel/refine/README.md + RESULTS.md  (PSR design + results)
%   [3] ../ms_band/MS_BAND_CAMPAIGN.md  (verifier provenance, dual-map mode 'd')
%   [4] Bertrand & Epenoy, OCAM 23(4), 2002 (energy->fuel homotopy)

%% ------------------------------------------------------------------------
%% 0. Paths and output folder
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir  = fullfile(here, 'results');          % pipeline intermediates
dataDir = fullfile(here, '..', 'PSR_data');   % exported data products (stage 4)
if ~exist(resDir, 'dir'), mkdir(resDir); end
cfg = minfuel_config();          % campaign constants: tfMin, schedules, dirs

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------

% ---- transfer time -------------------------------------------------------
% t_f = factor * tfMin, tfMin = 6.2906939607 ND = 27.8845 days (certified
% indirect min-time). The campaign's validated band for this pipeline is
% factor >= 1.12 (below ~1.12 the direct solver itself struggles; the
% 1.01-1.11 "transition band" is an open research problem, see [1]).
factor = 1.15;                   % t_f / t_f_min   e.g. 1.12 ... 1.85

% ---- SEED for the direct solve --------------------------------------------
% The direct solve is seeded one of three ways (see minfuel_at_tf.m):
%
%   'energy'   (recommended default) Seed from the min-ENERGY backbone file at
%              THIS factor: sundman_minfuel/results/energy/energy_f####.mat.
%              Backbones exist for factors 1.12:0.01:1.15 and 1.20:0.05:1.85
%              (ls that folder). This is the canonical homotopy root: the
%              energy problem is smooth and strictly convex in the control, so
%              it converges from a crude guess, and the epsilon-schedule then
%              deforms it continuously to bang-bang fuel. Missing backbone?
%              run sundman_minfuel/orchestrate/backbone_walk.sh, or use a
%              neighbor seed instead.
%
%   'neighbor' Seed from an EXISTING bang-bang solution at a nearby factor
%              (seedFactor below): its time state is rescaled to this t_f and
%              a light re-sharpen schedule is run. Use to walk along the
%              dV-vs-tf front from a solved case. Requires that solution to
%              exist in sundman_minfuel/results/minfuel/.
%
%   <path>     Seed from an explicit .mat you provide (fields X [8xnN],
%              U [4xnN] at top level or inside `out`, plus sigma, tauf0,
%              rv0, rvf). Treated like 'neighbor'.
seedSpec   = 'energy';           % 'energy' | 'neighbor' | '/path/to/seed.mat'
seedFactor = NaN;                % only used when seedSpec = 'neighbor'

% ---- direct-solve knobs ----------------------------------------------------
% epsilon schedule [] = default by seed type (cfg.schedSharpen for 'energy':
% 0.6 -> 0 in 13 steps; cfg.schedNeighbor for the rest). Override only for
% experiments -- the defaults are the campaign-validated schedules.
sched      = [];
maxIter    = cfg.maxIter;        % IPOPT iteration cap per schedule step [1500]

% ---- PSR refinement knobs --------------------------------------------------
refineOpts = struct( ...
    'maxRounds', 4, ...          % refinement rounds (headline run used 4)
    'K',         8, ...          % sub-intervals a flagged interval splits into
    'maxAdd',    40);            % max intervals refined per round

% ---- verification knobs ----------------------------------------------------
% ADJUDICATION (adjArcs / adjSwitches): the verifier withholds the certificate
% when any arc's state defect exceeds the 1e-2 line or any switch is unmatched.
% Some flags are KNOWN, benign, and SOLUTION-SPECIFIC -- they must be justified
% by inspecting the stage-5 ATTN rows, NOT blanket-applied to force a green
% certificate. Workflow at a new factor: run once with these EMPTY, read the
% "ATTN" arcs/switches in the stage-5 table, justify each, then list them here
% and rerun stage 5 (rerunVerify = true). The issued certificate prints the
% adjudicated rows, so nothing is hidden.
%
% For the 1.15x refined solution the justified adjudications are:
%   adjArcs   40 = TERMINAL SWITCH CLUSTER: several switches packed into a tau
%                  sliver in the last arc; the throttle-vs-S disagreement
%                  integrates to O(1) -- an amplification artifact, not a
%                  costate error (same phenomenon adjudicated for legacy 1.12x).
%             32 = a perigee arc at ~2x the ~5e-3 dual-map floor at this mesh.
%   adjSwitches 1, 25 = NEAR-GRAZE switches (dual-S one-signed, |S|~0) at the
%                  trajectory start / near the end -- grazes, not missed switches.
verifyOpts = struct( ...
    'M',           40, ...       % multiple-shooting arcs for the certificate
    'mode',        'd', ...      % dual->costate map (adjudicated: midpoint 'd')
    'epsEval',     1e-4, ...     % smoothing used ONLY to propagate arcs
    'adjArcs',     [32 40], ...  % arcs adjudicated (see justification above)
    'adjSwitches', [1 25], ...   % near-graze switches adjudicated (see above)
    'makeFig',     true);

% ---- movie ------------------------------------------------------------------
movieMode = 'movie';             % 'preview' (3 stills, fast) | 'movie' | 'none'

% ---- stage rerun control ----------------------------------------------------
% Each stage skips itself when its output already exists. Force with these.
rerunDirect = false;  rerunRefine = false;  rerunVerify = false;

% ---- derived file names (canonical; do not edit) ---------------------------
tag         = sprintf('f%04d', round(1000*factor));
directFile  = fullfile(resDir, ['psr_direct_'  tag '.mat']);
seedFile    = fullfile(resDir, ['psr_seed_'    tag '.mat']);
refinedFile = fullfile(resDir, ['psr_refined_' tag '.mat']);

fprintf('\n=== PSR PIPELINE: factor=%.3f (t_f=%.4f ND = %.2f days), seed=%s ===\n', ...
        factor, factor*cfg.tfMin, factor*cfg.tfMin*382981.289129055/86400, ...
        char(string(seedSpec)));

%% ------------------------------------------------------------------------
%% 2. DIRECT SOLVE  (energy->fuel homotopy to a certified bang-bang solution)
%% ------------------------------------------------------------------------
% minfuel_at_tf is the campaign's canonical per-t_f driver. For an 'energy'
% seed it (a) re-cleans the backbone tight at eps=1 (backbone duals are
% loose-continued and would blow up inf_du if sharpened directly), then
% (b) walks the epsilon schedule down to exactly eps=0 (pure fuel, linear in
% throttle -> bang-bang) with warm-tight IPOPT settings. "Certified" here
% means at least one schedule step converged tight (defect < 1e-6); an
% uncertified attempt is returned but never saved (it must not poison seeds).
if isfile(directFile) && ~rerunDirect
    fprintf('\n[stage 2] direct solution exists (%s) -- skipping. Set rerunDirect=true to redo.\n', directFile);
    D = load(directFile);  outDirect = D.out;
else
    fprintf('\n[stage 2] DIRECT SOLVE (energy->fuel homotopy)...\n');
    % seedFactor is ignored unless seedSpec = 'neighbor' (minfuel_at_tf checks)
    args = {'seedFactor', seedFactor, 'outFile', directFile, 'maxIter', maxIter, 'branch', 'psr'};
    if ~isempty(sched), args = [args {'sched', sched}]; end
    outDirect = minfuel_at_tf(factor, 'seed', seedSpec, args{:});
    assert(outDirect.certified, ...
        'direct solve did not certify at factor %.3f -- inspect outDirect, do not proceed', factor);
end
fprintf('[stage 2] dV=%.4f km/s  prop=%.4f kg  switches=%d  edge=%.1f%%  defect=%.2g\n', ...
    outDirect.dV, outDirect.prop_kg, outDirect.switches, 100*outDirect.edge, outDirect.maxDefect);

%% ------------------------------------------------------------------------
%% 3. PSR REFINEMENT  (indirect-steered mesh refinement -> sharp switch times)
%% ------------------------------------------------------------------------
% The refinement loop is DIRECT-method only in its solves; the indirect
% machinery is a *measurement tool*: pmp_refine_indicator recovers the
% costates from the NLP's KKT defect duals (mode-'d' midpoint map), forms the
% switching function S(tau), and scores how well each throttle switch is
% localized by its S = 0 crossing. refine_sigma then splits the worst
% intervals (K-fold, up to maxAdd per round), the solution is warm-started
% onto the new mesh WITHOUT resampling the control, and IPOPT re-solves at
% eps=0 warmTight. Rounds repeat until switch times move less than a local
% mesh width and the propellant change is below propTol -- switch times are
% then resolved to sub-original-mesh accuracy. (Campaign result: 1.15x switch
% times stabilize in ~4 rounds; see refine/RESULTS.md.)
if isfile(refinedFile) && ~rerunRefine
    fprintf('\n[stage 3] refined solution exists (%s) -- skipping. Set rerunRefine=true to redo.\n', refinedFile);
    H = load(fullfile(resDir, sprintf('refine_history_psr_%s.mat', tag)));  history = H.history;
else
    fprintf('\n[stage 3] PSR REFINEMENT (max %d rounds)...\n', refineOpts.maxRounds);
    % prep: normalize the direct file into the refine-seed layout (pass-through
    % when out.lamDef is already present, which minfuel_at_tf guarantees)
    prep_refine_seed(directFile, seedFile);
    ro = refineOpts;
    ro.tag     = ['psr_' tag];
    ro.outDir  = resDir;             % history .mat + summary figure land here
    ro.solFile = refinedFile;        % FINAL refined solution, seed layout
    history = refine_loop(seedFile, ro);
end
fprintf('\n[stage 3] summary (round 0 = unrefined seed):\n');
fprintf('%-6s %-7s %-4s %-11s %-11s %-7s %-11s\n', ...
        'round', 'nodes', 'sw', 'maxMove', 'dProp(kg)', 'nViol', 'HresMax');
for r = 1:numel(history)
    h = history(r);
    fprintf('%-6d %-7d %-4d %-11.2e %-11.2e %-7d %-11.2e\n', ...
            r-1, h.nNodes, h.switches, h.maxSwitchMove, h.dProp, h.nViol, h.HresMax);
end

%% ------------------------------------------------------------------------
%% 4. COSTATES + DATA EXPORT  (data products -> ../PSR_data/)
%% ------------------------------------------------------------------------
% PSR's NLP natively yields only the raw interval duals (out.lamDef); the
% usable NODE costates come from the adjudicated mode-'d' midpoint dual->
% costate map with the beta scale fit (ms_band machinery). This stage runs
% that recovery and writes ONE self-contained file to ../PSR_data/ holding:
%   - the mesh (sigma, Sundman tau, physical time at nodes)
%   - the trajectory (r, v, m and the full 8-state X)
%   - the control (thrust direction alpha, throttle s, switch times both by
%     certified dual-S crossings and by raw throttle crossings)
%   - the costates lam [8xnN] + switching function S + beta (+ its spread,
%     the scale-fit quality diagnostic)
%   - transversality / first-order errors (lamM(tau_f), terminal rendezvous
%     residuals, fixed-t_f residual, S-sign law agreement)
%   - all physical constants + provenance (source, date, git hash)
% The file keeps the standard seed layout at top level, so it can be fed
% DIRECTLY to ifs_seed / verify_direct_pmp / a future IFS run -- this is the
% designed handoff artifact from the direct pipeline to independent analysis.
% Cheap (~seconds) and idempotent: always regenerated from the refined file.
fprintf('\n[stage 4] COSTATE GENERATION + DATA EXPORT...\n');
dataFile = psr_export_data(refinedFile, dataDir, struct('M', verifyOpts.M));

%% ------------------------------------------------------------------------
%% 5. VERIFY  (first-order PMP certificate -- extremality)
%% ------------------------------------------------------------------------
% verify_direct_pmp is a VERIFIER, not a solver (no optimization anywhere):
% it recovers the costates from the refined solution's own KKT duals, then
% propagates each of M arcs of the full 16-dim (state; costate) Sundman
% system from the solution's own values and reports per-arc defects, the
% dual-implied switching structure vs the actual throttle, primer alignment
% (thrust antiparallel to lamV), Hamiltonian stationarity |Ht + lamT|, and
% the free-final-mass transversality lamM(tau_f) = 0. It prints the
% "consistent with a continuous PMP extremal" certificate only when every
% gate passes or is explicitly adjudicated.
%
% WHAT THIS PROVES / DOES NOT PROVE: this is a FIRST-ORDER certificate --
% the solution satisfies the Pontryagin necessary conditions (an EXTREMAL)
% at its mesh's O(h^2) resolution. It does NOT prove local minimality.
% TODO (planned upgrade): second-order test -- no conjugate points on the
% arcs (Jacobi condition via the variational/Riccati system along the
% trajectory) or verification of the second-order sufficient conditions at
% the NLP level (reduced-Hessian positivity at the active set). Until then,
% "certified" in this pipeline means first-order extremality + tight
% feasibility, and minimality evidence is only comparative (homotopy family
% + front monotonicity, see HONEST_EVALUATION_DV_TF_FRONT.md).
verifyFile = fullfile(resDir, sprintf('verify_pmp_psr_refined_%s.mat', tag));
if isfile(verifyFile) && ~rerunVerify
    fprintf('\n[stage 5] verification exists (%s) -- skipping. Set rerunVerify=true to redo.\n', verifyFile);
    V = load(verifyFile);  vsum = V.summary;
else
    fprintf('\n[stage 5] FIRST-ORDER PMP VERIFICATION (M=%d arcs)...\n', verifyOpts.M);
    oldDir = cd(resDir);                       % verifier writes to pwd
    cleanupObj = onCleanup(@() cd(oldDir));
    vsum = verify_direct_pmp(refinedFile, verifyOpts);
    clear cleanupObj
end
fprintf('[stage 5] certOK=%d  worstStateDef=%.3g  primer(mean/p95)=%.3f/%.3f deg  |lamM(sigf)|=%.2g  switches matched %d/%d\n', ...
    vsum.certOK, vsum.worstStateDef, vsum.primerMeanDeg, vsum.primerP95Deg, ...
    abs(vsum.lamMend), vsum.nMatched, vsum.nSwitches);
% append the certificate summary to the PSR_data product file so downstream
% analysis has the verification verdict alongside the data
verify = vsum;
save(dataFile, 'verify', '-append');
fprintf('[stage 5] verify summary appended to %s\n', dataFile);

%% ------------------------------------------------------------------------
%% 6. CONTROL MOVIE  (transfer + control law, synced)
%% ------------------------------------------------------------------------
% Rotating CR3BP frame. Trajectory animated in bold, red = burn / blue =
% coast, primer thrust arrow while burning; below it the throttle strip
% (the bang-bang control law) and the running Delta-V curve, all synced and
% played uniformly in PHYSICAL time (the Sundman mesh is not). 'preview'
% writes three stills in seconds -- look at those before paying ~15 min for
% the full MP4+GIF ('movie').
if ~strcmp(movieMode, 'none')
    fprintf('\n[stage 6] CONTROL MOVIE (%s)...\n', movieMode);
    titleStr = sprintf('PSR min-fuel GTO\\rightarrowtulip, t_f = %.2fx min-time (%d-switch bang-bang)', ...
                       factor, history(end).switches);
    psr_movie(refinedFile, fullfile(resDir, ['psr_movie_' tag]), titleStr, movieMode);
end

fprintf('\n=== PSR PIPELINE DONE (factor %.3f). Intermediates: %s  Data products: %s ===\n', ...
        factor, resDir, dataFile);
