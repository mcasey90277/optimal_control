% RUN_PSR_BATCH  Sweep the PSR pipeline over many t_f factors.
%
% The batch sibling of run_psr.m: instead of one factor, run the full PSR
% pipeline for a VECTOR of t_f factors (or automatically for every factor that
% has a min-energy seed on disk), at a single homotopy endpoint epsMin. Each
% factor runs the same stages as run_psr -- direct energy->fuel solve, PSR
% refinement (bang-bang only), costate + data export, and the IPOPT-native
% local-minimality certificate -- with a per-factor try/catch so one failure
% does not kill the sweep, and a summary table at the end.
%
% Movies are OFF by default (a batch of movies is hours); enable per taste.
% Every factor is RESUMABLE: stages skip when their output already exists.
%
% Edit section 1, then run. Outputs: PSR/results/ (intermediates), ../PSR_data/
% (data products), and a printed + saved summary table.
%
% REFERENCES: run_psr.m (the single-factor driver whose stages this mirrors);
%   gen_energy_seed.m (make an energy seed for a factor you don't have yet;
%   usable band ~[1.12x, 1.95x] -- see LOW_THRUST_MINFUEL_CAMPAIGN.md).

%% ------------------------------------------------------------------------
%% 0. Paths
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
P.resDir  = fullfile(here, 'results');
P.dataDir = fullfile(here, '..', 'PSR_data');
if ~exist(P.resDir, 'dir'), mkdir(P.resDir); end
cfg = minfuel_config();

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------

% ---- WHICH t_f factors to sweep -------------------------------------------
% Either an explicit vector, e.g.  factors = [1.12 1.15 1.20 1.30];
% or the string 'energy' to run EVERY factor that has an energy seed on disk
% (sundman_minfuel/results/energy/energy_f####.mat -- currently the usable
% band ~1.12x..1.95x). 'energy' is the "run all t_f's we have seeds for" mode.
factors = 'energy';              % vector of factors, or 'energy' (all seeds)

% ---- homotopy endpoint (same meaning as run_psr's epsMin) ------------------
% 0 (default) -> bang-bang; >0 -> smoother eps-optimal control. Applied to
% every factor in the sweep. epsMin>0 auto-skips the bang-bang-only stages
% (refinement, PMP verify) exactly as in run_psr.
epsMin = 0;

% ---- seed policy -----------------------------------------------------------
% 'energy' (recommended for a sweep) -- each factor seeds from its OWN energy
% backbone. (Other run_psr seed modes exist but chaining them across a sweep is
% a different workflow; 'energy' is the clean per-factor batch mode.)
seedSpec = 'energy';

% ---- knobs (same as run_psr) ----------------------------------------------
P.sched      = [];               % [] = campaign default schedule, truncated at epsMin
P.maxIter    = cfg.maxIter;
P.refineOpts = struct('maxRounds',4,'K',8,'maxAdd',40);
P.verifyOpts = struct('M',40,'mode','d','epsEval',1e-4, ...
                      'adjArcs',[],'adjSwitches',[],'makeFig',false);
% MOVIE per factor: 'none' | 'preview' (3 fast stills) | 'movie' (full MP4+GIF,
% ~15 min each). Default 'none' -- a full-movie batch is hours; set 'movie' to
% render one per factor (same control movie run_psr makes), 'preview' for quick
% stills.
P.movieMode  = 'none';
% runVerify: also run the first-order PMP verifier (verify_direct_pmp) and append
% its `verify` summary to each bang-bang data product, so the .mat matches
% run_psr's exactly (the core psr_data product + ipoptCert are always written).
% Off by default -- verify is ~5-10 min/factor. (No effect for epsMin>0.)
P.runVerify  = false;
P.rerunDirect = false;  P.rerunRefine = false;

% pass sweep-wide settings into P (used by the per-factor runner below)
P.epsMin = epsMin;  P.seedSpec = seedSpec;  P.cfg = cfg;

%% ------------------------------------------------------------------------
%% 2. Resolve the factor list
%% ------------------------------------------------------------------------
if ischar(factors) || isstring(factors)
    assert(strcmpi(factors,'energy'), 'factors must be a vector or ''energy''');
    d = dir(fullfile(cfg.dirs.energy, 'energy_f*.mat'));
    fv = arrayfun(@(e) cfg.fparse(e.name), d);
    factors = sort(fv(~isnan(fv)));
    fprintf('run_psr_batch: %d factors with energy seeds: %s\n', numel(factors), mat2str(factors,4));
else
    factors = sort(factors(:).');
end
assert(~isempty(factors), 'no factors to run');

%% ------------------------------------------------------------------------
%% 3. Sweep
%% ------------------------------------------------------------------------
nF = numel(factors);
res(nF) = struct('factor',[],'ok',[],'dV',[],'prop',[],'switches',[], ...
    'edge',[],'defect',[],'certLocalMin',[],'dataFile',[],'err','');
fprintf('\n=== PSR BATCH: %d factors, epsMin=%.3g, seed=%s ===\n', nF, epsMin, seedSpec);
for i = 1:nF
    f = factors(i);
    fprintf('\n########## [%d/%d] factor %.3f ##########\n', i, nF, f);
    r = struct('factor',f,'ok',false,'dV',NaN,'prop',NaN,'switches',NaN, ...
               'edge',NaN,'defect',NaN,'certLocalMin',NaN,'dataFile','','err','');
    try
        r = run_one_factor(f, P, r);
    catch ME
        r.err = ME.message;
        fprintf('  FACTOR %.3f FAILED: %s\n', f, ME.message);
    end
    res(i) = r;
end

%% ------------------------------------------------------------------------
%% 4. Summary
%% ------------------------------------------------------------------------
fprintf('\n=== PSR BATCH SUMMARY (epsMin=%.3g) ===\n', epsMin);
fprintf('%-8s %-4s %-9s %-9s %-4s %-7s %-9s %-8s %s\n', ...
    'factor','ok','dV(km/s)','prop(kg)','sw','edge%','defect','certLM','note');
for i = 1:nF
    r = res(i);
    fprintf('%-8.3f %-4d %-9.4f %-9.4f %-4d %-7.1f %-9.1e %-8s %s\n', ...
        r.factor, r.ok, r.dV, r.prop, r.switches, 100*r.edge, r.defect, ...
        certstr(r.certLocalMin), r.err);
end
sumFile = fullfile(P.dataDir, sprintf('psr_batch_summary_minEps%s.mat', ...
                   strrep(sprintf('%g',epsMin),'.','p')));
if ~exist(P.dataDir,'dir'), mkdir(P.dataDir); end
save(sumFile, 'res', 'factors', 'epsMin', 'seedSpec');
fprintf('\nsaved summary: %s\n', sumFile);

% ============================================================================
function r = run_one_factor(factor, P, r)
% One factor through the PSR pipeline (mirrors run_psr stages 2-6). Returns the
% result row r. Resumable: stages skip when their output exists.
cfg = P.cfg;
eTag = strrep(sprintf('%g', P.epsMin), '.', 'p');
tag  = sprintf('f%04d_minEps%s', round(1000*factor), eTag);
directFile  = fullfile(P.resDir, ['psr_direct_'  tag '.mat']);
seedFile    = fullfile(P.resDir, ['psr_seed_'    tag '.mat']);
refinedFile = fullfile(P.resDir, ['psr_refined_' tag '.mat']);
bangBang = (P.epsMin == 0);

% --- stage 2: direct solve (energy->fuel homotopy, endpoint epsMin) ---------
if isfile(directFile) && ~P.rerunDirect
    D = load(directFile);  outDirect = D.out;
    fprintf('  [2] direct exists -- skip\n');
else
    if isempty(P.sched)
        if strcmp(P.seedSpec,'energy'), base = cfg.schedSharpen; else, base = cfg.schedNeighbor; end
    else
        base = P.sched;
    end
    effSched = [base(base > P.epsMin), P.epsMin];
    fprintf('  [2] direct solve (endpoint eps=%.3g)...\n', P.epsMin);
    outDirect = minfuel_at_tf(factor, 'seed', P.seedSpec, 'sched', effSched, ...
        'outFile', directFile, 'maxIter', P.maxIter, 'branch', 'psr');
    assert(outDirect.certified, 'direct solve did not certify');
end
r.dV = outDirect.dV;  r.prop = outDirect.prop_kg;  r.edge = outDirect.edge;
r.defect = outDirect.maxDefect;  r.switches = outDirect.switches;

% --- stage 3: PSR refinement (bang-bang only) -------------------------------
if bangBang
    if isfile(refinedFile) && ~P.rerunRefine
        finalSol = refinedFile;  fprintf('  [3] refined exists -- skip\n');
    else
        fprintf('  [3] PSR refinement...\n');
        prep_refine_seed(directFile, seedFile);
        ro = P.refineOpts;  ro.tag = ['psr_' tag];  ro.outDir = P.resDir;  ro.solFile = refinedFile;
        h = refine_loop(seedFile, ro);
        if ~isempty(h), r.switches = h(end).switches; end
        finalSol = refinedFile;
    end
else
    finalSol = directFile;  fprintf('  [3] smooth eps>0 -- refinement N/A\n');
end

% --- stage 4: costate + data export -----------------------------------------
fprintf('  [4] export...\n');
r.dataFile = psr_export_data(finalSol, P.dataDir, struct('M',P.verifyOpts.M,'epsMin',P.epsMin,'quiet',true));

% --- stage 5: first-order PMP verify (bang-bang only; optional, matches run_psr)
if bangBang && P.runVerify
    fprintf('  [5] first-order PMP verify...\n');
    oldDir = cd(P.resDir);  cu = onCleanup(@() cd(oldDir));   % verifier writes to pwd
    verify = verify_direct_pmp(finalSol, P.verifyOpts); %#ok<NASGU>
    clear cu
    save(r.dataFile, 'verify', '-append');
end

% --- stage 5b: IPOPT-native local-min certificate ---------------------------
fprintf('  [5b] IPOPT local-min certificate...\n');
ic = psr_ipopt_certify(finalSol, struct('eps',P.epsMin,'verbose',false));
r.certLocalMin = ic.certLocalMin;
ipoptCert = ic; %#ok<NASGU>
save(r.dataFile, 'ipoptCert', '-append');

% --- stage 6: movie (optional) ----------------------------------------------
if ~strcmp(P.movieMode,'none')
    fprintf('  [6] movie (%s)...\n', P.movieMode);
    ttl = sprintf('PSR GTO\\rightarrowtulip t_f=%.2fx (eps=%.3g)', factor, P.epsMin);
    psr_movie(finalSol, fullfile(P.resDir, ['psr_movie_' tag]), ttl, P.movieMode);
end

r.ok = true;
fprintf('  DONE factor %.3f: dV=%.4f prop=%.4f sw=%d certLM=%s\n', ...
        factor, r.dV, r.prop, r.switches, certstr(r.certLocalMin));
end

% ============================================================================
function s = certstr(c)
if isnan(c), s = '-'; elseif c, s = 'YES'; else, s = 'no'; end
end
