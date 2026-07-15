function row = psr_run_one(factor, opts)
% PSR_RUN_ONE  Run the PSR pipeline for ONE t_f factor (the batch unit of work).
%
% The per-factor pipeline extracted so it can be called (a) in a loop by
% run_psr_batch.m, and (b) as a standalone `matlab -batch` process by
% psr_batch.sh -- the crash-robust shell walker. Running each factor in its own
% process is what lets the sporadic UNCATCHABLE CasADi/IPOPT MEX FATAL crash
% (~1 in 10 solves) kill only that factor instead of the whole sweep (a MATLAB
% try/catch cannot catch a MEX FATAL). Mirrors run_psr.m's stages 2-6.
%
% INPUTS:
%   factor - t_f / t_f_min [scalar]
%   opts   - (optional) struct; all fields have run_psr-matching defaults:
%            epsMin      homotopy endpoint [0]  (0=bang-bang, >0=smooth)
%            seedSpec    'energy' | 'neighbor' | <path>  ['energy']
%            movieMode   'none' | 'preview' | 'movie'  ['none']
%            runVerify   also append first-order PMP verify (bang-bang) [false]
%            maxIter     IPOPT cap  [cfg.maxIter]
%            sched       eps schedule override [] (default, truncated at epsMin)
%            refineOpts  struct maxRounds/K/maxAdd [4/8/40]
%            verifyOpts  struct for verify_direct_pmp [M=40,mode='d',...]
%            rerunDirect / rerunRefine  force those stages [false]
%            resDir      intermediates dir  [PSR/results]
%            dataDir     data-product dir   [../PSR_data]
%
% OUTPUTS:
%   row - result struct: factor, ok, dV, prop, switches, edge, defect,
%         certLocalMin, dataFile, err. ALSO saved to
%         <dataDir>/psr_result_f####_minEps#.mat (var `row`) so the shell walker
%         can collect a summary after crashes.
%
% REFERENCES: run_psr.m (single-factor driver), run_psr_batch.m (in-process
%   sweep), psr_batch.sh (crash-robust per-process sweep).

here = fileparts(mfilename('fullpath'));
addpath(here);  setup_paths();
addpath(fullfile(here, '..', 'sundman_minfuel'));   % insertion_states (single-source; PSR vendors the rest)
cfg = minfuel_config();

if nargin < 2, opts = struct(); end
d = @(f,v) getfield_default(opts, f, v);
epsMin     = d('epsMin', 0);
seedSpec   = d('seedSpec', 'energy');
movieMode  = d('movieMode', 'none');
runVerify  = d('runVerify', false);
maxIter    = d('maxIter', cfg.maxIter);
sched      = d('sched', []);
refineOpts = d('refineOpts', struct('maxRounds',4,'K',8,'maxAdd',40));
verifyOpts = d('verifyOpts', struct('M',40,'mode','d','epsEval',1e-4, ...
                    'adjArcs',[],'adjSwitches',[],'makeFig',false));
rerunDirect = d('rerunDirect', false);
rerunRefine = d('rerunRefine', false);
resDir     = d('resDir', fullfile(here,'results'));
dataDir    = d('dataDir', fullfile(here,'..','PSR_data'));
if ~exist(resDir,'dir'), mkdir(resDir); end
if ~exist(dataDir,'dir'), mkdir(dataDir); end

% ---- INSERTION POINT (edit here to retarget) ---------------------------------
insertion = 'campaign';        % tulip: 'campaign'|'maxydot'|'apoapsis'  (elfo: 'nearest'|'apolune'|'perilune')
% insertion = 'maxydot';       % uncomment to use the max-ydot point (needs a matching energy seed)
% insertion = 'apoapsis';      % uncomment to use the slowest/apoapsis point (needs a matching seed)
[rv0, rvf, insMeta] = insertion_states('tulip', insertion);   % <TGT> = 'tulip' or 'elfo'

eTag = strrep(sprintf('%g', epsMin), '.', 'p');
tag  = sprintf('f%04d_minEps%s', round(1000*factor), eTag);
directFile  = fullfile(resDir, ['psr_direct_'  tag '.mat']);
seedFile    = fullfile(resDir, ['psr_seed_'    tag '.mat']);
refinedFile = fullfile(resDir, ['psr_refined_' tag '.mat']);
bangBang    = (epsMin == 0);

row = struct('factor',factor,'ok',false,'dV',NaN,'prop',NaN,'switches',NaN, ...
    'edge',NaN,'defect',NaN,'certLocalMin',NaN,'dataFile','','err','');
resultFile = fullfile(dataDir, sprintf('psr_result_%s.mat', tag));

fprintf('\n=== PSR_RUN_ONE factor %.3f (epsMin=%.3g, seed=%s) ===\n', factor, epsMin, char(string(seedSpec)));

% --- stage 2: direct solve --------------------------------------------------
if isfile(directFile) && ~rerunDirect
    D = load(directFile);  outDirect = D.out;  fprintf('[2] direct exists -- skip\n');
else
    if isempty(sched)
        if strcmp(seedSpec,'energy'), base = cfg.schedSharpen; else, base = cfg.schedNeighbor; end
    else
        base = sched;
    end
    effSched = [base(base > epsMin), epsMin];
    fprintf('[2] direct solve (endpoint eps=%.3g)...\n', epsMin);
    % drift guard: locate the same seed minfuel_at_tf is about to load and
    % confirm it matches the declared insertion point ('energy' / explicit-file
    % seedSpec cases; 'neighbor' isn't wired through this driver's opts, so it
    % is left to minfuel_at_tf's own error -- see run_psr.m for the full case).
    if strcmpi(seedSpec, 'energy')
        guardSeedFile = fullfile(cfg.dirs.energy, cfg.fname('energy', factor));
    elseif ~strcmpi(seedSpec, 'neighbor') && ischar(seedSpec) && isfile(seedSpec)
        guardSeedFile = seedSpec;
    else
        guardSeedFile = '';
    end
    if ~isempty(guardSeedFile) && isfile(guardSeedFile)
        S = load(guardSeedFile, 'rvf', 'rv0');
        assert(norm(S.rvf(:).' - rvf) < 1e-10 && norm(S.rv0(:).' - rv0) < 1e-10, ...
            'insertion:drift', ['seed endpoints differ from the declared %s insertion ' ...
            '(rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
            insMeta.label, norm(S.rvf(:).'-rvf), norm(S.rv0(:).'-rv0));
    end
    outDirect = minfuel_at_tf(factor, 'seed', seedSpec, 'sched', effSched, ...
        'outFile', directFile, 'maxIter', maxIter, 'branch', 'psr');
    assert(outDirect.certified, 'direct solve did not certify at factor %.3f', factor);
end
row.dV = outDirect.dV;  row.prop = outDirect.prop_kg;  row.edge = outDirect.edge;
row.defect = outDirect.maxDefect;  row.switches = outDirect.switches;

% --- stage 3: PSR refinement (bang-bang only) -------------------------------
if bangBang
    if isfile(refinedFile) && ~rerunRefine
        finalSol = refinedFile;  fprintf('[3] refined exists -- skip\n');
    else
        fprintf('[3] PSR refinement...\n');
        prep_refine_seed(directFile, seedFile);
        ro = refineOpts;  ro.tag = ['psr_' tag];  ro.outDir = resDir;  ro.solFile = refinedFile;
        h = refine_loop(seedFile, ro);
        if ~isempty(h), row.switches = h(end).switches; end
        finalSol = refinedFile;
    end
else
    finalSol = directFile;  fprintf('[3] smooth eps>0 -- refinement N/A\n');
end

% --- stage 4: costate + data export -----------------------------------------
fprintf('[4] export...\n');
row.dataFile = psr_export_data(finalSol, dataDir, struct('M',verifyOpts.M,'epsMin',epsMin,'quiet',true));

% --- stage 5: first-order PMP verify (bang-bang only; optional) -------------
if bangBang && runVerify
    fprintf('[5] first-order PMP verify...\n');
    oldDir = cd(resDir);  cu = onCleanup(@() cd(oldDir));
    verify = verify_direct_pmp(finalSol, verifyOpts); %#ok<NASGU>
    clear cu
    save(row.dataFile, 'verify', '-append');
end

% --- stage 5b: IPOPT-native local-min certificate ---------------------------
fprintf('[5b] IPOPT local-min certificate...\n');
ic = psr_ipopt_certify(finalSol, struct('eps',epsMin,'verbose',false));
row.certLocalMin = ic.certLocalMin;
ipoptCert = ic; %#ok<NASGU>
save(row.dataFile, 'ipoptCert', '-append');

% --- stage 6: movie (optional) ----------------------------------------------
if ~strcmp(movieMode,'none')
    fprintf('[6] movie (%s)...\n', movieMode);
    ttl = sprintf('PSR GTO\\rightarrowtulip t_f=%.2fx (eps=%.3g)', factor, epsMin);
    psr_movie(finalSol, fullfile(resDir, ['psr_movie_' tag]), ttl, movieMode);
end

row.ok = true;
save(resultFile, 'row');
fprintf('=== DONE factor %.3f: dV=%.4f prop=%.4f sw=%d certLM=%d  (row -> %s) ===\n', ...
        factor, row.dV, row.prop, row.switches, row.certLocalMin, resultFile);
end

% ---------------------------------------------------------------------------
function v = getfield_default(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
