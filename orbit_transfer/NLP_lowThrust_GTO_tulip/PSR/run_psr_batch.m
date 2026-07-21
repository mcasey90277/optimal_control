% RUN_PSR_BATCH  Sweep the PSR pipeline over many t_f factors (in one process).
%
% The in-process batch sibling of run_psr.m: run the full PSR pipeline
% (psr_run_one) for a VECTOR of t_f factors, or factors='energy' to auto-run
% every factor that has a min-energy seed on disk, at a single homotopy endpoint
% epsMin. Each factor writes the SAME data product as run_psr
% (PSR_data/psr_data_tf####_minEps#.mat, incl. the ipoptCert), optionally a
% movie, and a per-factor result row; a summary table prints + saves at the end.
%
% *** CRASH NOTE ***  This runs every factor in ONE MATLAB process, so the
% sporadic UNCATCHABLE CasADi/IPOPT MEX FATAL crash (~1 in 10 solves) kills the
% WHOLE sweep (a try/catch cannot catch a MEX FATAL). For a long or unattended
% sweep use the shell walker **psr_batch.sh** instead -- it runs each factor in
% its own process so a crash kills only that factor. Both are resumable (stages
% skip when outputs exist) and share psr_run_one / psr_collect_summary.
%
% Edit section 1, then run. (`./psr_batch.sh <epsMin> energy` is the robust
% terminal equivalent.)

%% ------------------------------------------------------------------------
%% 0. Paths
%% ------------------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
dataDir = fullfile(here, '..', 'PSR_data');
cfg = minfuel_config();

%% ------------------------------------------------------------------------
%% 1. PARAMETERS  (edit this section only)
%% ------------------------------------------------------------------------

% ---- WHICH t_f factors to sweep -------------------------------------------
% Either an explicit vector, e.g.  factors = [1.12 1.15 1.20 1.30];
% or the string 'energy' to run EVERY factor that has an energy seed on disk
% (currently the usable band ~1.12x..1.95x). 'energy' = "all t_f's we have seeds for".
factors = 'energy';

% ---- homotopy endpoint (same meaning as run_psr's epsMin) ------------------
% 0 (default) -> bang-bang; >0 -> smoother eps-optimal control. Applied to every
% factor; epsMin>0 auto-skips the bang-bang-only stages, exactly as in run_psr.
epsMin = 0;

% ---- other params (same meaning as run_psr) -------------------------------
seedSpec  = 'energy';            % per-factor seed ('energy' is the batch mode)
movieMode = 'none';              % 'none' | 'preview' | 'movie'  (per factor)
runVerify = false;               % also append first-order PMP verify (slow)
rerunDirect = false;  rerunRefine = false;

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
    factors = sort(double(factors(:)).');
end
assert(~isempty(factors), 'no factors to run');

%% ------------------------------------------------------------------------
%% 3. Sweep  (psr_run_one per factor; try/catch keeps a CATCHABLE failure from
%%           killing the sweep -- but a MEX FATAL still kills all; use the shell)
%% ------------------------------------------------------------------------
opts = struct('epsMin',epsMin, 'seedSpec',seedSpec, 'movieMode',movieMode, ...
    'runVerify',runVerify, 'rerunDirect',rerunDirect, 'rerunRefine',rerunRefine, ...
    'dataDir',dataDir);
eTag = strrep(sprintf('%g', epsMin), '.', 'p');
nF = numel(factors);
fprintf('\n=== PSR BATCH (in-process): %d factors, epsMin=%.3g ===\n', nF, epsMin);
for i = 1:nF
    f = factors(i);
    fprintf('\n########## [%d/%d] factor %.3f ##########\n', i, nF, f);
    try
        psr_run_one(f, opts);
    catch ME
        fprintf('  FACTOR %.3f FAILED: %s\n', f, ME.message);
        tag = sprintf('f%04d_minEps%s', round(1000*f), eTag);
        row = struct('factor',f,'ok',false,'dV',NaN,'prop',NaN,'switches',NaN, ...
            'edge',NaN,'defect',NaN,'certLocalMin',NaN,'dataFile','','err',ME.message); %#ok<NASGU>
        save(fullfile(dataDir, sprintf('psr_result_%s.mat', tag)), 'row');
    end
end

%% ------------------------------------------------------------------------
%% 4. Summary
%% ------------------------------------------------------------------------
psr_collect_summary(epsMin, dataDir);
