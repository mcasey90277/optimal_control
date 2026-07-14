function row = elfo_run_one(factor, opts)
% ELFO_RUN_ONE  Run the GTO->ELFO min-fuel homotopy for ONE t_f factor (batch unit).
%
% The per-factor unit of work for the ELFO min-fuel tf-grid campaign, extracted
% so it can be called (a) as a standalone `matlab -batch` process by elfo_batch.sh
% (each factor in its OWN process, so an UNCATCHABLE CasADi/IPOPT MEX FATAL crash
% kills only that factor, not the sweep -- a try/catch cannot catch a MEX FATAL),
% and (b) directly. The ELFO analog of PSR/psr_run_one.m; wraps gen_elfo_minfuel.
%
% INPUTS:
%   factor - t_f / tfMin (tfMin = tulip min-time 6.2907 ND, a shared scale) [scalar]
%   opts   - (optional) struct:
%            .epsMin    homotopy endpoint [0]  (0 = bang-bang fuel, >0 = smooth)
%            .maxIter   IPOPT cap (tight) [2000]
%            .looseIter IPOPT cap (loose probe) [500]
%            .resDir    seeds + results dir [elfo/results]
%            .rerun     ignore an existing result row and re-solve [false]
%
% OUTPUTS:
%   row - [1x1] struct: factor, tf, tf_days, ok, epsReached, epsFloor, dV, prop,
%         switches, edge, defect, ipoptStatus, dataFile, err. ALSO saved to
%         <resDir>/elfo_result_f####_minEps#.mat (var `row`) so the shell walker
%         can build a summary after crashes.
%
% REFERENCES:
%   [1] gen_elfo_minfuel.m (homotopy core); [2] run_elfo_minfuel.m (interactive
%       entry); [3] PSR/psr_run_one.m (tulip analog); [4] elfo_batch.sh.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

if nargin < 2, opts = struct(); end
d = @(f,v) getdef(opts, f, v);
epsMin    = d('epsMin', 0);
maxIter   = d('maxIter', 2000);
looseIter = d('looseIter', 500);
resDir    = d('resDir', fullfile(here,'results'));
rerun     = d('rerun', false);
if ~exist(resDir,'dir'), mkdir(resDir); end

tf   = factor * cfg.tfMin;
eTag = strrep(sprintf('%g', epsMin), '.', 'p');
tag  = sprintf('f%04d_minEps%s', round(1000*factor), eTag);
resultFile = fullfile(resDir, sprintf('elfo_result_%s.mat', tag));

row = struct('factor',factor,'tf',tf,'tf_days',tf*p.tStar/86400,'ok',false, ...
    'epsReached',false,'epsFloor',NaN,'dV',NaN,'prop',NaN,'switches',NaN, ...
    'edge',NaN,'defect',NaN,'ipoptStatus','','dataFile','','err','');

% --- resumable: an existing row means this factor is done --------------------
if isfile(resultFile) && ~rerun
    L = load(resultFile, 'row');
    if isfield(L,'row'), row = L.row; fprintf('elfo_run_one f=%.3f: row exists -- skip\n', factor); return; end
end

% --- resolve the factor-keyed energy seed (base-seed fallback near 1.20x) ----
seed = fullfile(resDir, sprintf('energy_elfo_f%04d.mat', round(1000*factor)));
if ~isfile(seed)
    base = fullfile(resDir, 'energy_elfo_freetf.mat');
    if isfile(base)
        B = load(base,'X');  if abs(B.X(8,end) - tf) < 0.02, seed = base; end
    end
end
if ~isfile(seed)
    row.err = sprintf('no ELFO energy seed for factor %.3f (build via gen_elfo_energy_tfsweep)', factor);
    save(resultFile,'row');  fprintf('elfo_run_one f=%.3f: %s\n', factor, row.err);  return
end

fprintf('\n=== ELFO_RUN_ONE factor %.3f (epsMin=%g, seed=%s) ===\n', factor, epsMin, seed);

% --- solve: energy->fuel homotopy (gen_elfo_minfuel errors id 'minfuel:stuck'
%     at the sharpening wall -- catch it and record the eps-floor) -----------
try
    outFile = gen_elfo_minfuel(struct('seedFile',seed,'target','ELFO','epsMin',epsMin, ...
        'maxIter',maxIter,'looseIter',looseIter,'resume',~rerun));
    L  = load(outFile);
    ss = L.U(4,:);  mf = L.X(7,end);
    row.ok         = true;
    row.epsReached = (L.epsilon <= epsMin + 1e-9);
    row.switches   = sum(abs(diff(ss>0.5)));
    row.edge       = mean(ss>0.95 | ss<0.05);
    row.dV         = p.c*log(1/mf)*p.lStar/p.tStar;
    row.prop       = p.m0kg*(1-mf);
    row.defect     = L.out.maxDefect;
    if isfield(L.out,'ipoptStatus'), row.ipoptStatus = L.out.ipoptStatus; end
    row.dataFile   = outFile;
catch ME
    row.err = ME.message;
    if strcmp(ME.identifier, 'minfuel:stuck')
        row.epsReached = false;
        tok = regexp(ME.message, 'eps=([0-9.]+)', 'tokens', 'once');
        if ~isempty(tok), row.epsFloor = str2double(tok{1}); end
        fprintf('elfo_run_one f=%.3f: sharpening wall at eps=%.4f\n', factor, row.epsFloor);
    else
        fprintf('elfo_run_one f=%.3f: ERROR %s\n', factor, ME.message);
    end
end

save(resultFile, 'row');
fprintf('=== elfo_run_one f=%.3f: ok=%d epsReached=%d dV=%.4f sw=%d (row -> %s) ===\n', ...
        factor, row.ok, row.epsReached, row.dV, row.switches, resultFile);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
