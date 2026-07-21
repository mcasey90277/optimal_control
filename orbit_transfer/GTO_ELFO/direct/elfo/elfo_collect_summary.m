function res = elfo_collect_summary(epsMin, resDir)
% ELFO_COLLECT_SUMMARY  Build the ELFO min-fuel tf-grid map from per-factor rows.
%
% Scans <resDir>/elfo_result_f####_minEps<e>.mat (each holding one `row` struct
% saved by elfo_run_one) for the given epsMin, prints the tf-grid convergence
% table (which factors reach eps=0 vs stall, switches, the dV-time front), and
% saves <resDir>/elfo_batch_summary_minEps<e>_<insertionLabel>.mat. Used by
% elfo_batch.sh. The ELFO analog of PSR/psr_collect_summary.m.
%
% Each `row` (since Task 4) carries rv0, rvf, insertion (= insMeta.label); rows
% saved before Task 4 lack these and are backfilled (NaN(1,6) / '') so older
% cached rows still aggregate. The summary's top-level rv0/rvf/insertion are
% taken from the first row (a batch shares one declared criterion).
%
% INPUTS:
%   epsMin - homotopy endpoint the sweep ran at [scalar]
%   resDir - elfo results directory [default elfo/results relative to this file]
% OUTPUTS:
%   res - [1xN] struct array of the collected rows, sorted by factor

here = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(resDir), resDir = fullfile(here, 'results'); end
eTag = strrep(sprintf('%g', epsMin), '.', 'p');

d = dir(fullfile(resDir, sprintf('elfo_result_f*_minEps%s.mat', eTag)));
res = struct('factor',{},'tf',{},'tf_days',{},'ok',{},'epsReached',{},'epsFloor',{}, ...
    'dV',{},'prop',{},'switches',{},'edge',{},'defect',{},'ipoptStatus',{},'dataFile',{},'err',{}, ...
    'rv0',{},'rvf',{},'insertion',{});
for k = 1:numel(d)
    S = load(fullfile(resDir, d(k).name), 'row');
    if isfield(S,'row')
        r = S.row;
        % backward-compat: rows saved before Task 4 lack rv0/rvf/insertion --
        % backfill so the struct-array concatenation below does not error.
        if ~isfield(r,'rv0'),       r.rv0 = NaN(1,6);   end
        if ~isfield(r,'rvf'),       r.rvf = NaN(1,6);   end
        if ~isfield(r,'insertion'), r.insertion = '';   end
        res(end+1) = r; %#ok<AGROW>
    end
end
if isempty(res)
    fprintf('elfo_collect_summary: no result rows for epsMin=%.3g in %s\n', epsMin, resDir);
    return
end
[~, ord] = sort([res.factor]);  res = res(ord);

fprintf('\n=== ELFO MIN-FUEL tf-GRID MAP (epsMin=%.3g, %d factors) ===\n', epsMin, numel(res));
fprintf('%-7s %-6s %-4s %-8s %-9s %-9s %-4s %-7s %-9s %s\n', ...
    'factor','tf(d)','ok','epsRch','dV(km/s)','prop(kg)','sw','edge%','defect','note');
for k = 1:numel(res)
    r = res(k);
    if ~r.ok
        er = '-';
    elseif r.epsReached
        er = 'YES';
    else
        er = sprintf('no@%.3g', r.epsFloor);
    end
    fprintf('%-7.3f %-6.2f %-4d %-8s %-9.4f %-9.4f %-4d %-7.1f %-9.1e %s\n', ...
        r.factor, r.tf_days, r.ok, er, r.dV, r.prop, r.switches, 100*r.edge, r.defect, r.err);
end
% endpoints/insertion for the summary: taken from the first row (all rows in a
% batch share one declared criterion in the current design; empty/NaN for a
% batch built entirely from pre-Task-4 legacy rows).
rv0 = res(1).rv0;  rvf = res(1).rvf;  insertion = res(1).insertion; %#ok<NASGU>
insTag = insertion;  if isempty(insTag), insTag = 'unknown'; end
sumFile = fullfile(resDir, sprintf('elfo_batch_summary_minEps%s_%s.mat', eTag, insTag));
save(sumFile, 'res', 'epsMin', 'rv0', 'rvf', 'insertion');
fprintf('\nsaved summary: %s\n', sumFile);
end
