function res = psr_collect_summary(epsMin, dataDir)
% PSR_COLLECT_SUMMARY  Build the PSR batch summary from per-factor result rows.
%
% Scans <dataDir>/psr_result_f####_minEps<e>.mat (each holding one `row` struct
% saved by psr_run_one) for the given epsMin, prints a summary table, and saves
% <dataDir>/psr_batch_summary_minEps<e>_<insertionLabel>.mat. Used by
% psr_batch.sh (after a crash-robust per-process sweep) and by run_psr_batch.m.
%
% Each `row` (since Task 4) carries rv0, rvf, insertion (= insMeta.label); rows
% saved before Task 4 lack these and are backfilled (NaN(1,6) / '') so older
% cached rows still aggregate. The summary's top-level rv0/rvf/insertion are
% taken from the first row (a batch shares one declared criterion).
%
% INPUTS:
%   epsMin  - homotopy endpoint the sweep ran at [scalar]
%   dataDir - PSR_data directory [default ../PSR_data relative to this file]
% OUTPUTS:
%   res - [1xN] struct array of the collected rows, sorted by factor

here = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(dataDir), dataDir = fullfile(here, '..', 'PSR_data'); end
eTag = strrep(sprintf('%g', epsMin), '.', 'p');

d = dir(fullfile(dataDir, sprintf('psr_result_f*_minEps%s.mat', eTag)));
res = struct('factor',{},'ok',{},'dV',{},'prop',{},'switches',{}, ...
    'edge',{},'defect',{},'certLocalMin',{},'dataFile',{},'err',{}, ...
    'rv0',{},'rvf',{},'insertion',{});
for k = 1:numel(d)
    S = load(fullfile(dataDir, d(k).name), 'row');
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
    fprintf('psr_collect_summary: no result rows for epsMin=%.3g in %s\n', epsMin, dataDir);
    return
end
[~, ord] = sort([res.factor]);  res = res(ord);

fprintf('\n=== PSR BATCH SUMMARY (epsMin=%.3g, %d factors) ===\n', epsMin, numel(res));
fprintf('%-8s %-4s %-9s %-9s %-4s %-7s %-9s %-8s %s\n', ...
    'factor','ok','dV(km/s)','prop(kg)','sw','edge%','defect','certLM','note');
for k = 1:numel(res)
    r = res(k);
    if isnan(r.certLocalMin), lm='-'; elseif r.certLocalMin, lm='YES'; else, lm='no'; end
    fprintf('%-8.3f %-4d %-9.4f %-9.4f %-4d %-7.1f %-9.1e %-8s %s\n', ...
        r.factor, r.ok, r.dV, r.prop, r.switches, 100*r.edge, r.defect, lm, r.err);
end
% endpoints/insertion for the summary: taken from the first row (all rows in a
% batch share one declared criterion in the current design; empty/NaN for a
% batch built entirely from pre-Task-4 legacy rows).
rv0 = res(1).rv0;  rvf = res(1).rvf;  insertion = res(1).insertion; %#ok<NASGU>
insTag = insertion;  if isempty(insTag), insTag = 'unknown'; end
sumFile = fullfile(dataDir, sprintf('psr_batch_summary_minEps%s_%s.mat', eTag, insTag));
save(sumFile, 'res', 'epsMin', 'rv0', 'rvf', 'insertion');
fprintf('\nsaved summary: %s\n', sumFile);
end
