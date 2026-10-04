%% RECERTIFY_CHUNK_JOB  (2026-10-04)  One chunk of the status-layer backfill's
% re-certification (recertify_candidates): the legacy verdict-0 candidates,
% the rib stop points and the unmatched primaries listed by
% backfill_v2_job's harvest, items CHUNK:NCHUNK:end, under the current
% certifier. Saves after every item and resumes: run it again after a kill.
%
% Launched by batch/run_recertify.sh N (one process per chunk), never in the
% interactive session. Paths follow this file's location.
%
% INPUTS:  none (environment: CHUNK, NCHUNK, BACKFILL_OUT [results/
%          library_70mN_24x48_v2])
% OUTPUTS: ~/RECERT_<CHUNK>_VERDICT.txt (one line); <outDir>/recert_<CHUNK>.mat
%          (+ .log)
k = str2double(getenv('CHUNK'));  n = str2double(getenv('NCHUNK'));
verdictF = fullfile(getenv('HOME'), sprintf('RECERT_%s_VERDICT.txt', getenv('CHUNK')));
ind = fileparts(fileparts(mfilename('fullpath')));
outDir = strtrim(getenv('BACKFILL_OUT'));
if isempty(outDir), outDir = fullfile(ind, 'results', 'library_70mN_24x48_v2'); end

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    assert(isfinite(k) && isfinite(n) && k >= 1 && k <= n, 'recertify_chunk_job:chunk', ...
           'set CHUNK and NCHUNK (1 <= CHUNK <= NCHUNK), got "%s" of "%s"', getenv('CHUNK'), getenv('NCHUNK'));
    t0 = tic;
    items = recertify_candidates(fullfile(outDir, 'harvest.mat'), k, n, fullfile(outDir, sprintf('recert_%d.mat', k)), ...
                                 struct('pool', capped_pool(1)));
    st = arrayfun(@statusOf, items);
    msg = sprintf(['RECERT CHUNK %d/%d FINISHED: %d items | threw %d | moved %d | status 4: %d, 3: %d, 2: %d, 1: %d, ' ...
                   'below floor: %d | %.1f h'], k, n, numel(items), nnz(~cellfun(@isempty, {items.err})), nnz([items.moved]), ...
                  nnz(st == 4), nnz(st == 3), nnz(st == 2), nnz(st == 1), nnz(st < 1), toc(t0)/3600);
catch ME
    msg = sprintf('RECERT CHUNK %s FAILED: %s | %s', getenv('CHUNK'), ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);

function s = statusOf(it)
% STATUSOF  An item's status (NaN when it threw or carries none).
% INPUTS: it (one item).  OUTPUTS: s (double).
s = NaN;
if isempty(it.err) && isfield(it.C, 'status') && ~isempty(it.C.status), s = double(it.C.status); end
end
