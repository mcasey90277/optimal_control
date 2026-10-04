%% BACKFILL_V2_JOB  (2026-10-04)  The optimality-status layer for the library
% of record WITHOUT a rebuild (spec 6): library_70mN_24x48_merged ->
% library_70mN_24x48_v2. The record is read, never written; adoption is a
% separate decision on this job's verdict (the merge's pattern).
%
% Stages (environment variable BACKFILL_STAGE, default harvest):
%   harvest  -- backfill_status_layer('harvest'): candidates on disk, the
%               primaries' junctions by root, the re-certification work
%               list. No solve; minutes. Prints the counts and the launch
%               line for the chunks.
%   (then)      batch/run_recertify.sh N  -- N recertify_chunk_job processes
%   assemble -- backfill_status_layer('assemble') (refuses while a primary
%               lacks junctions), compare_phase_catalogs against the record,
%               and the audit's PRIMARIES pass (fly every primary from its
%               junctions). Writes the verdict with the alternatives' audit
%               marked pending.
%   (then)      batch/run_recertify.sh N audit  -- the alternatives' audit
%               in N chunks (audit_status_layer .idxAlt = CHUNK:NCHUNK:end),
%               one process each: thousands of full re-certifications do not
%               fit one serial process
%   verdict  -- reads every audit chunk and rewrites the verdict.
%   audit    -- one audit chunk (CHUNK, NCHUNK from the environment); what
%               run_recertify.sh launches with its second argument 'audit'.
%
% BACKFILL_OUT overrides the output folder (default
% results/library_70mN_24x48_v2). Paths follow this file's location, so the
% job runs against the checkout it lives in.
%
% Launch with nohup, never in the interactive session:
%   BACKFILL_STAGE=harvest nohup /Applications/MATLAB_R2026a.app/bin/matlab \
%       -batch "run('<this file>')" > ~/backfill_v2_harvest.out 2>&1 &
%
% INPUTS:  none (environment: BACKFILL_STAGE, BACKFILL_OUT, CHUNK, NCHUNK)
% OUTPUTS: <verdictF> (one line); <outDir>/{harvest, recert_*, audit_*,
%          comparison_with_record, costate_catalog_dro_tulip_70mN}.mat
stage = lower(strtrim(getenv('BACKFILL_STAGE')));  if isempty(stage), stage = 'harvest'; end
ind = fileparts(fileparts(mfilename('fullpath')));
RES = fullfile(ind, 'results');
recordMat = fullfile(RES, 'library_70mN_24x48_merged', 'costate_catalog_dro_tulip_70mN.mat');
outDir = strtrim(getenv('BACKFILL_OUT'));  if isempty(outDir), outDir = fullfile(RES, 'library_70mN_24x48_v2'); end
verdictF = fullfile(getenv('HOME'), 'BACKFILL_V2_VERDICT.txt');
if strcmp(stage, 'audit'), verdictF = fullfile(getenv('HOME'), sprintf('BACKFILL_V2_AUDIT_%s_VERDICT.txt', getenv('CHUNK'))); end

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    switch stage
        case 'harvest'
            r48 = fullfile(RES, 'reproduce_70mN_24x48', 'round_01');
            r24 = fullfile(RES, 'reproduce_70mN_24x24', 'round_01');
            sources = [{fullfile(r48, 'arrival_sheet_70mN_nA48.mat')}, ribFiles(r48), {fullfile(r48, 'fine_rib_direct_holes.mat')}, ...
                       {fullfile(r24, 'arrival_sheet_70mN_nA24.mat')}, ribFiles(r24), {fullfile(r24, 'fine_rib_direct_holes.mat')}, ...
                       {fullfile(RES, 'sheet96_resolution_test', 'arrival_sheet_70mN_nA96.mat')}, ...
                       {fullfile(RES, 'departure_rib96_test', 'rib96_col28.mat')}, ...
                       {fullfile(RES, 'departure_rib96_test', 'rib96_col28_ckpt.mat')}];
            missing = sources(~cellfun(@isfile, sources));
            assert(isempty(missing), 'backfill_v2_job:sources', 'missing source(s): %s', strjoin(missing, ', '));
            t0 = tic;
            H = backfill_status_layer('harvest', recordMat, sources, outDir);
            for k = 1:numel(sources), fprintf('  %5d  %s\n', H.nBySource(k), sources{k}); end
            nItems = numel(H.needRecert) + numel(H.stops) + numel(H.unmatched);
            msg = sprintf(['HARVEST FINISHED: %d candidates from %d sources | needRecert %d | stops %d | unmatched ' ...
                           'primaries %d of %d | %d items to re-certify | %.0f s | next: batch/run_recertify.sh N ' ...
                           '(BACKFILL_OUT=%s)'], numel(H.cands), numel(sources), numel(H.needRecert), numel(H.stops), ...
                          numel(H.unmatched), nnz(~cellfun(@isempty, H.primJ)) + numel(H.unmatched), nItems, toc(t0), outDir);

        case 'assemble'
            c2 = backfill_status_layer('assemble', recordMat, outDir);
            catMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
            C = compare_phase_catalogs(catMat, recordMat);
            save(fullfile(outDir, 'comparison_with_record.mat'), 'C');
            msg = ['BACKFILL ASSEMBLED: ' summaryLine(c2, recordMat, C) ...
                   ' | audit PENDING: batch/run_recertify.sh N audit, then BACKFILL_STAGE=verdict'];

        case 'audit'
            k = str2double(getenv('CHUNK'));  n = str2double(getenv('NCHUNK'));
            assert(isfinite(k) && isfinite(n) && k >= 1 && k <= n, 'backfill_v2_job:chunk', 'set CHUNK and NCHUNK (1 <= CHUNK <= NCHUNK)');
            catMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
            L = load(catMat);  c2 = L.(char(fieldnames(L)));
            idx = k:n:numel(c2.alternatives);
            if isempty(idx) && k ~= 1, error('backfill_v2_job:emptyChunk', 'chunk %d of %d holds no alternative', k, n); end
            % chunk 1 also carries the PRIMARIES pass (audit_status_layer has no
            % "no alternatives" option: an empty .idxAlt means all of them)
            A = audit_status_layer(catMat, struct('idxAlt', idx, 'skipPrimaries', k ~= 1, ...
                                                  'out', fullfile(outDir, sprintf('audit_alt_%d.mat', k)), 'pool', capped_pool()));
            msg = sprintf('BACKFILL AUDIT CHUNK %d/%d: %d alternatives%s, %d ok / %d bad', k, n, numel(idx), ...
                          tern(k == 1, ' + the primaries', ''), A.nOk, A.nBad);

        case 'verdict'
            catMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
            L = load(catMat);  c2 = L.(char(fieldnames(L)));
            Cm = load(fullfile(outDir, 'comparison_with_record.mat'));
            nOk = 0;  nBad = 0;  seen = [];  nPrim = 0;
            for f = dir(fullfile(outDir, 'audit_alt_*.mat'))'
                Q = load(fullfile(outDir, f.name));
                assert(isfield(Q, 'audit'), 'backfill_v2_job:auditChunk', '%s is a partial (crash salvage): rerun that chunk', f.name);
                assert(strcmp(Q.audit.contentKey, catalog_content_key(c2)) && strcmp(Q.audit.altContentKey, alternatives_content_key(c2)), ...
                       'backfill_v2_job:auditKey', '%s audited another catalog', f.name);
                nOk = nOk + Q.audit.nOk;  nBad = nBad + Q.audit.nBad;
                seen = [seen, [Q.audit.altRows.k]];  nPrim = nPrim + numel(Q.audit.primRows);
            end
            nAlt = numel(c2.alternatives);
            msg = ['BACKFILL VERDICT: ' summaryLine(c2, recordMat, Cm.C) ...
                   sprintf(' | audit %d ok / %d bad | primaries audited %d | alternatives audited %d of %d', ...
                           nOk, nBad, nPrim, numel(unique(seen)), nAlt)];

        otherwise
            error('backfill_v2_job:stage', 'BACKFILL_STAGE must be harvest | assemble | audit | verdict, got ''%s''', stage);
    end
catch ME
    msg = sprintf('BACKFILL %s FAILED: %s | %s', upper(stage), ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);

function s = summaryLine(c2, recordMat, C)
% SUMMARYLINE  Entries, alternatives by status, re-certified, moved, content
% key, comparison with the record.
% INPUTS: c2 (v2 catalog); recordMat (char); C (compare_phase_catalogs).
% OUTPUTS: s (char).
L = load(recordMat);  rec = L.(char(fieldnames(L)));
st = [c2.alternatives.status];
s = sprintf(['%d entries | alternatives %d (status 4: %d, 3: %d, 2: %d, 1: %d) | re-certified %d, moved %d, ' ...
             'primaries moved %d, primaries not 4 on re-polish %d | content key equal: %s | vs record: %d of %d agree, ' ...
             'missing %d, slower %d, faster %d'], ...
            nnz(c2.sheets(1).has_solution), numel(st), nnz(st == 4), nnz(st == 3), nnz(st == 2), nnz(st == 1), ...
            c2.status_layer.nRecert, c2.status_layer.nMoved, numel(c2.status_layer.movedPrimaries), ...
            numel(c2.status_layer.primNot4), tern(strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'yes', 'NO'), ...
            C.nAgree, C.nRef, C.nOnlyRef, C.nNewSlower, C.nNewFaster);
end

function f = ribFiles(d)
% RIBFILES  The rib files of a build round, sorted.  INPUTS: d (folder).
% OUTPUTS: f {1 x n} full paths.
L = dir(fullfile(d, 'fine_rib_col*.mat'));
f = sort(fullfile(d, {L.name}));
end

function s = tern(c, a, b)
% TERN  Ternary.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
