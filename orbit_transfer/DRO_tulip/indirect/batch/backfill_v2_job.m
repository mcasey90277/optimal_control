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
%   verdict  -- reads every audit chunk and rewrites the verdict, which
%               starts with ONE decision token (backfill_clean_verdict):
%               CLEAN only if audit nBad == 0, coverage complete, primNot4
%               and movedPrimaries empty, content key equal, and the
%               comparison with the record shows no change; else NOT CLEAN
%               with the failing conditions named.
%   audit    -- one audit chunk (CHUNK, NCHUNK from the environment); what
%               run_recertify.sh launches with its second argument 'audit'.
%   relabel  -- after a NOT CLEAN audit whose only problems are statuses that
%               flip at a threshold: relabel_borderline on every audit chunk's
%               altRows (global index .k). Refuses unless every audit_alt_*.mat
%               is a final save bound to THIS catalog (both content keys), no
%               alternative is audited twice, and no pre-relabel copy or
%               audit_round1/ exists yet. Copies the catalog to
%               costate_catalog_dro_tulip_70mN_prerelabel.mat, saves the
%               relabelled one back (same file, same variable name), moves
%               audit_*.mat, audit_*_batch.out and audit_driver.log into
%               audit_round1/ (they no longer match the alternatives key),
%               and writes ~/BACKFILL_V2_RELABEL_VERDICT.txt. Then re-run the
%               audit (borderline rows: the stored status is a lower bound).
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
%          comparison_with_record, costate_catalog_dro_tulip_70mN}.mat;
%          relabel: <outDir>/costate_catalog_dro_tulip_70mN_prerelabel.mat,
%          <outDir>/audit_round1/
stage = lower(strtrim(getenv('BACKFILL_STAGE')));  if isempty(stage), stage = 'harvest'; end
ind = fileparts(fileparts(mfilename('fullpath')));
RES = fullfile(ind, 'results');
recordMat = fullfile(RES, 'library_70mN_24x48_merged', 'costate_catalog_dro_tulip_70mN.mat');
outDir = strtrim(getenv('BACKFILL_OUT'));  if isempty(outDir), outDir = fullfile(RES, 'library_70mN_24x48_v2'); end
verdictF = fullfile(getenv('HOME'), 'BACKFILL_V2_VERDICT.txt');
if strcmp(stage, 'audit'), verdictF = fullfile(getenv('HOME'), sprintf('BACKFILL_V2_AUDIT_%s_VERDICT.txt', getenv('CHUNK'))); end
if strcmp(stage, 'relabel'), verdictF = fullfile(getenv('HOME'), 'BACKFILL_V2_RELABEL_VERDICT.txt'); end

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
            % each chunk file must be the FINAL save (not crash salvage) of an
            % audit of THIS catalog: both content keys must match, or it is stale
            for f = dir(fullfile(outDir, 'audit_alt_*.mat'))'
                Q = load(fullfile(outDir, f.name));
                assert(isfield(Q, 'audit'), 'backfill_v2_job:auditChunk', '%s is a partial (crash salvage): rerun that chunk', f.name);
                assert(strcmp(Q.audit.contentKey, catalog_content_key(c2)) && strcmp(Q.audit.altContentKey, alternatives_content_key(c2)), ...
                       'backfill_v2_job:auditKey', '%s audited another catalog', f.name);
                nOk = nOk + Q.audit.nOk;  nBad = nBad + Q.audit.nBad;
                seen = [seen, [Q.audit.altRows.k]];  nPrim = nPrim + numel(Q.audit.primRows);
            end
            % COVERAGE: every alternative audited exactly once, every primary audited
            nAlt = numel(c2.alternatives);  nPrimCat = nnz(c2.sheets(1).has_solution);
            dupl = unique(seen(arrayfun(@(x) nnz(seen == x) > 1, seen)));
            gap = setdiff(1:nAlt, seen);
            coverageOk = numel(seen) == nAlt && numel(unique(seen)) == nAlt && nPrim == nPrimCat;
            if ~coverageOk
                error('backfill_v2_job:coverage', ['audit coverage incomplete: %d alternative row(s) audited for %d ' ...
                      'alternatives (%d never, e.g. %s; %d more than once, e.g. %s); %d primaries audited of %d'], ...
                      numel(seen), nAlt, numel(gap), mat2str(gap(1:min(5, end))), numel(dupl), ...
                      mat2str(dupl(1:min(5, end))), nPrim, nPrimCat);
            end
            Lr = load(recordMat);  rec = Lr.(char(fieldnames(Lr)));
            V = struct('nBad', nBad, 'coverageOk', coverageOk, 'primNot4', c2.status_layer.primNot4, ...
                       'movedPrimaries', c2.status_layer.movedPrimaries, ...
                       'contentKeyEqual', strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'cmp', Cm.C);
            [token, whyNot] = backfill_clean_verdict(V);
            msg = ['BACKFILL VERDICT: ' token tern(isempty(whyNot), '', [' (' whyNot ')']) ' | ' summaryLine(c2, recordMat, Cm.C) ...
                   sprintf(' | audit %d ok / %d bad | primaries audited %d | alternatives audited %d of %d', ...
                           nOk, nBad, nPrim, numel(unique(seen)), nAlt)];

        case 'relabel'
            catMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
            preMat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN_prerelabel.mat');
            arch = fullfile(outDir, 'audit_round1');
            assert(~isfile(preMat), 'backfill_v2_job:relabelled', '%s exists: this catalog was relabelled already', preMat);
            assert(~isfolder(arch), 'backfill_v2_job:archive', '%s exists: refusing to mix audit rounds', arch);
            L = load(catMat);  vn = char(fieldnames(L));  c2 = L.(vn);
            ck = catalog_content_key(c2);  ak = alternatives_content_key(c2);
            chunks = dir(fullfile(outDir, 'audit_alt_*.mat'));
            assert(~isempty(chunks), 'backfill_v2_job:noAudit', 'no audit_alt_*.mat in %s: run the audit first', outDir);
            rows = struct('k', {}, 'ok', {}, 'why', {}, 'statusNow', {}, 'moved', {});
            for f = chunks'
                Q = load(fullfile(outDir, f.name));
                assert(isfield(Q, 'audit'), 'backfill_v2_job:auditChunk', '%s is a partial (crash salvage): rerun that chunk', f.name);
                assert(strcmp(Q.audit.contentKey, ck) && strcmp(Q.audit.altContentKey, ak), ...
                       'backfill_v2_job:auditKey', '%s audited another catalog: refusing to relabel from it', f.name);
                ar = Q.audit.altRows;
                if ~isempty(ar), rows = [rows, rmfield(ar, setdiff(fieldnames(ar), fieldnames(rows)))]; end
            end
            seen = [rows.k];
            assert(numel(unique(seen)) == numel(seen), 'backfill_v2_job:duplicate', 'an alternative is audited in two chunks');
            [c3, info] = relabel_borderline(c2, rows);
            assert(copyfile(catMat, preMat), 'backfill_v2_job:copy', 'could not copy %s', catMat);
            S = struct();  S.(vn) = c3;  save(catMat, '-struct', 'S');
            mkdir(arch);
            moved = [dir(fullfile(outDir, 'audit_*.mat')); dir(fullfile(outDir, 'audit_*_batch.out')); dir(fullfile(outDir, 'audit_driver.log'))];
            for f = moved', movefile(fullfile(outDir, f.name), fullfile(arch, f.name)); end
            st = [c3.alternatives.status];
            rl = strjoin(arrayfun(@(m) sprintf('k%d %d->%d (re-audit %d)', info.rows(m, 1), info.rows(m, 2), info.rows(m, 4), ...
                                              info.rows(m, 3)), 1:size(info.rows, 1), 'UniformOutput', false), ', ');
            nr = strjoin(arrayfun(@(x) sprintf('k%d (%s)', x.k, x.why), info.notRelabelled, 'UniformOutput', false), '; ');
            msg = sprintf(['BACKFILL RELABEL: %d borderline of %d alternatives [%s] | NOT relabelled %d%s | status counts ' ...
                           '4: %d, 3: %d, 2: %d, 1: %d | pre-relabel copy %s | %d audit file(s) moved to audit_round1/ | ' ...
                           'next: batch/run_recertify.sh N audit, then BACKFILL_STAGE=verdict'], ...
                          size(info.rows, 1), numel(st), rl, numel(info.notRelabelled), tern(isempty(nr), '', [' [' nr ']']), ...
                          nnz(st == 4), nnz(st == 3), nnz(st == 2), nnz(st == 1), preMat, numel(moved));

        otherwise
            error('backfill_v2_job:stage', 'BACKFILL_STAGE must be harvest | assemble | audit | verdict | relabel, got ''%s''', stage);
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
             'primaries moved %d, primaries not 4 on re-polish %d, re-certification errors %d | content key equal: %s | ' ...
             'vs record: %d of %d agree, ' ...
             'missing %d, slower %d, faster %d'], ...
            nnz(c2.sheets(1).has_solution), numel(st), nnz(st == 4), nnz(st == 3), nnz(st == 2), nnz(st == 1), ...
            c2.status_layer.nRecert, c2.status_layer.nMoved, numel(c2.status_layer.movedPrimaries), ...
            numel(c2.status_layer.primNot4), c2.status_layer.nRecertErrors, tern(strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'yes', 'NO'), ...
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
