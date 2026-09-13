function out = run_costate_library(opts)
%% Purpose:
%
%   THE ENTRY SCRIPT for building a DRO-to-tulip costate library at a given
%   resolution. One call, from the operating point and the grid, to a
%   packaged and audited catalog.
%
%     run_costate_library(struct('nD', 24, 'nA', 24))
%
%   It runs the campaign the way doc/CAMPAIGN_DISCIPLINE.md says to, because
%   the 24x24 library was built the other way and lost about a day:
%
%     1. ARRIVAL SHEET from the saved arcs, re-scanned at this grid's levels
%        (crossings_from_arc). Refining the arrival axis needs no new
%        continuation -- the arcs already stored every crossing there is.
%     2. CALIBRATE: walk ONE departure column and measure it. Every budget
%        and watchdog below is sized from that measurement rather than from
%        a guess. The 24x24 plan was built on 25-40 s per certification when
%        the real figure was 2.5 minutes, and every downstream number was
%        wrong by that factor.
%     3. RIBS through a work_queue: one column per unit, workers claim the
%        next free one. Never a fixed list per worker -- columns run from 90
%        minutes to over five hours, so a static split leaves workers idle
%        behind their own slow column and loses whole lists when one dies.
%     4. PACKAGE, AUDIT, SWEEP via the existing chain, with blockers rather
%        than aborts.
%
%   Stages are switchable and every one is resumable from its artifacts, so
%   a re-run continues rather than restarting.
%
%  ASSUMPTIONS / NOTES:
%
% • The coarse grid is a SUBSET of a doubled one (k/12 = 2k/24), so raising
%   the resolution ADDS cells; existing certified entries are reused.
% • Rib workers are separate MATLAB processes, launched and verified by
%   run_campaign_workers.sh. This function can launch them (.launch) or
%   print the command for a human to run.
%
%% Inputs:
%
%  opts                     struct (optional)
%   .nD [24] .nA [24]          the grid
%   .sA0 [0.0754]              arrival phase origin
%   .outDir ['results_fine']   everything this run writes
%   .nWorkers [4]              rib workers
%   .run                       stage switches: .sheet .calibrate .ribs
%                              .package .audit .sweep (all true)
%   .launch [false]            actually spawn the rib workers
%   .unitSec []                skip calibration and use this
%
%% Outputs:
%
%  out                      struct                  .sheet .queue .catalog
%                                                   .unitSec .cmd (the
%                                                   launch command) .blockers
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
nD = d('nD', 24);  nA = d('nA', 24);  sA0 = d('sA0', 0.0754);
outDir = d('outDir', fullfile(here, 'results_fine'));
nWorkers = d('nWorkers', 4);
run_ = d('run', struct());
on = @(f) ~isfield(run_, f) || run_.(f);
if ~isfolder(outDir), mkdir(outDir); end

out = struct('sheet', '', 'queue', fullfile(outDir, 'ribq'), 'catalog', '', ...
             'unitSec', NaN, 'cmd', '', 'blockers', {{}});
sheetMat = fullfile(outDir, sprintf('arrival_sheet_70mN_nA%d.mat', nA));
out.sheet = sheetMat;
fprintf('COSTATE LIBRARY: %d departure x %d arrival, out %s\n', nD, nA, outDir);

%% 1. ARRIVAL SHEET -- from the saved arcs, no new continuation
if on('sheet') && ~isfile(sheetMat)
    build_arrival_sheet(struct('nA', nA, 'sA0', sA0, 'rescan', true, 'out', sheetMat));
end
assert(isfile(sheetMat), 'run_costate_library:noSheet', ...
       'the arrival sheet is missing: %s (stage 1 makes it)', sheetMat);
L_ = load(sheetMat);  S = L_.S;
cols = find(isfinite(S.TF));
safe_report(@() fprintf('1. sheet: %d of %d arrival phases certified\n', numel(cols), numel(S.sA)), 'sheet');

%% 2. CALIBRATE -- measure ONE unit before sizing anything
ribOut = @(j) fullfile(outDir, sprintf('fine_rib_col%02d.mat', j));
out.unitSec = d('unitSec', NaN);
if on('calibrate') && ~isfinite(out.unitSec)
    doneAlready = cols(arrayfun(@(j) isfile(ribOut(j)), cols));
    if ~isempty(doneAlready)
        % a column already on disk means the campaign has run before; take
        % the conservative default rather than re-walking one to time it
        out.unitSec = 3600;
        fprintf('2. calibration: %d column(s) already on disk; unit assumed %.0f s\n', ...
                numel(doneAlready), out.unitSec);
    else
        j0 = cols(1);
        fprintf('2. calibration: walking column %d to measure one unit ...\n', j0);
        tC = tic;
        build_ribs(sheetMat, struct('only', j0, 'direction', -1, 'nD', nD, ...
                                    'nPts', nD-1, 'wallSec', 900, 'out', ribOut(j0)));
        out.unitSec = toc(tC);
        fprintf('2. one column took %.0f s; budgets sized from THIS, not a guess\n', out.unitSec);
    end
end
if ~isfinite(out.unitSec), out.unitSec = 3600; end

%% 3. RIBS -- a work queue, one column per unit
if on('ribs')
    work_queue('init', out.queue, cols(:).', ribOut);
    st = work_queue('status', out.queue);
    fprintf('3. rib queue: %d unit(s), %d already done\n', st.n, st.nDone);
    jobFile = fullfile(outDir, 'rib_unit_job.m');
    writeRibJob(jobFile, sheetMat, outDir, out.queue, nD);
    out.cmd = sprintf('%s %s %d %d %s', ...
        fullfile(fileparts(fileparts(here)), 'costate_common', 'run_campaign_workers.sh'), ...
        jobFile, nWorkers, round(out.unitSec), outDir);
    if st.nTodo > 0
        if d('launch', false)
            fprintf('3. launching %d worker(s) ...\n', nWorkers);
            system(sprintf('%s', out.cmd));
        else
            fprintf(['3. %d unit(s) to walk. Launch them with:\n\n    %s\n\n' ...
                     '   then re-run this function to package.\n'], st.nTodo, out.cmd);
            return
        end
    end
end

%% 4. PACKAGE / AUDIT / SWEEP -- the existing chain, blockers not aborts
if on('package') || on('audit') || on('sweep')
    co = struct('outDir', outDir, 'run', struct( ...
        'arcs', false, 'sheet', false, 'ribs', false, ...
        'package', on('package'), 'audit', on('audit'), 'sweep', on('sweep'), ...
        'pictures', true, 'deliverable', false));
    fprintf('4. package / audit / sweep via build_70mN_library ...\n');
    % the chain is a SCRIPT and reads chainOverrides from the caller's
    % workspace; base is the only workspace a function can hand it one in
    assignin('base', 'chainOverrides', co);
    evalin('base', 'build_70mN_library');
end
out.catalog = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
end

% ------------------------------------------------------------------------
function writeRibJob(jobFile, sheetMat, outDir, qDir, nD)
% WRITERIBJOB  Emit the per-worker job script the launcher runs. It reads
% WORKER_TAG from the environment and walks whatever column the queue hands
% it.  INPUTS: jobFile; sheetMat; outDir; qDir; nD.  OUTPUTS: none.
txt = sprintf([ ...
 '%%%% RIB_UNIT_JOB  One rib worker. Written by run_costate_library.\n' ...
 'here = pwd; cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(here);\n' ...
 'addpath(''%s'', ''%s'');\n' ...
 'capped_pool(1);\n' ...
 'tag = getenv(''WORKER_TAG''); if isempty(tag), tag = ''w0''; end\n' ...
 'ribOut = @(j) fullfile(''%s'', sprintf(''fine_rib_col%%02d.mat'', j));\n' ...
 'unitFcn = @(j, beat) build_ribs(''%s'', struct(''only'', j, ''direction'', -1, ...\n' ...
 '        ''nD'', %d, ''nPts'', %d, ''wallSec'', 900, ''out'', ribOut(j)));\n' ...
 'campaign_worker(''%s'', ''%s'', tag, unitFcn, struct(''logFile'', ...\n' ...
 '        fullfile(''%s'', [''worker_'' tag ''.log''])));\n'], ...
 fileparts(sheetMat), fullfile(fileparts(fileparts(fileparts(sheetMat))), 'costate_common'), ...
 outDir, sheetMat, nD, nD-1, qDir, fullfile(outDir, 'hb'), outDir);
fid = fopen(jobFile, 'w');  fprintf(fid, '%s', txt);  fclose(fid);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
