function out = run_costate_library(opts)
%% Purpose:
%
%   THE ENTRY SCRIPT for building a DRO-to-tulip costate library at a given
%   resolution. One call, from the operating point and the grid, to a
%   packaged and audited catalog.
%
%     run_costate_library(struct('nD', 24, 'nA', 24))          % by resolution
%     run_costate_library(struct('sD', (0:23)/24, ...           % or by phase
%                                'sA', 0.0754 + (0:23)/24))
%
%   THE PHASES ARE CHOSEN IN ONE PLACE, section 0 below, and printed before
%   anything runs. Either give the two phase VECTORS outright, or give a
%   resolution and let them be derived. Both axes are resolved together --
%   the departure grid used to be set by nD deep inside the rib stage while
%   the arrival grid was set by nA in the sheet stage, so there was no one
%   place to look.
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
%   THE ORBITS AND THE ENGINE -- what the library is OF:
%   .tauDRO [1.0]              departure DRO period, ND (sets its size)
%   .NpTulip [7]               arrival tulip petal count (SETS ITS PERIOD:
%                              tau = 2*pi*(Np-2)/(Np-1), it is not free)
%   .pmTulip [-1]              tulip branch, -1 is the southern one
%   .thrustN [0.070] .ispS [900] .m0kg [150]      the engine
%
%   THE GRID -- give EITHER the phases or the resolution:
%   .sD  [(0:nD-1)/nD]         departure phases, fractions of the DRO period
%   .sA  [sA0 + (0:nA-1)/nA]   arrival phases, fractions of the tulip period
%   .nD [24] .nA [24]          resolution, used only when .sD/.sA are absent
%   .sA0 [0.0754]              arrival origin, used only to derive .sA
%   .onlyA []                  walk ribs on these arrival phases only
%                              (values or indices into .sA; default all
%                              certified ones)
%   .outDir ['results_fine']   everything this run writes (made absolute)
%   .nWorkers [4]              rib workers
%   .run                       stage switches: .sheet .calibrate .ribs
%                              .package .audit .sweep (all true)
%   .launch [false]            actually spawn the rib workers
%   .maxAtt [3] .staleSec [1800]   the ONE ownership/retry policy, handed
%                              to the queue, the workers and the monitor
%
%  LIFECYCLE. This function is safe to call repeatedly while workers run:
%  re-opening the queue keeps their claims and the attempt counts (the
%  first version wiped both, so "re-run to package" would have freed
%  columns still being walked). It returns out.state:
%    'launched'  workers were spawned; call again later to package
%    'pending'   units remain and nothing was launched (.launch false, or
%                workers are still on them); the command is in out.cmd
%    'blocked'   units were RETIRED after repeated failure; the library
%                cannot be complete until they are reset or dropped
%    'packaged'  every unit's artifact exists and stage 4 ran
%  Packaging NEVER runs while a unit is open: a rib file being written is
%  not a rib.
%   .unitSec []                skip calibration and use this
%
%% Outputs:
%
%  out                      struct                  .state (above) .sD .sA
%                                                   .orbits .engine .sheet
%                                                   .queue .unitSec .cmd
%                                                   .blockers (why it is not
%                                                   'packaged') .catalog
%                                                   .sheet .queue .catalog
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
outDir = d('outDir', fullfile(here, 'results_fine'));
nWorkers = d('nWorkers', 4);
run_ = d('run', struct());
on = @(f) ~isfield(run_, f) || run_.(f);
if ~isfolder(outDir), mkdir(outDir); end
% ABSOLUTE from here on: the generated worker job is executed by run(),
% which changes directory to the job's folder, so a relative outDir would
% mean something else inside the worker
outDir = absPath(outDir);
maxAtt = d('maxAtt', 3);  staleSec = d('staleSec', 1800);

%% ========================================================================
%  0. WHAT THIS LIBRARY IS OF, and over which phases. The ONLY place the
%     orbits, the engine and the phases are chosen -- all four used to be
%     implicit: the orbits and engine came from arclength_arrival's own
%     defaults, the arrival grid from nA in the sheet stage, the departure
%     grid from nD deep inside the rib stage.
%     A phase is a fraction of its own orbit's period: sD of the DRO's, sA
%     of the tulip's. The builders underneath need UNIFORM grids, so an
%     explicit vector is checked for that and refused by name if it is not.
%% ========================================================================
orbits = struct('tauDRO', d('tauDRO', 1.0), 'NpTulip', d('NpTulip', 7), ...
                'pmTulip', d('pmTulip', -1));
engine = struct('thrustN', d('thrustN', 0.070), 'ispS', d('ispS', 900), 'm0kg', d('m0kg', 150));
tauTulip = 2*pi*(orbits.NpTulip - 2)/(orbits.NpTulip - 1);      % LOCKED by Np
sA0 = d('sA0', 0.0754);
nD = d('nD', 24);  nA = d('nA', 24);
sD = d('sD', (0:nD-1)/nD);
sA = d('sA', sA0 + (0:nA-1)/nA);
sD = sD(:).';  sA = sA(:).';
nD = numel(sD);  nA = numel(sA);
checkUniform(sD, 'departure (.sD)');
checkUniform(sA, 'arrival (.sA)');
sA0 = sA(1);

out = struct('state', 'pending', 'sheet', '', 'queue', fullfile(outDir, 'ribq'), ...
             'catalog', '', 'unitSec', NaN, 'cmd', '', 'blockers', {{}});
sheetMat = fullfile(outDir, sprintf('arrival_sheet_70mN_nA%d.mat', nA));
out.sheet = sheetMat;
tStar = 382981.289129055;
fprintf('COSTATE LIBRARY\n');
fprintf('  departure    : DRO, tau = %.4f ND (%.3f d)\n', orbits.tauDRO, orbits.tauDRO*tStar/86400);
fprintf('  arrival      : %d-petal tulip, branch %+d, tau = %.4f ND (%.3f d, LOCKED by Np)\n', ...
        orbits.NpTulip, orbits.pmTulip, tauTulip, tauTulip*tStar/86400);
fprintf('  engine       : %.0f mN, Isp %g s, %g kg\n', engine.thrustN*1000, engine.ispS, engine.m0kg);
fprintf('  %d departure x %d arrival phases\n', nD, nA);
fprintf('  departure sD : %s\n', phaseList(sD));
fprintf('  arrival   sA : %s\n', phaseList(sA));
fprintf('  writing to   : %s\n', outDir);
out.sD = sD;  out.sA = sA;  out.orbits = orbits;  out.engine = engine;

%% 1. ARRIVAL SHEET -- from the saved arcs, no new continuation
if on('sheet') && ~isfile(sheetMat)
    build_arrival_sheet(struct('nA', nA, 'sA0', sA0, 'rescan', true, 'out', sheetMat, ...
        'thrustN', engine.thrustN, 'ispS', engine.ispS, 'm0kg', engine.m0kg, ...
        'tauDRO', orbits.tauDRO, 'NpTulip', orbits.NpTulip, 'pmTulip', orbits.pmTulip));
end
assert(isfile(sheetMat), 'run_costate_library:noSheet', ...
       'the arrival sheet is missing: %s (stage 1 makes it)', sheetMat);
L_ = load(sheetMat);  S = L_.S;
% THE SHEET ON DISK DECIDES: a reused sheet built at another operating point
% would silently make this a library of something else
if isfield(S, 'problem')
    P = S.problem;
    wantP = [engine.thrustN, engine.ispS, engine.m0kg, orbits.tauDRO, orbits.NpTulip, orbits.pmTulip];
    haveP = [P.thrustN, P.ispS, P.m0kg, P.tauDRO, P.NpTulip, P.pmTulip];
    assert(all(abs(haveP - wantP) <= 1e-9*max(1, abs(wantP))), 'run_costate_library:identity', ...
        ['the arrival sheet on disk is a DIFFERENT problem.\n' ...
         '  on disk   : %.1f mN, Isp %g s, %g kg, DRO tau %g, %d-petal tulip branch %+d\n' ...
         '  asked for : %.1f mN, Isp %g s, %g kg, DRO tau %g, %d-petal tulip branch %+d\n' ...
         'Point .outDir somewhere else, or delete %s.'], ...
        haveP(1)*1000, haveP(2), haveP(3), haveP(4), haveP(5), haveP(6), ...
        wantP(1)*1000, wantP(2), wantP(3), wantP(4), wantP(5), wantP(6), sheetMat);
end
cols = find(isfinite(S.TF));
onlyA = d('onlyA', []);
if ~isempty(onlyA)
    if all(onlyA == round(onlyA)) && all(onlyA >= 1) && all(onlyA <= numel(S.sA))
        want = onlyA(:).';                       % indices
    else
        want = arrayfun(@(v) find(abs(S.sA - mod(v,1)) < 1e-9, 1), onlyA(:).'); % values
    end
    cols = intersect(cols, want);
    fprintf('   .onlyA restricts the ribs to arrival phases %s\n', phaseList(S.sA(cols)));
end
safe_report(@() fprintf('1. sheet: %d of %d arrival phases certified; %d will get ribs\n', ...
                        nnz(isfinite(S.TF)), numel(S.sA), numel(cols)), 'sheet');
missingA = setdiff(1:numel(S.sA), find(isfinite(S.TF)));
if ~isempty(missingA)
    % these columns are part of the REQUESTED grid and will be absent from
    % the library; say so here rather than letting them vanish from the count
    out.blockers{end+1} = sprintf('%d arrival phase(s) have no certified seed and get no rib: %s', ...
                                  numel(missingA), phaseList(S.sA(missingA)));
end
if isempty(cols)
    out.state = 'blocked';
    out.blockers{end+1} = 'no arrival column is eligible for ribs (none certified, or .onlyA excluded them all)';
    fprintf('BLOCKED: %s\n', out.blockers{end});
    return
end

%% 2. CALIBRATE -- measure ONE unit before sizing anything
ribOut = @(j) fullfile(outDir, sprintf('fine_rib_col%02d.mat', j));
out.unitSec = d('unitSec', NaN);
if on('calibrate') && ~isfinite(out.unitSec)
    doneAlready = cols(arrayfun(@(j) isfile(ribOut(j)), cols));
    calMat = fullfile(outDir, 'calibration.mat');
    if isfile(calMat)
        cal = load(calMat);
        out.unitSec = cal.unitSec;
        fprintf('2. calibration: column %d measured %.0f s on %s (reused)\n', cal.column, cal.unitSec, cal.measured);
    elseif ~isempty(doneAlready)
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
        % persisted with what it measured: a later call must not replace a
        % measurement with the 3600 s guess just because a rib file exists
        cal = struct('unitSec', out.unitSec, 'column', j0, 'nD', nD, 'measured', char(datetime('now')));
        save(fullfile(outDir, 'calibration.mat'), '-struct', 'cal');
        fprintf('2. one column took %.0f s; budgets sized from THIS, not a guess\n', out.unitSec);
    end
end
if ~isfinite(out.unitSec), out.unitSec = 3600; end

%% 3. RIBS -- a work queue, one column per unit
policy = struct('staleSec', staleSec, 'maxAtt', maxAtt);
codeRoots = {here, fullfile(fileparts(fileparts(here)), 'costate_common')};
if on('ribs')
    % 'init' OPENS an existing queue (claims and attempts kept) and only
    % creates when there is none -- see work_queue
    work_queue('init', out.queue, cols(:).', ribOut);
    st = work_queue('status', out.queue, policy);
    fprintf('3. rib queue: %d unit(s): %d done, %d running, %d stale, %d retired, %d to do\n', ...
            st.n, st.nDone, st.nRunning, st.nStale, st.nRetired, st.nTodo);
    if st.nRetired > 0
        out.blockers{end+1} = sprintf('%d column(s) RETIRED after %d failed attempts: %s (work_queue reset, or drop them)', ...
                                      st.nRetired, maxAtt, mat2str(st.retired));
    end
    jobFile = fullfile(outDir, 'rib_unit_job.m');
    writeRibJob(jobFile, sheetMat, outDir, out.queue, nD, codeRoots, policy);
    launcher = fullfile(codeRoots{2}, 'run_campaign_workers.sh');
    out.cmd = strjoin({shq(launcher), shq(jobFile), num2str(nWorkers), ...
                       num2str(round(out.unitSec)), shq(outDir)}, ' ');
    if st.nOpen > 0
        if d('launch', false)
            fprintf('3. launching %d worker(s) ...\n', nWorkers);
            [rc, txt] = system(out.cmd);
            fprintf('%s', txt);
            if rc ~= 0
                out.state = 'blocked';
                out.blockers{end+1} = sprintf('launcher exited %d -- see above; nothing packaged', rc);
                return
            end
            out.state = 'launched';
            fprintf('3. workers launched. Call this function again when they finish to package.\n');
        else
            out.state = 'pending';
            fprintf(['3. %d unit(s) to walk. Launch them with:\n\n    %s\n\n' ...
                     '   then call this function again to package.\n'], st.nOpen, out.cmd);
        end
        return
    end
    if st.nRunning > 0
        out.state = 'pending';
        fprintf('3. %d column(s) still being walked by live workers: %s. Packaging waits.\n', ...
                st.nRunning, mat2str(st.running));
        return
    end
    if ~st.complete
        out.state = 'blocked';
        fprintf('3. BLOCKED: %s\n', strjoin(out.blockers, '\n   '));
        return
    end
end

%% 4. PACKAGE / AUDIT / SWEEP -- the existing chain, blockers not aborts
if on('package') || on('audit') || on('sweep')
    co = struct('outDir', outDir, 'run', struct( ...
        'arcs', false, 'sheet', false, 'ribs', false, ...
        'package', on('package'), 'audit', on('audit'), 'sweep', on('sweep'), ...
        'pictures', d('pictures', true), 'deliverable', false));
    fprintf('4. package / audit / sweep via build_70mN_library ...\n');
    % the chain is a SCRIPT that reads chainOverrides from the workspace it
    % runs in; base is the only workspace a function can hand it one in.
    % Whatever was there before is put back afterwards, success or not.
    prev = [];  hadPrev = evalin('base', 'exist(''chainOverrides'', ''var'')') == 1;
    if hadPrev, prev = evalin('base', 'chainOverrides'); end
    restore = onCleanup(@() restoreBase(hadPrev, prev));
    assignin('base', 'chainOverrides', co);
    evalin('base', 'build_70mN_library');
    if evalin('base', 'exist(''chainBlockers'', ''var'')') == 1
        cb = evalin('base', 'chainBlockers');
        if iscell(cb), out.blockers = [out.blockers, cb(:).']; end
    end
end
cat = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
if isfile(cat) && on('package')
    out.catalog = cat;  out.state = 'packaged';
else
    out.catalog = '';
    if on('package'), out.blockers{end+1} = sprintf('no catalog was written to %s', cat); end
end
if ~isempty(out.blockers)
    fprintf('BLOCKERS:\n');  fprintf('  - %s\n', out.blockers{:});
end
end

% ------------------------------------------------------------------------
function restoreBase(hadPrev, prev)
% RESTOREBASE  Put the caller's chainOverrides back (or remove ours).
% INPUTS: hadPrev logical; prev.  OUTPUTS: none.
if hadPrev, assignin('base', 'chainOverrides', prev);
else,       evalin('base', 'clear chainOverrides'); end
end

% ------------------------------------------------------------------------
function p = absPath(p)
% ABSPATH  Absolute, normalised path.  INPUTS: p.  OUTPUTS: p.
if ~(startsWith(p, '/') || (numel(p) > 1 && p(2) == ':')), p = fullfile(pwd, p); end
p = char(java.io.File(p).getCanonicalPath());
end

% ------------------------------------------------------------------------
function q = shq(s)
% SHQ  Quote one shell argument (POSIX single quotes).  INPUTS: s.
% OUTPUTS: q.
q = ['''' strrep(s, '''', '''\''''') ''''];
end

% ------------------------------------------------------------------------
function q = mlq(s)
% MLQ  Quote one MATLAB char literal.  INPUTS: s.  OUTPUTS: q.
q = ['''' strrep(s, '''', '''''') ''''];
end

% ------------------------------------------------------------------------
function checkUniform(v, name)
% CHECKUNIFORM  The builders underneath assume a uniform phase grid; say so
% by name rather than producing a sheet whose levels do not line up.
% INPUTS: v; name.  OUTPUTS: none.
assert(numel(v) >= 2, 'run_costate_library:grid', '%s grid needs at least two phases', name);
dv = diff(v);
assert(max(abs(dv - dv(1))) < 1e-9, 'run_costate_library:grid', ...
       ['the %s grid must be UNIFORM (the sheet and rib builders step by a ' ...
        'fixed 1/n). Got spacings %.6g .. %.6g'], name, min(dv), max(dv));
end

% ------------------------------------------------------------------------
function s = phaseList(v)
% PHASELIST  A phase vector as text, abbreviated in the middle when long.
% INPUTS: v.  OUTPUTS: s char.
if numel(v) <= 8
    s = strjoin(compose('%.4f', v), ' ');
else
    s = sprintf('%s ... %s  (%d, step %.4f)', strjoin(compose('%.4f', v(1:3)), ' '), ...
                strjoin(compose('%.4f', v(end-1:end)), ' '), numel(v), v(2)-v(1));
end
end

% ------------------------------------------------------------------------
function writeRibJob(jobFile, sheetMat, outDir, qDir, nD, codeRoots, policy)
% WRITERIBJOB  Emit the per-worker job script the launcher runs. It reads
% WORKER_TAG from the environment and walks whatever column the queue hands
% it. Paths are absolute and quoted as MATLAB literals; the code roots are
% the DRIVER'S own, not something inferred from the output directory (the
% first version derived a costate_common that did not exist). The worker's
% heartbeat is passed to the walker as its progress callback, so a claim
% stays fresh for as long as the column takes.
% INPUTS: jobFile; sheetMat; outDir; qDir; nD; codeRoots cell; policy.
% OUTPUTS: none.
roots = strjoin(cellfun(@mlq, codeRoots, 'UniformOutput', false), ', ');
txt = sprintf([ ...
 '%%%% RIB_UNIT_JOB  One rib worker. Written by run_costate_library.\n' ...
 'here = pwd; cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(here);\n' ...
 'addpath(%s);\n' ...
 'capped_pool(1);\n' ...
 'tag = getenv(''WORKER_TAG'');\n' ...
 'assert(~isempty(tag), ''rib_unit_job:tag'', ''WORKER_TAG is not set: run this through run_campaign_workers.sh'');\n' ...
 'ribOut = @(j) fullfile(%s, sprintf(''fine_rib_col%%02d.mat'', j));\n' ...
 'unitFcn = @(j, beat) build_ribs(%s, struct(''only'', j, ''direction'', -1, ...\n' ...
 '        ''nD'', %d, ''nPts'', %d, ''wallSec'', 900, ''out'', ribOut(j), ''progress'', beat));\n' ...
 'campaign_worker(%s, %s, tag, unitFcn, struct(''logFile'', ...\n' ...
 '        fullfile(%s, [''worker_'' tag ''.log'']), ''staleSec'', %d, ''maxAtt'', %d, ''idleSec'', %d));\n'], ...
 roots, mlq(outDir), mlq(sheetMat), nD, nD-1, mlq(qDir), mlq(fullfile(outDir, 'hb')), mlq(outDir), ...
 policy.staleSec, policy.maxAtt, 2*policy.staleSec);
tmp = [jobFile '.part'];
fid = fopen(tmp, 'w');
assert(fid >= 0, 'run_costate_library:job', 'cannot write %s', tmp);
fprintf(fid, '%s', txt);  fclose(fid);
[okMv, msgMv] = movefile(tmp, jobFile);
assert(okMv, 'run_costate_library:job', 'cannot publish %s: %s', jobFile, msgMv);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
