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
%   .extraRibFiles {}          rib files from ANOTHER round of this library
%                              (e.g. the previously swept family's) handed
%                              to the packager alongside this round's: per
%                              cell the packager keeps the fastest certified
%                              point, so an earlier family's certified
%                              points survive where the new spine's rib did
%                              not beat them. Validated against this
%                              campaign's problem identity by the packager.
%   .outDir ['results_fine']   everything this run writes (made absolute)
%   .nWorkers [4]              rib workers
%   .run                       stage switches: .sheet .ribs
%                              .package .audit .sweep (all true)
%   .launch [false]            actually spawn the rib workers
%   .maxAtt [3]                attempts after which a column is retired
%   .staleSec [1800]           heartbeat age reported as STALLED (alarm)
%   .hangSec [2700]            heartbeat silence after which the launcher's
%                              supervisor KILLS a worker. The heartbeat
%                              ticks after every capped STAGE of a
%                              certification (certify_root .progress), and
%                              the largest stage cap is 900 s, so 2700 s is
%                              three stages of silence. (Pass 3 found the
%                              beat used to fire once per certification,
%                              which runs seven capped stages, 4500 s.)
%   .foreignPattern            a pgrep pattern; if any such process is
%                              alive the launch is refused (workers from an
%                              older launcher hold no lock, so the queue
%                              cannot see them) ['fine_ribs_range_job']
%
%  LIFECYCLE. Calling this function while workers run is safe: 'open' is
%  read-only (it validates the unit set against the queue on disk), and a
%  column is owned by a kernel lock its worker holds, which no call here
%  can take. It returns out.state:
%    'launched'  the campaign SUPERVISOR was started: it keeps the workers
%                alive and runs the finalize job (this function, packaging
%                stages on) when every column is published
%    'pending'   units remain and nothing was launched (.launch false, or
%                workers are on them); the launch command is in out.cmd
%    'blocked'   units were RETIRED after repeated failure, or foreign
%                workers are alive; nothing was launched or packaged
%    'packaged'  every column's artifact exists, none is held, and THIS
%                call wrote the catalog
%    'failed'    the packaging chain threw; out.blockers has the message
%    'measured'  an audit/sweep-only call ran its stages (no packaging)
%  Packaging runs only when every column is published and no lock is held.
%
%% Outputs:
%
%  out                      struct                  .state (above) .sD .sA
%                                                   .orbits .engine .sheet
%                                                   .queue .cmd .receipt
%                                                   .blockers (why it is not
%                                                   'packaged') .catalog
%                                                   .sheet .queue .catalog
%                                                   .cmd (the
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
maxAtt = d('maxAtt', 3);  staleSec = d('staleSec', 1800);  hangSec = d('hangSec', 2700);
% ONE CONTROLLER AT A TIME. Manifest creation, queue creation, launching
% and packaging are serialised by a lock on the campaign directory held
% for the whole call: two entry-point processes on one directory could
% otherwise both "create" the manifest, both launch, or both package.
% (Astra pass 3, 3.2.) Workers do not take this lock; units have their own.
CL = unit_lock('try', fullfile(outDir, 'campaign.lock'));
assert(CL.held, 'run_costate_library:busy', ...
       'another run_costate_library (or the finish job) is active on %s; wait for it', outDir);
unlockCampaign = onCleanup(@() unit_lock('release', CL.file, CL.token));
invocationId = char(java.util.UUID.randomUUID());

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
% ONLY THE FULL 1/n LATTICE IS BUILT. The builders take n and an origin,
% not a vector, so any other uniform vector (a partial span, a different
% step, a repeated or descending phase) would be printed here and then
% silently replaced by the lattice. Refuse it instead. (Astra pass 2, J.)
assert(abs(sD(2) - sD(1) - 1/nD) < 1e-9 && abs(sA(2) - sA(1) - 1/nA) < 1e-9, ...
       'run_costate_library:grid', ...
       ['only the full periodic lattice is supported: .sD must step by 1/%d and .sA by 1/%d ' ...
        '(got %.6g and %.6g). Change .nD/.nA, not the vectors.'], nD, nA, sD(2)-sD(1), sA(2)-sA(1));
sA0 = sA(1);

out = struct('state', 'pending', 'sheet', '', 'queue', fullfile(outDir, 'ribq'), ...
             'catalog', '', 'cmd', '', 'blockers', {{}});
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
% THE CAMPAIGN MANIFEST: what this output directory IS a library of, written
% on the first call and checked on every later one, so a re-run with other
% orbits, engine, phases or policy cannot quietly adopt the artifacts here.
% (Astra pass 2, ruling on the pushbacks: "one immutable configuration
% file in a campaign-specific directory".)
manifest = struct('orbits', orbits, 'engine', engine, 'sD', sD, 'sA', sA, ...
                  'policy', struct('maxAtt', maxAtt, 'staleSec', staleSec, 'hangSec', hangSec), ...
                  'codeRevision', gitRevision(here), 'created', char(datetime('now')));
manFile = fullfile(outDir, 'campaign_manifest.mat');
if isfile(manFile)
    M = load(manFile);
    for f_ = {'orbits', 'engine', 'sD', 'sA'}
        assert(isequaln(M.(f_{1}), manifest.(f_{1})), 'run_costate_library:manifest', ...
               ['%s differs from the campaign manifest in %s (created %s). This directory is a ' ...
                'library of ONE problem and grid; use another .outDir.'], f_{1}, outDir, M.created);
    end
    if ~isequal(M.policy, manifest.policy)
        fprintf('   NOTE: policy differs from the manifest (maxAtt/staleSec/hangSec); the manifest''s stands for this campaign\n');
        maxAtt = M.policy.maxAtt;  staleSec = M.policy.staleSec;  hangSec = M.policy.hangSec;
    end
    if ~strcmp(M.codeRevision, manifest.codeRevision)
        fprintf('   NOTE: code revision now %s, manifest was written at %s\n', manifest.codeRevision, M.codeRevision);
    end
else
    tmpM = sprintf('%s.%s.part', manFile, char(java.util.UUID.randomUUID()));
    save(tmpM, '-struct', 'manifest');  publish_atomic(tmpM, manFile);
    fprintf('  manifest     : written (%s)\n', manFile);
end

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
assert(isfield(S, 'problem') && isstruct(S.problem), 'run_costate_library:identity', ...
       'the sheet %s carries no problem identity; rebuild it with build_arrival_sheet', sheetMat);
if true
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
% THE SHEET MUST BE THE GRID ASKED FOR: its arrival phases are the requested
% .sA, and it was certified at the requested departure origin
wrapd = @(x) abs(mod(x + 0.5, 1) - 0.5);
assert(numel(S.sA) == nA && max(wrapd(S.sA(:).' - sA)) < 1e-9, 'run_costate_library:grid', ...
       'the sheet %s holds arrival phases %s, not the requested %s', sheetMat, ...
       phaseList(S.sA(:).'), phaseList(sA));
if isfield(S, 'problem') && isfield(S.problem, 'sD')
    assert(wrapd(S.problem.sD - sD(1)) < 1e-9, 'run_costate_library:grid', ...
           'the sheet was certified at departure phase %.6f, but .sD starts at %.6f', S.problem.sD, sD(1));
end
cols = find(isfinite(S.TF));
onlyA = d('onlyA', []);
foreignPat = d('foreignPattern', 'fine_ribs_range_job');
if ~isempty(onlyA)
    if all(onlyA == round(onlyA)) && all(onlyA >= 1) && all(onlyA <= numel(S.sA))
        want = onlyA(:).';                       % indices
    else
        want = arrayfun(@(v) find(abs(S.sA - mod(v,1)) < 1e-9, 1), onlyA(:).'); % values
    end
    excluded = setdiff(find(isfinite(S.TF)), want);
    cols = intersect(cols, want);
    fprintf('   .onlyA restricts the ribs to arrival phases %s\n', phaseList(S.sA(cols)));
    if ~isempty(excluded)
        out.blockers{end+1} = sprintf('.onlyA EXCLUDES %d certified arrival column(s) from this library: %s', ...
                                      numel(excluded), mat2str(excluded));
    end
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

%% 2. (calibration REMOVED) -- the watchdog is sized from the solver's wall
%  cap, not from a column walked outside the queue's ownership (Astra pass
%  2: an unowned calibration solve could duplicate a live worker's column)
ribOut = @(j) fullfile(outDir, sprintf('fine_rib_col%02d.mat', j));
ribList = arrayfun(ribOut, cols, 'UniformOutput', false);

%% 3. RIBS -- a work queue, one column per unit
policy = struct('staleSec', staleSec, 'maxAtt', maxAtt, 'hangSec', hangSec);
codeRoots = {here, fullfile(fileparts(fileparts(here)), 'costate_common')};
ribSpec = @(j) struct('nPts', nD - 1, 'col', j, 'sA', S.sA(j), 'nD', nD, 'problem', S.problem);
if on('ribs')
    % EXISTING OUTPUTS ARE VALIDATED BEFORE THE QUEUE CAN CALL THEM DONE.
    % A file with the right name from another grid, another problem, or a
    % save that died half-way would otherwise be "done" to the queue and
    % "not a rib" to the packager, forever (Astra pass 3, 3.3/3.5). An
    % invalid one is quarantined (renamed) under the campaign lock, and the
    % column is walked again.
    for k = 1:numel(cols)
        f = ribOut(cols(k));
        if ~isfile(f), continue, end
        [okR, whyR] = rib_validate(f, ribSpec(cols(k)));
        if ~okR
            qf = sprintf('%s.invalid.%s', f, char(java.util.UUID.randomUUID()));
            movefile(f, qf);
            out.blockers{end+1} = sprintf('column %d: existing rib file was NOT a valid rib for this unit (%s); quarantined as %s and re-queued', ...
                                          cols(k), whyR, qf);
            fprintf('3. %s\n', out.blockers{end});
        end
    end
    % 'init' OPENS an existing queue read-only (validating the unit set) and
    % creates one only when none exists -- see work_queue
    work_queue('init', out.queue, cols(:).', ribOut);
    % the exact artifact set, for anything outside MATLAB that waits on it
    fidE = fopen(fullfile(out.queue, 'expected_units.txt'), 'w');
    if fidE >= 0, fprintf(fidE, '%s\n', ribList{:}); fclose(fidE); end
    st = work_queue('status', out.queue, policy);
    fprintf('3. rib queue: %d unit(s): %d done, %d running, %d abandoned, %d retired, %d to do\n', ...
            st.n, st.nDone, st.nRunning, st.nAbandoned, st.nRetired, st.nTodo);
    if st.nRetired > 0
        out.blockers{end+1} = sprintf('%d column(s) RETIRED after %d failed attempts: %s (work_queue reset, or drop them)', ...
                                      st.nRetired, maxAtt, mat2str(st.retired));
    end
    jobFile = fullfile(outDir, 'rib_unit_job.m');
    writeRibJob(jobFile, sheetMat, outDir, out.queue, nD, codeRoots, policy);
    % THE SUPERVISOR owns the campaign from here: it keeps nWorkers alive
    % (relaunching through run_campaign_workers.sh within a budget), and
    % when every column is published runs the finalize job -- this same
    % function with the packaging stages on -- exactly once.
    finJob = fullfile(outDir, 'finalize_job.m');
    writeFinalizeJob(finJob, outDir, nD, nA, sA0, codeRoots, policy, d('extraRibFiles', {}));
    supervisor = fullfile(codeRoots{2}, 'campaign_supervisor.sh');
    out.cmd = strjoin({'nohup', shq(supervisor), shq(jobFile), num2str(nWorkers), ...
                       num2str(round(hangSec)), shq(outDir), shq(out.queue), num2str(maxAtt), shq(finJob), ...
                       '>', shq(fullfile(outDir, 'supervisor.out')), '2>&1 &'}, ' ');
    if st.nOpen > 0
        if d('launch', false)
            % workers from an older launcher hold no lock: the queue would
            % hand out a column one of them is walking (it did, 2026-09-13)
            pat = d('foreignPattern', 'fine_ribs_range_job');
            [~, fp] = system(sprintf('pgrep -f %s', shq(pat)));
            if ~isempty(strtrim(fp))
                out.state = 'blocked';
                out.blockers{end+1} = sprintf('%d foreign worker process(es) match "%s"; drain them before launching queue workers', ...
                                              numel(strsplit(strtrim(fp))), pat);
                fprintf('3. BLOCKED: %s\n', out.blockers{end});
                return
            end
            fprintf('3. starting the campaign supervisor for %d worker(s) ...\n', nWorkers);
            [rc, txt] = system(out.cmd);
            if rc ~= 0
                out.state = 'blocked';
                out.blockers{end+1} = sprintf('could not start the supervisor (%d): %s', rc, strtrim(txt));
                return
            end
            out.state = 'launched';
            fprintf(['3. supervisor started (log: %s). It keeps %d workers alive and runs the\n' ...
                     '   finalize job when every column is published; verdict in %s.\n'], ...
                    fullfile(outDir, 'supervisor.log'), nWorkers, fullfile(outDir, 'SUPERVISOR_VERDICT.txt'));
        else
            out.state = 'pending';
            fprintf(['3. %d unit(s) to walk. Start the supervised campaign with:\n\n    %s\n\n' ...
                     '   (it launches the workers, replaces any that die, and packages when done)\n'], st.nOpen, out.cmd);
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

%% 3b. COMPLETION BARRIER -- enforced whether or not stage 3 ran
% (with .run.ribs off, stage 3 is skipped, and packaging used to proceed on
% whatever happened to be on disk). A rib artifact counts only if it LOADS
% and holds a rib: a column half-way through an older, non-atomic save
% exists as a file but is not a rib.
if on('package')
    % foreign workers are excluded from PACKAGING too, not only from
    % launching: one could still be writing a column
    [~, fp] = system(sprintf('pgrep -f %s', shq(foreignPat)));
    if ~isempty(strtrim(fp))
        out.state = 'pending';
        out.blockers{end+1} = sprintf('%d foreign worker process(es) match "%s" -- nothing packaged while they run', ...
                                      numel(strsplit(strtrim(fp))), foreignPat);
        fprintf('3b. %s\n', out.blockers{end});
        return
    end
    % ONE reason per column, parallel to `good`, so the two masks below
    % index the same domain (the first version indexed a list of only the
    % short/invalid columns by a mask over all columns and threw on exactly
    % the path meant to name an invalid column -- Astra pass 3, 3.5)
    good = false(size(cols));  isShort = false(size(cols));  reason = repmat({''}, size(cols));
    for k = 1:numel(cols)
        [good(k), why, ri] = rib_validate(ribList{k}, ribSpec(cols(k)));
        if good(k) && ~ri.complete
            isShort(k) = true;
            reason{k} = sprintf('col %d: %d of %d points (%s)', cols(k), ri.nPts, nD - 1, ri.stop);
        elseif ~good(k)
            reason{k} = sprintf('col %d: %s', cols(k), why);
        end
    end
    if ~all(good)
        out.state = 'pending';
        out.blockers{end+1} = sprintf('%d of %d rib column(s) are not a publishable rib: %s -- nothing packaged', ...
                                      nnz(~good), numel(cols), strjoin(reason(~good), '; '));
        fprintf('3b. %s\n', out.blockers{end});
        return
    end
    if any(isShort)
        % a SHORT column is real data (the walker stalled where it says),
        % but it is not full coverage and the caller must see that
        out.blockers{end+1} = sprintf('%d column(s) are shorter than %d points: %s', ...
                                      nnz(isShort), nD - 1, strjoin(reason(isShort), '; '));
    end
    if isfolder(out.queue) && isfile(fullfile(out.queue, 'queue.mat'))
        stq = work_queue('status', out.queue, policy);
        if stq.nHeld > 0
            out.state = 'pending';
            out.blockers{end+1} = sprintf('%d unit lock(s) still held by live processes: %s -- nothing packaged', ...
                                          stq.nHeld, mat2str(stq.units(stq.held)));
            fprintf('3b. %s\n', out.blockers{end});
            return
        end
    end
    fprintf('3b. all %d rib columns published, none held\n', numel(cols));
end

%% 4. PACKAGE / AUDIT / SWEEP -- the existing chain, blockers not aborts
if on('package') || on('audit') || on('sweep')
    % the chain is told the GRID, the SHEET and the EXACT RIB FILES. It used
    % to be told only outDir, and its own parameter block is the 12 x 12
    % library's: default sheet name, results/arrival_rib*.mat, nD = 12.
    co = struct('outDir', outDir, 'run', struct( ...
        'arcs', false, 'sheet', false, 'ribs', false, ...
        'package', on('package'), 'audit', on('audit'), 'sweep', on('sweep'), ...
        'pictures', d('pictures', true), 'deliverable', false), ...
        'grid', struct('nA', nA, 'nD', nD, 'sA0', sA0, 'sD0', sD(1)), ...
        'sheetFile', sheetMat, 'engine', engine, 'orbits', orbits, 'invocationId', invocationId);
    extra = d('extraRibFiles', {});
    if ischar(extra), extra = {extra}; end
    extra = extra(cellfun(@isfile, extra));
    co.ribFiles = [ribList(:); extra(:)].';
    if ~isempty(extra), fprintf('4. packaging with %d extra rib file(s) from another round\n', numel(extra)); end
    % the chain sets figure defaults and closes figures; those are process-
    % global, so they are saved here, outside its clearvars, and restored
    dfv = get(0, 'DefaultFigureVisible');
    restoreFig = onCleanup(@() set(0, 'DefaultFigureVisible', dfv));
    fprintf('4. package / audit / sweep via build_70mN_library ...\n');
    % the chain is a SCRIPT, and a script runs in the workspace of whatever
    % calls it -- so it runs inside runChain's own workspace, not base. The
    % first version handed it chainOverrides through base, where the
    % script's opening clearvars wiped the caller's variables (it did, in
    % the shared session, on 2026-09-13).
    % the script also cd's; the way back is held HERE, out of its reach
    home_ = pwd;  backHome = onCleanup(@() cd(home_));
    try
        cb = runChain(co);
    catch ME
        out.state = 'failed';
        out.blockers{end+1} = sprintf('packaging chain threw: %s (%s)', ME.message, ME.identifier);
        fprintf('4. %s\n', out.blockers{end});
        return
    end
    clear backHome restoreFig
    if iscell(cb), out.blockers = [out.blockers, cb(:).']; end
end
catFile = fullfile(outDir, 'costate_catalog_dro_tulip_70mN.mat');
% 'packaged' means THIS call wrote the catalog: the chain leaves a RECEIPT
% stamped with this call's invocation id beside the catalog (an mtime
% could be met by an older file or a concurrent writer)
rcFile = fullfile(outDir, 'catalog_receipt.mat');
out.catalog = '';
if on('package')
    if isfile(rcFile)
        rc = load(rcFile);
        if strcmp(rc.invocationId, invocationId) && isfile(catFile)
            out.catalog = catFile;  out.state = 'packaged';
            out.receipt = rc;
        else
            out.blockers{end+1} = sprintf('the catalog receipt in %s is from another invocation (%s); not this call''s product', outDir, rc.invocationId);
        end
    else
        out.blockers{end+1} = sprintf('no catalog receipt was written in %s by this call', outDir);
    end
elseif on('audit') || on('sweep')
    out.state = 'measured';                        % audit/sweep ran; nothing was packaged by design
end
if ~isempty(out.blockers)
    fprintf('BLOCKERS:\n');  fprintf('  - %s\n', out.blockers{:});
end
end

% ------------------------------------------------------------------------
function chainBlockers = runChain(chainOverrides)
% RUNCHAIN  Run the chain SCRIPT in this function's private workspace, with
% chainOverrides in scope; the script assigns chainBlockers as its last act.
% Nothing else may live here: the script opens with clearvars -except
% chainOverrides, so any other variable (a default, a cleanup handle) would
% be wiped. The caller restores the working directory.
% INPUTS: chainOverrides struct.  OUTPUTS: chainBlockers cell.
build_70mN_library;
end

% ------------------------------------------------------------------------
function r = gitRevision(here)
% GITREVISION  Short git hash of the code, or 'unknown'.  INPUTS: here.
% OUTPUTS: r char.
[st, txt] = system(sprintf('cd %s && git rev-parse --short HEAD 2>/dev/null', shq(here)));
if st == 0, r = strtrim(txt); else, r = 'unknown'; end
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
function writeFinalizeJob(jobFile, outDir, nD, nA, sA0, codeRoots, policy, extraRibs)
% WRITEFINALIZEJOB  The job the supervisor runs once every column is
% published: this function again, with only the packaging stages on, on
% the same directory and grid. It exits non-zero unless it packaged.
% INPUTS: jobFile; outDir; nD; nA; sA0; codeRoots; policy; extraRibs cellstr.
% OUTPUTS: none.
roots = strjoin(cellfun(@mlq, codeRoots, 'UniformOutput', false), ', ');
if ischar(extraRibs), extraRibs = {extraRibs}; end
% DOUBLE braces: struct('f', {a, b}) builds a struct ARRAY; struct('f', {{a, b}})
% stores the cell (the round-3 finalizer failed on exactly this)
extraLit = ['{{' strjoin(cellfun(@mlq, extraRibs(:).', 'UniformOutput', false), ', ') '}}'];
txt = sprintf([ ...
 '%%%% FINALIZE_JOB  Package, audit and sweep the library. Written by run_costate_library.\n' ...
 'here = pwd; cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(here);\n' ...
 'addpath(%s);\n' ...
 'out = run_costate_library(struct(''outDir'', %s, ''nD'', %d, ''nA'', %d, ''sA0'', %.10g, ...\n' ...
 '    ''maxAtt'', %d, ''staleSec'', %d, ''hangSec'', %d, ''launch'', false, ''extraRibFiles'', %s, ...\n' ...
 '    ''run'', struct(''sheet'', false, ''ribs'', true, ''package'', true, ''audit'', true, ''sweep'', true)));\n' ...
 'fprintf(''FINALIZE: state %%s\\n'', out.state);\n' ...
 'if ~strcmp(out.state, ''packaged''), error(''finalize_job:notPackaged'', ''state %%s: %%s'', out.state, strjoin(out.blockers, '' | '')); end\n'], ...
 roots, mlq(outDir), nD, nA, sA0, policy.maxAtt, policy.staleSec, policy.hangSec, extraLit);
tmp = sprintf('%s.%s.part', jobFile, char(java.util.UUID.randomUUID()));
fid = fopen(tmp, 'w');
assert(fid >= 0, 'run_costate_library:job', 'cannot write %s', tmp);
fprintf(fid, '%s', txt);  fclose(fid);
publish_atomic(tmp, jobFile);
end

% ------------------------------------------------------------------------
function writeRibJob(jobFile, sheetMat, outDir, qDir, nD, codeRoots, policy)
% WRITERIBJOB  Emit the per-worker job script the launcher runs. It reads
% WORKER_TAG from the environment and walks whatever column the queue hands
% it. Paths are absolute and quoted as MATLAB literals; the code roots are
% the DRIVER'S own, not something inferred from the output directory (the
% first version derived a costate_common that did not exist). The worker's
% heartbeat is passed to the walker as its progress callback. The unit
% writes to the attempt-specific tmpOut the queue hands it; the WORKER
% validates that file (rib_validate) and publishes it under the lock.
% INPUTS: jobFile; sheetMat; outDir; qDir; nD; codeRoots cell; policy.
% OUTPUTS: none.
roots = strjoin(cellfun(@mlq, codeRoots, 'UniformOutput', false), ', ');
txt = sprintf([ ...
 '%%%% RIB_UNIT_JOB  One rib worker. Written by run_costate_library.\n' ...
 'here = pwd; cd(''/Users/msc/Desktop/proj7/external/pumpkynPie''); startup(); cd(here);\n' ...
 'addpath(%s);\n' ...
 'pool = capped_pool(1);\n' ...
 'assert(~isempty(pool) && isvalid(pool), ''rib_unit_job:pool'', ''no parallel pool: this worker cannot fence its solves and will not start'');\n' ...
 'tag = getenv(''WORKER_TAG'');\n' ...
 'assert(~isempty(tag), ''rib_unit_job:tag'', ''WORKER_TAG is not set: run this through run_campaign_workers.sh'');\n' ...
 'S_ = load(%s);  S_ = S_.S;\n' ...
 'ckpt = @(j) fullfile(%s, sprintf(''fine_rib_col%%02d.mat.ckpt'', j));   %% the UNIT''s checkpoint\n' ...
 'unitFcn = @(j, beat, tmpOut) build_ribs(%s, struct(''only'', j, ''direction'', -1, ...\n' ...
 '        ''nD'', %d, ''nPts'', %d, ''wallSec'', 900, ''out'', tmpOut, ''progress'', beat, ''checkpoint'', ckpt(j)));\n' ...
 'spec = @(j) struct(''nPts'', %d, ''col'', j, ''sA'', S_.sA(j), ''nD'', %d, ''problem'', S_.problem);\n' ...
 'campaign_worker(%s, %s, tag, unitFcn, struct(''logFile'', ...\n' ...
 '        fullfile(%s, [''worker_'' tag ''.log'']), ''staleSec'', %d, ''maxAtt'', %d, ...\n' ...
 '        ''validateFcn'', @(f, j) rib_validate(f, spec(j))));\n'], ...
 roots, mlq(sheetMat), mlq(outDir), mlq(sheetMat), nD, nD-1, nD-1, nD, mlq(qDir), mlq(fullfile(outDir, 'hb')), mlq(outDir), ...
 policy.staleSec, policy.maxAtt);
tmp = sprintf('%s.%s.part', jobFile, char(java.util.UUID.randomUUID()));
fid = fopen(tmp, 'w');
assert(fid >= 0, 'run_costate_library:job', 'cannot write %s', tmp);
fprintf(fid, '%s', txt);  fclose(fid);
publish_atomic(tmp, jobFile);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
