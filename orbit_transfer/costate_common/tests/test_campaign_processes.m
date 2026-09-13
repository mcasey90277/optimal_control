function ok = test_campaign_processes()
%% Purpose:
%
%   The guarantees that only mean anything across PROCESSES, tested with
%   real MATLAB worker processes launched by run_campaign_workers.sh
%   against a temporary queue (about 3 minutes; each worker is a MATLAB
%   start). Astra's pass-2 review made this a precondition for any
%   unattended run: "concurrency is the core feature being introduced, and
%   sequential tests cannot reach the races".
%
%   Units: 1,2,4,6 take 5 s; 3 takes 40 s and its owner is KILLED -9 mid-
%   unit; 5 always throws. Every attempt appends a line to the unit's
%   attempts log, so duplicate computation is visible as a count.
%
%     1. the launcher reports every worker READY (attached to the queue);
%     2. units 1,2,4,6 are computed EXACTLY ONCE across three workers (the
%        lock is exclusive across processes);
%     3. unit 3's owner is killed; another worker takes it over (the
%        kernel freed the lock) and it ends DONE; its attempts log shows
%        exactly two starts;
%     4. unit 5 is RETIRED after exactly three attempts, and a .failed
%        file is left as evidence, no .part;
%     5. the queue reaches FINISHED (all done or retired, no lock held) and
%        campaign_status reports the killed worker as not done;
%     6. no unit's artifact was written by a process without the lock: no
%        stray .part files remain.
%
%% Inputs:  none
%% Outputs: ok [logical]  (also throws when false)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
cc = fileparts(fileparts(mfilename('fullpath')));      % costate_common
addpath(cc);
T = fullfile(tempdir, ['cptest_' char(java.util.UUID.randomUUID())]);
mkdir(T);  q = fullfile(T, 'q');
% the evidence is kept when the test fails: the first run of this test
% deleted the worker logs that would have said why the workers died
keep = struct('ok', false);
cleanup = onCleanup(@() cleanupAll(T, keep));
outF = @(id) fullfile(T, sprintf('unit%02d.done', id));
work_queue('create', q, 1:6, outF);

% the unit function, as a file the workers can call
fid = fopen(fullfile(T, 'cptest_unit.m'), 'w');
fprintf(fid, ['function cptest_unit(id, beat, tmpOut, T)\n' ...
    'fid = fopen(fullfile(T, sprintf(''unit%%02d.attempts'', id)), ''a''); fprintf(fid, ''%%s start %%d\\n'', getenv(''WORKER_TAG''), matlabProcessID); fclose(fid);\n' ...
    'if id == 5, error(''cptest:always'', ''unit 5 always fails''); end\n' ...
    'dur = 5; if id == 3, dur = 40; end\n' ...
    'for k = 1:dur, pause(1); beat(); end\n' ...
    'fid = fopen(tmpOut, ''w''); fprintf(fid, ''unit %%d by %%s'', id, getenv(''WORKER_TAG'')); fclose(fid);\n' ...
    'end\n']);
fclose(fid);
job = fullfile(T, 'job.m');
fid = fopen(job, 'w');
fprintf(fid, ['addpath(''%s'', ''%s'');\n' ...
    'tag = getenv(''WORKER_TAG'');\n' ...
    'unitFcn = @(id, beat, tmpOut) cptest_unit(id, beat, tmpOut, ''%s'');\n' ...
    'campaign_worker(''%s'', ''%s'', tag, unitFcn, struct(''logFile'', fullfile(''%s'', [''worker_'' tag ''.log'']), ...\n' ...
    '    ''maxAtt'', 3, ''staleSec'', 60, ''validateFcn'', @(f) deal(isfile(f) && numel(fileread(f)) > 0, ''empty result'')));\n'], ...
    cc, T, T, q, fullfile(T, 'hb'), T);
fclose(fid);

% ---- 1. launch ---------------------------------------------------------------
launcher = fullfile(cc, 'run_campaign_workers.sh');
[rc, txt] = system(sprintf('"%s" "%s" 3 60 "%s"', launcher, job, T));
fprintf('%s', txt);
ok = chk(ok, rc == 0 && contains(txt, 'LAUNCH OK') && numel(strfind(txt, ': READY')) == 3, ...
         'the launcher reports all three workers READY (attached to the queue)');
if ~ok
    fprintf('--- launcher.log ---\n%s\n', fileread(fullfile(T, 'launcher.log')));
    w = dir(fullfile(T, 'launch_*', 'worker_*.out'));
    for m = 1:numel(w)
        fprintf('--- %s (%d bytes) ---\n%s\n', w(m).name, w(m).bytes, fileread(fullfile(w(m).folder, w(m).name)));
    end
    fprintf('--- job.m ---\n%s\n', fileread(job));
    error('test_campaign_processes:launch', 'workers did not start; evidence kept in %s', T);
end

% ---- 3. kill unit 3's owner mid-unit ----------------------------------------------
t0 = tic;  owner3 = fullfile(q, '3.owner');
while ~isfile(owner3) && toc(t0) < 120, pause(1); end
ok = chk(ok, isfile(owner3), 'unit 3 was claimed');
pause(8);                                            % well inside its 40 s
parts = strsplit(strtrim(fileread(owner3)), '|');  pid3 = str2double(parts{3});  tag3 = parts{2};
system(sprintf('kill -9 %d', pid3));
fprintf('  killed unit 3''s owner %s (pid %d) mid-unit\n', tag3, pid3);

% ---- wait for the campaign to finish --------------------------------------------
t0 = tic;
while toc(t0) < 300
    st = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
    if st.finished, break, end
    pause(5);
end
fprintf('  finished after %.0f s: %d done, %d retired, %d held\n', toc(t0), st.nDone, st.nRetired, st.nHeld);
ok = chk(ok, st.finished && st.nDone == 5 && st.nRetired == 1 && st.nHeld == 0, ...
         'the queue FINISHED: 5 done, 1 retired, no lock held');

% ---- 2. exactly-once ------------------------------------------------------------------
starts = @(id) numel(regexp(fileread(fullfile(T, sprintf('unit%02d.attempts', id))), 'start', 'match'));
n = arrayfun(starts, 1:6);
ok = chk(ok, isequal(n([1 2 4 6]), [1 1 1 1]), sprintf('units 1,2,4,6 were computed EXACTLY ONCE (starts: %s)', mat2str(n)));
ok = chk(ok, n(3) == 2 && isfile(outF(3)), 'unit 3 was taken over after its owner was killed and finished (2 starts)');
ok = chk(ok, n(5) == 3 && ismember(5, st.retired), 'unit 5 was RETIRED after exactly 3 attempts');

% ---- 4, 6. evidence and no strays ------------------------------------------------------
parts_ = dir(fullfile(T, '*.part'));  failed_ = dir(fullfile(T, '*.failed'));
ok = chk(ok, isempty(parts_), sprintf('no stray .part files (%d)', numel(parts_)));
ok = chk(ok, isempty(failed_) || all(contains({failed_.name}, 'unit05')), 'the only .failed evidence belongs to unit 5');

% ---- 5. the monitor sees the killed worker ---------------------------------------------
hbs = dir(fullfile(T, 'hb', '*.hb'));  tags = erase({hbs.name}, '.hb');
% the idle survivors notice the finished queue on their next poll (15 s),
% so give them that long before reading their final state
t0 = tic;
while toc(t0) < 90
    S = campaign_status(q, fullfile(T, 'hb'), tags, struct('staleSec', 60, 'maxAtt', 3, 'quiet', true));
    if nnz(strcmp({S.workers.state}, 'done')) >= 2, break, end
    pause(3);
end
k3 = strcmp({S.workers.tag}, tag3);
ok = chk(ok, any(k3) && ~strcmp(S.workers(k3).state, 'done'), ...
         sprintf('the killed worker %s is reported %s, not done', tag3, S.workers(k3).state));
ok = chk(ok, nnz(strcmp({S.workers.state}, 'done')) == 2, 'the two surviving workers exited cleanly (done)');

if ok, fprintf('TEST_CAMPAIGN_PROCESSES: ALL PASS\n');  keep.ok = true;  cleanup = onCleanup(@() cleanupAll(T, keep));
else,  fprintf('TEST_CAMPAIGN_PROCESSES: FAIL (evidence kept in %s)\n', T);  error('test_campaign_processes:fail', 'FAILED'); end
end

function cleanupAll(T, keep)
% CLEANUPALL  Kill any worker still running against T; remove T only after
% a pass.  INPUTS: T; keep struct (.ok).  OUTPUTS: none.
% the job path is in the workers' ENVIRONMENT, not their argv, so match on
% the tag prefix the launcher gives them instead
system(sprintf('pkill -9 -f "CAMPAIGN_JOB=%s" 2>/dev/null; pkill -9 -f "%s" 2>/dev/null', fullfile(T, 'job.m'), fullfile(T, 'job.m')));
pause(1);
if keep.ok, try rmdir(T, 's'); catch, end, end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
