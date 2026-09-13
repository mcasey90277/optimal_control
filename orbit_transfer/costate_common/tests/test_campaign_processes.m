function ok = test_campaign_processes()
%% Purpose:
%
%   The guarantees that only mean anything across PROCESSES, tested with
%   real MATLAB worker processes launched by run_campaign_workers.sh
%   against temporary queues (about six minutes). Astra's pass-2 review
%   made this a precondition for any unattended run; its pass-3 review
%   named the smallest cases that close the most risk, and they are here.
%
%   PHASE 0  registry across processes: this process holds a lock, clears
%            its functions, probes its own lock (the sequence that used to
%            drop the kernel lock); a second MATLAB process must still be
%            refused.
%   PHASE 1  six units, three workers. 1,2,4,6 take 5 s. 3 takes 40 s and
%            writes its temporary result after 2 s; its owner is killed -9
%            after that write, so the takeover must not publish the dead
%            owner's file. 5 writes a partial result then throws, every
%            time. Checks: all READY; 1,2,4,6 computed exactly once; 3
%            taken over and published by the REPLACEMENT; 5 retired after
%            exactly 3 attempts with .failed evidence; the only .part
%            debris belongs to the killed attempt; queue finished with no
%            lock held; exit records: survivors 0, killed one non-zero.
%   PHASE 2  one unit, maxAtt = 1, its owner killed mid-unit: the unit is
%            RETIRED (terminal), nobody starts it again, queue finished.
%   PHASE 3  one unit whose worker beats once then goes silent; hangSec =
%            60: the launcher's per-worker supervisor must kill it, record
%            the reason and a non-zero exit, and the unit must be ABANDONED
%            (claimable), with the monitor alarming that no live worker
%            remains.
%   PHASE 4  four units under campaign_supervisor with two workers; one
%            worker is killed -9: the supervisor must launch a replacement,
%            see every unit published, run the finalizer exactly once, and
%            release its lock.
%
%% Inputs:  none
%% Outputs: ok [logical]  (also throws when false; evidence kept on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
cc = fileparts(fileparts(mfilename('fullpath')));      % costate_common
addpath(cc);
T = fullfile(tempdir, ['cptest_' char(java.util.UUID.randomUUID())]);
mkdir(T);
keep = struct('ok', false);
cleanup = onCleanup(@() cleanupAll(T, keep));
launcher = fullfile(cc, 'run_campaign_workers.sh');
M = '/Applications/MATLAB_R2026a.app/bin/matlab';

% ---- PHASE 0: the registry across processes ------------------------------------
lf = fullfile(T, 'reg.lock');
L = unit_lock('try', lf);
clear functions                                        % wipes ordinary persistent state
pr = unit_lock('probe', lf);                           % the sequence that used to free the lock
probeOut = fullfile(T, 'other_got_it.txt');
fid = fopen(fullfile(T, 'probe_other.m'), 'w');
fprintf(fid, ['addpath(''%s''); Tt = unit_lock(''try'', ''%s''); fid = fopen(''%s'', ''w''); ' ...
              'fprintf(fid, ''%%d'', Tt.held); fclose(fid); if Tt.held, unit_lock(''release'', ''%s'', Tt.token); end\n'], ...
        cc, lf, probeOut, lf);
fclose(fid);
system(sprintf('%s -batch "run(''%s'')" > %s 2>&1', M, fullfile(T, 'probe_other.m'), fullfile(T, 'probe_other.log')));
ok = chk(ok, L.held && pr.mine && unit_lock('holds', lf, L.token) && strcmp(strtrim(fileread(probeOut)), '0'), ...
         'PHASE 0: after clear functions + own probe, a second PROCESS is still refused the lock');
unit_lock('release', lf, L.token);

% ---- the unit function, as a file the workers can call ----------------------------
fid = fopen(fullfile(T, 'cptest_unit.m'), 'w');
fprintf(fid, ['function cptest_unit(id, beat, tmpOut, T)\n' ...
    'fid = fopen(fullfile(T, sprintf(''unit%%02d.attempts'', id)), ''a''); fprintf(fid, ''%%s start %%d\\n'', getenv(''WORKER_TAG''), matlabProcessID); fclose(fid);\n' ...
    'if id == 5, fid = fopen(tmpOut, ''w''); fprintf(fid, ''partial''); fclose(fid); error(''cptest:always'', ''unit 5 always fails''); end\n' ...
    'if id == 8, beat(); pause(300); return, end\n' ...          % PHASE 3: silent after one beat
    'dur = 5; if id == 3 || id == 7, dur = 40; end; if id >= 9, dur = 12; end\n' ...
    'for k = 1:dur\n' ...
    '    pause(1); beat();\n' ...
    '    if k == 2, fid = fopen(tmpOut, ''w''); fprintf(fid, ''unit %%d by %%s (early)'', id, getenv(''WORKER_TAG'')); fclose(fid); end\n' ...
    'end\n' ...
    'fid = fopen(tmpOut, ''w''); fprintf(fid, ''unit %%d by %%s'', id, getenv(''WORKER_TAG'')); fclose(fid);\n' ...
    'end\n']);
fclose(fid);
writeJob = @(name, q, maxAtt) writeJobFile(fullfile(T, name), cc, T, q, maxAtt);

% ---- PHASE 1 -------------------------------------------------------------------------
q = fullfile(T, 'q');  outF = @(id) fullfile(T, sprintf('unit%02d.done', id));
work_queue('create', q, 1:6, outF);
job = writeJob('job.m', q, 3);
[rc, txt] = system(sprintf('"%s" "%s" 3 60 "%s"', launcher, job, T));
fprintf('%s', txt);
ok = chk(ok, rc == 0 && contains(txt, 'LAUNCH OK') && numel(strfind(txt, ': READY')) == 3, ...
         'PHASE 1: the launcher reports all three workers READY');
if ~ok, dumpEvidence(T, job); error('test_campaign_processes:launch', 'workers did not start; evidence kept in %s', T); end

t0 = tic;  owner3 = fullfile(q, '3.owner');
while ~isfile(owner3) && toc(t0) < 120, pause(1); end
ok = chk(ok, isfile(owner3), 'unit 3 was claimed');
parts = strsplit(strtrim(fileread(owner3)), '|');  tag3 = parts{2};  pid3 = str2double(parts{3});
tok3 = extractAfter(parts{1}, ':');
% kill AFTER its early temporary write (the write is at 2 s; wait for the file)
t0 = tic;  early3 = fullfile(T, sprintf('unit03.done.%s.part', tok3));
while ~isfile(early3) && toc(t0) < 60, pause(0.5); end
ok = chk(ok, isfile(early3), 'unit 3''s owner wrote its temporary result');
system(sprintf('kill -9 %d', pid3));
fprintf('  killed unit 3''s owner %s (pid %d) after its temporary write\n', tag3, pid3);

t0 = tic;
while toc(t0) < 300
    st = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
    if st.finished, break, end
    pause(5);
end
fprintf('  finished after %.0f s: %d done, %d retired, %d held\n', toc(t0), st.nDone, st.nRetired, st.nHeld);
ok = chk(ok, st.finished && st.nDone == 5 && st.nRetired == 1 && st.nHeld == 0, ...
         'the queue FINISHED: 5 done, 1 retired, no lock held');
starts = @(id) numel(regexp(fileread(fullfile(T, sprintf('unit%02d.attempts', id))), 'start', 'match'));
n = arrayfun(starts, 1:6);
ok = chk(ok, isequal(n([1 2 4 6]), [1 1 1 1]), sprintf('units 1,2,4,6 were computed EXACTLY ONCE (starts: %s)', mat2str(n)));
ok = chk(ok, n(3) == 2 && isfile(outF(3)) && ~contains(fileread(outF(3)), tag3) && ~contains(fileread(outF(3)), 'early'), ...
         'unit 3 was taken over and PUBLISHED BY THE REPLACEMENT, not from the dead owner''s file');
ok = chk(ok, n(5) == 3 && ismember(5, st.retired), 'unit 5 was RETIRED after exactly 3 attempts');
failed_ = dir(fullfile(T, 'unit05.done.*.failed'));
ok = chk(ok, numel(failed_) == 3 && all(arrayfun(@(f) strcmp(fileread(fullfile(f.folder, f.name)), 'partial'), failed_)), ...
         'each of unit 5''s three attempts left its partial result as .failed evidence');
parts_ = dir(fullfile(T, '*.part'));
ok = chk(ok, all(contains({parts_.name}, tok3)), sprintf('the only .part debris belongs to the killed attempt (%d file(s))', numel(parts_)));
ld = dir(fullfile(T, 'launch_*'));  ld = ld([ld.isdir]);  ldir = fullfile(ld(end).folder, ld(end).name);
t0 = tic;  while numel(dir(fullfile(ldir, 'exit_*'))) < 3 && toc(t0) < 120, pause(2); end
ex = dir(fullfile(ldir, 'exit_*'));
rcs = nan(1, numel(ex));  reasons = cell(1, numel(ex));
for k = 1:numel(ex)
    t = strsplit(strtrim(fileread(fullfile(ex(k).folder, ex(k).name))), ' ');
    rcs(k) = str2double(t{1});  reasons{k} = strjoin(t(2:end), ' ');
end
killedRow = contains({ex.name}, tag3);
ok = chk(ok, numel(ex) == 3 && all(rcs(~killedRow) == 0) && rcs(killedRow) ~= 0, ...
         sprintf('exit records: survivors rc 0, killed worker rc %d (reaped by its supervisor)', rcs(killedRow)));
hbs = dir(fullfile(T, 'hb', '*.hb'));  tags = erase({hbs.name}, '.hb');
S = campaign_status(q, fullfile(T, 'hb'), tags, struct('staleSec', 60, 'maxAtt', 3, 'quiet', true));
k3 = strcmp({S.workers.tag}, tag3);
ok = chk(ok, any(k3) && ~strcmp(S.workers(k3).state, 'done') && nnz(strcmp({S.workers.state}, 'done')) == 2, ...
         'the monitor reports the killed worker not done and the survivors done');

% ---- PHASE 2: a kill on the FINAL attempt retires the unit ---------------------------
q2 = fullfile(T, 'q2');  outF2 = @(id) fullfile(T, sprintf('unit%02d.done', id));
work_queue('create', q2, 7, outF2);
job2 = writeJob('job2.m', q2, 1);
[rc, txt] = system(sprintf('"%s" "%s" 1 60 "%s"', launcher, job2, T));
ok = chk(ok, rc == 0 && contains(txt, 'LAUNCH OK'), 'PHASE 2: one worker READY');
owner7 = fullfile(q2, '7.owner');  t0 = tic;
while ~isfile(owner7) && toc(t0) < 120, pause(1); end
parts = strsplit(strtrim(fileread(owner7)), '|');  pid7 = str2double(parts{3});
pause(5);  system(sprintf('kill -9 %d', pid7));
t0 = tic;
while toc(t0) < 120
    st2 = work_queue('status', q2, struct('maxAtt', 1, 'staleSec', 60));
    if ~st2.held(1), break, end
    pause(2);
end
ok = chk(ok, ismember(7, st2.retired) && st2.finished && st2.nAbandoned == 0, ...
         'a unit whose owner died on its LAST attempt is RETIRED and the queue is finished (no limbo)');
ok = chk(ok, starts(7) == 1, 'and nobody started it again');

% ---- PHASE 3: the supervisor kills a silent worker -----------------------------------
q3 = fullfile(T, 'q3');  outF3 = @(id) fullfile(T, sprintf('unit%02d.done', id));
work_queue('create', q3, 8, outF3);
job3 = writeJob('job3.m', q3, 3);
[rc, txt] = system(sprintf('"%s" "%s" 1 60 "%s"', launcher, job3, T));
ok = chk(ok, rc == 0 && contains(txt, 'LAUNCH OK'), 'PHASE 3: one worker READY');
ld = dir(fullfile(T, 'launch_*'));  ld = ld([ld.isdir]);  ldir3 = fullfile(ld(end).folder, ld(end).name);
t0 = tic;
while isempty(dir(fullfile(ldir3, 'exit_*'))) && toc(t0) < 200, pause(3); end
ex3 = dir(fullfile(ldir3, 'exit_*'));
if isempty(ex3), rec = 'NO EXIT RECORD'; else, rec = strtrim(fileread(fullfile(ex3(1).folder, ex3(1).name))); end
ok = chk(ok, ~isempty(ex3) && contains(rec, 'no heartbeat') && ~startsWith(rec, '0 '), ...
         sprintf('the SUPERVISOR killed the silent worker after ~%.0f s and recorded it: "%s"', toc(t0), rec));
st3 = work_queue('status', q3, struct('maxAtt', 3, 'staleSec', 60));
ok = chk(ok, ismember(8, st3.abandoned) && st3.nHeld == 0, 'its unit is ABANDONED (lock freed by the kill) and claimable');
hbs3 = dir(fullfile(T, 'hb', 'w1_*.hb'));  tags3 = erase({hbs3.name}, '.hb');
S3 = campaign_status(q3, fullfile(T, 'hb'), tags3(end), struct('staleSec', 60, 'maxAtt', 3, 'quiet', true));
ok = chk(ok, any(contains(S3.alarm, 'NO LIVE WORKER')), 'and the monitor alarms that no live worker remains');
% ---- PHASE 4: the SUPERVISOR replaces a killed worker and finalizes once ------------
q4 = fullfile(T, 'q4');  outF4 = @(id) fullfile(T, sprintf('unit%02d.done', id));
work_queue('create', q4, 9:12, outF4);                     % four units of 5 s (ids 9..12)
fidE = fopen(fullfile(q4, 'expected_units.txt'), 'w'); fprintf(fidE, '%s\n', outF4(9), outF4(10), outF4(11), outF4(12)); fclose(fidE);
job4 = writeJob('job4.m', q4, 3);
fin = fullfile(T, 'finalize.m');
fid = fopen(fin, 'w'); fprintf(fid, 'fid = fopen(''%s'', ''a''); fprintf(fid, ''finalized %%s\\n'', char(datetime(''now''))); fclose(fid);\n', fullfile(T, 'FINALIZED')); fclose(fid);
sup = fullfile(cc, 'campaign_supervisor.sh');
system(sprintf('POLL_SEC=5 nohup "%s" "%s" 2 60 "%s" "%s" 3 "%s" > "%s" 2>&1 &', sup, job4, T, q4, fin, fullfile(T, 'supervisor.out')));
% wait for the first launch to be READY, then kill one worker
t0 = tic;
while numel(dir(fullfile(T, 'launch_*', 'pid_*'))) < 5 && toc(t0) < 120, pause(1); end   % 3+1+1 earlier, +2 now = 7
t0 = tic;  while numel(dir(fullfile(T, 'launch_*', 'pid_*'))) < 7 && toc(t0) < 120, pause(1); end
pids4 = dir(fullfile(T, 'launch_*', 'pid_*'));  [~, order] = sort([pids4.datenum]);  pids4 = pids4(order(end-1:end));
pause(8);
victim = str2double(strtrim(fileread(fullfile(pids4(1).folder, pids4(1).name))));
system(sprintf('kill -9 %d', victim));
fprintf('  killed supervised worker pid %d\n', victim);
t0 = tic;
while ~isfile(fullfile(T, 'SUPERVISOR_VERDICT.txt')) && toc(t0) < 240, pause(3); end
sv = ''; if isfile(fullfile(T, 'SUPERVISOR_VERDICT.txt')), sv = fileread(fullfile(T, 'SUPERVISOR_VERDICT.txt')); end
nLaunch = numel(dir(fullfile(T, 'launch_*', 'pid_*')));
ok = chk(ok, contains(sv, 'FINISHED') && contains(sv, 'finalizer exited 0'), ...
         sprintf('PHASE 4: the supervisor saw the campaign through and ran the finalizer (%s)', strtrim(regexprep(sv, '=== [^\n]* ===\s*', ''))));
ok = chk(ok, nLaunch >= 8, sprintf('it launched a REPLACEMENT after the kill (%d worker launches in all)', nLaunch));
ok = chk(ok, all(arrayfun(@(id) isfile(outF4(id)), 9:12)), 'all four units were published');
ok = chk(ok, isfile(fullfile(T, 'FINALIZED')) && numel(regexp(fileread(fullfile(T, 'FINALIZED')), 'finalized', 'match')) == 1, ...
         'the finalizer ran exactly once');
ok = chk(ok, ~isfolder(fullfile(T, 'supervisor.lock')), 'the supervisor released its lock on exit');

pidsAll = dir(fullfile(T, 'launch_*', 'pid_*'));
alive = 0;
for k = 1:numel(pidsAll)
    p = str2double(strtrim(fileread(fullfile(pidsAll(k).folder, pidsAll(k).name))));
    [s_, ~] = system(sprintf('kill -0 %d 2>/dev/null', p));  alive = alive + (s_ == 0);
end
ok = chk(ok, alive == 0, 'no launched worker process is still alive');

if ok, fprintf('TEST_CAMPAIGN_PROCESSES: ALL PASS\n');  keep.ok = true;  cleanup = onCleanup(@() cleanupAll(T, keep));
else,  fprintf('TEST_CAMPAIGN_PROCESSES: FAIL (evidence kept in %s)\n', T);  error('test_campaign_processes:fail', 'FAILED'); end
end

function job = writeJobFile(job, cc, T, q, maxAtt)
% WRITEJOBFILE  A worker job for queue q.  INPUTS: job path; cc; T; q;
% maxAtt.  OUTPUTS: job path.
fid = fopen(job, 'w');
fprintf(fid, ['addpath(''%s'', ''%s'');\n' ...
    'tag = getenv(''WORKER_TAG'');\n' ...
    'unitFcn = @(id, beat, tmpOut) cptest_unit(id, beat, tmpOut, ''%s'');\n' ...
    'campaign_worker(''%s'', ''%s'', tag, unitFcn, struct(''logFile'', fullfile(''%s'', [''worker_'' tag ''.log'']), ...\n' ...
    '    ''maxAtt'', %d, ''staleSec'', 60, ''validateFcn'', @(f, id) deal(isfile(f) && numel(fileread(f)) > 0, ''empty result'')));\n'], ...
    cc, T, T, q, fullfile(T, 'hb'), T, maxAtt);
fclose(fid);
end

function dumpEvidence(T, job)
% DUMPEVIDENCE  Print what the launcher and workers left.  INPUTS: T; job.
% OUTPUTS: none.
if isfile(fullfile(T, 'launcher.log')), fprintf('--- launcher.log ---\n%s\n', fileread(fullfile(T, 'launcher.log'))); end
w = dir(fullfile(T, 'launch_*', 'worker_*.out'));
for m = 1:numel(w)
    fprintf('--- %s (%d bytes) ---\n%s\n', w(m).name, w(m).bytes, fileread(fullfile(w(m).folder, w(m).name)));
end
fprintf('--- job ---\n%s\n', fileread(job));
end

function cleanupAll(T, keep)
% CLEANUPALL  Kill every worker this test launched, BY THE PIDS THE LAUNCHER
% RECORDED (the job path is in their environment, invisible to pgrep -f),
% and their descendants; remove T only after a pass.
% INPUTS: T; keep struct (.ok).  OUTPUTS: none.
pids = dir(fullfile(T, 'launch_*', 'pid_*'));
for k = 1:numel(pids)
    p = strtrim(fileread(fullfile(pids(k).folder, pids(k).name)));
    system(sprintf('pkill -9 -P %s 2>/dev/null; kill -9 %s 2>/dev/null', p, p));
end
pause(1);
if keep.ok, try rmdir(T, 's'); catch, end, end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
