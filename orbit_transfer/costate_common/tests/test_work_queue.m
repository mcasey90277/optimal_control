function ok = test_work_queue()
%% Purpose:
%
%   Tests the campaign primitives against the SPECIFIC failures that cost a
%   day on the 24x24 library. Each check names the incident it exists for,
%   so a future change that reintroduces one fails here with the reason.
%
%     1. two workers never get the same unit (atomic claim);
%     2. a unit is DONE when its artifact exists -- never a separate flag
%        that can disagree with the results on disk;
%     3. a worker that DIES does not take its unit out of circulation: a
%        stale claim is reclaimed. (B and C were killed by their watchdogs
%        holding columns 13 and 16; nothing could pick them up, and nine
%        hours of walking was lost.)
%     4. the heartbeat tells NEVER STARTED apart from RUNNING. (Workers E,
%        F and G never launched because a shell loop globbed a bracket; the
%        monitor read the missing log as "alive, age 0 s".)
%     5. STALLED is distinguished from running and from done;
%     6. a reporting block cannot fail the work: safe_report catches the
%        exact expression that reported a completed, saved sheet as a failed
%        job, and says it is a bug in the report;
%     7. fmt_num renders the missing cases at the SAME width rather than
%        throwing -- NaN t_f and an empty count, both ordinary states of a
%        sparse grid.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
q = fullfile(tempdir, ['wqtest_' char(datetime('now','Format','yyyyMMdd_HHmmssSSS'))]);
mkdir(q);
cleanup = onCleanup(@() rmdir(q, 's'));
outF = @(id) fullfile(q, sprintf('unit%02d.done', id));
work_queue('init', q, 1:6, outF);

% ---- 1. atomic claiming -------------------------------------------------
a = work_queue('claim', q, 'A');  b = work_queue('claim', q, 'B');
ok = chk(ok, a.id ~= b.id && ~isnan(a.id) && ~isnan(b.id), ...
         sprintf('two workers claim different units (%d and %d)', a.id, b.id));

% ---- 2. done is the ARTIFACT --------------------------------------------
touchFile(outF(a.id));  work_queue('release', q, a.claim);
a2 = work_queue('claim', q, 'A');
ok = chk(ok, a2.id ~= a.id, sprintf('a finished unit is not handed out again (got %d, not %d)', a2.id, a.id));
s = work_queue('status', q);
ok = chk(ok, s.nDone == 1 && ismember(a.id, s.done), 'status reads DONE from the artifact, not a flag');

% ---- 3. a dead worker's unit returns to the queue ------------------------
ageFile(fullfile(b.claim, 'beat'));
c = work_queue('claim', q, 'C', struct('staleSec', 60));
ok = chk(ok, c.id == b.id, ...
         sprintf('a stale claim is reclaimed: C took unit %d back from the dead worker B', c.id));
s2 = work_queue('status', q);
ok = chk(ok, s2.nRunning >= 1, 'and it counts as running again');

% ---- 4, 5. the five heartbeat states ------------------------------------
hb = fullfile(q, 'hb');
campaign_heartbeat('beat', hb, 'W1', 'working');
campaign_heartbeat('done', hb, 'W2', 'finished');
campaign_heartbeat('fail', hb, 'W3', 'threw');
campaign_heartbeat('beat', hb, 'W4', 'old');   ageFile(fullfile(hb, 'W4.hb'));
r = campaign_heartbeat('read', hb, {'W1','W2','W3','W4','W5'}, struct('staleSec', 60));
want = {'running','done','failed','stalled','never'};
got = {r.state};
ok = chk(ok, isequal(got, want), sprintf('five states told apart: %s', strjoin(got, ' ')));
ok = chk(ok, strcmp(r(5).state, 'never') && contains(r(5).msg, 'never started'), ...
         'a worker that NEVER STARTED is not reported as alive');

% ---- 6. a report cannot fail the work -----------------------------------
okr = safe_report(@() char(string(NaN)), 'the fine-sheet summary');
ok = chk(ok, ~okr, 'safe_report catches the print that reported a saved sheet as failed');
okr2 = safe_report(@() fprintf(''), 'a clean report');
ok = chk(ok, okr2, 'and returns true when the report is fine');

% ---- 7. fixed-width missing values ---------------------------------------
w = 8;
a1 = fmt_num(NaN, w, 4);  a2 = fmt_num([], w, 4);  a3 = fmt_num(17.7976, w, 4);
ok = chk(ok, numel(a1) == w && numel(a2) == w && numel(a3) == w, ...
         sprintf('NaN, empty and a value all render at width %d ("%s" "%s" "%s")', w, a1, a2, a3));
ok = chk(ok, all(a1 == '-'), 'the missing case is dashes, not an error');

if ok, fprintf('TEST_WORK_QUEUE: ALL PASS\n'); else, fprintf('TEST_WORK_QUEUE: FAIL\n'); end
end

function touchFile(f)
% TOUCHFILE  Create an empty file.  INPUTS: f.  OUTPUTS: none.
fid = fopen(f, 'w');  fprintf(fid, 'x');  fclose(fid);
end

function ageFile(f)
% AGEFILE  Backdate a file so it reads as stale.  INPUTS: f.  OUTPUTS: none.
system(sprintf('touch -t 202001010000 "%s"', f));
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
