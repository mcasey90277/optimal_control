function ok = test_work_queue()
%% Purpose:
%
%   Tests the campaign primitives against the SPECIFIC failures that cost a
%   day on the 24x24 library, and the ones the Astra chain review of
%   2026-09-13 showed the first version still had. Each check names the
%   incident it exists for, so a change that reintroduces one fails here
%   with the reason. THIS FUNCTION THROWS on failure, so a batch job that
%   calls it cannot exit clean with a failed suite.
%
%     1. two workers never get the same unit (atomic claim);
%     2. a unit is DONE when its artifact exists -- never a separate flag;
%     3. a worker that DIES does not take its unit out of circulation: a
%        stale claim is reclaimed (B and C were killed holding 13 and 16);
%     4. the heartbeat tells NEVER STARTED apart from RUNNING (E, F, G
%        never launched; the monitor read the missing log as alive);
%     5. STALLED is distinguished from running and from done, and an
%        unreadable heartbeat is UNKNOWN, not running;
%     6. a report cannot fail the work (safe_report);
%     7. fmt_num renders missing values at the same width;
%     8. RE-OPENING the queue keeps live claims AND attempt counts (the
%        first version wiped both on every init, so "re-run to package"
%        would have freed columns still being walked and re-armed the
%        livelock the counter was built to stop);
%     9. a unit that fails maxAtt times is RETIRED and stays retired across
%        a re-open; 'reset' is the only way back, and it refuses while a
%        claim is fresh;
%    10. a superseded owner can neither beat nor release the new owner's
%        claim (ownership is a token, not a path);
%    11. status gives ONE state per unit: a live final attempt is running,
%        not retired, and the counts partition the units;
%    12. a claim won after another worker finished the unit is given back
%        (re-check under ownership).
%
%% Inputs:  none
%% Outputs: ok [logical]  (also throws when false)
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
ok = chk(ok, ~isempty(a.token) && ~strcmp(a.token, b.token), 'each claim carries its own owner token');

% ---- 2. done is the ARTIFACT --------------------------------------------
touchFile(outF(a.id));  rel = work_queue('release', q, a.claim, a.token);
ok = chk(ok, rel, 'the owner can release its claim');
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
ok = chk(ok, ismember(c.id, s2.running), 'and it counts as running again');

% ---- 10. the dead worker comes back: it must not touch C's claim ---------
okB = work_queue('beat', q, b.claim, b.token);
ok = chk(ok, ~okB, 'a superseded owner cannot refresh the new owner''s claim');
okR = work_queue('release', q, b.claim, b.token);
ok = chk(ok, ~okR && isfolder(c.claim), 'nor release it');
okC = work_queue('beat', q, c.claim, c.token);
ok = chk(ok, okC, 'while the live owner still can');

% ---- 8. re-open keeps claims and attempts --------------------------------
attBefore = work_queue('status', q).attempts;
work_queue('init', q, 1:6, outF);                      % what a re-run does
s3 = work_queue('status', q);
ok = chk(ok, isfolder(c.claim) && ismember(c.id, s3.running), ...
         're-opening the queue keeps a LIVE claim (a re-run to package must not free a column being walked)');
ok = chk(ok, isequal(s3.attempts, attBefore) && any(attBefore > 0), ...
         'and keeps the attempt counts (wiping them re-armed the livelock)');
o = work_queue('open', q, 1:8, outF);
ok = chk(ok, isequal(o.added, [7 8]) && numel(work_queue('status', q).units) == 8, ...
         'open adds new units without disturbing the old ones');
try
    work_queue('create', q, 1:8, outF);  created = true;
catch
    created = false;
end
ok = chk(ok, ~created, 'create refuses to overwrite an existing queue');

% ---- 9. retirement persists; reset is explicit ---------------------------
% burn unit a2 through its attempts: claim, release (failed), repeat
id9 = a2.id;  work_queue('release', q, a2.claim, a2.token);
for k = 1:3
    r = work_queue('claim', q, 'F', struct('maxAtt', 3));
    work_queue('release', q, r.claim, r.token);
    if r.id ~= id9, break, end
end
r = work_queue('claim', q, 'F', struct('maxAtt', 3));
ok = chk(ok, r.id ~= id9, sprintf('after 3 failed attempts unit %d is RETIRED and not handed out', id9));
if ~isnan(r.id), work_queue('release', q, r.claim, r.token); end
work_queue('init', q, 1:8, outF);
s9 = work_queue('status', q, struct('maxAtt', 3));
ok = chk(ok, ismember(id9, s9.retired), 'and it STAYS retired across a re-open');
try
    work_queue('reset', q, struct('staleSec', 60));  didReset = true;
catch
    didReset = false;
end
ok = chk(ok, ~didReset, 'reset refuses while a claim is fresh (a worker is on it)');
% corrupt attempt record fails CLOSED
corrupt = fullfile(q, sprintf('%d.att', 7));  fid = fopen(corrupt, 'w');  fprintf(fid, 'garbage');  fclose(fid);
s9b = work_queue('status', q, struct('maxAtt', 3));
ok = chk(ok, ismember(7, s9b.retired), 'an unreadable attempt record blocks the unit rather than re-opening it');
delete(corrupt);

% ---- 11. one state per unit -----------------------------------------------
s11 = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
ok = chk(ok, s11.nDone + s11.nRunning + s11.nStale + s11.nRetired + s11.nTodo == s11.n, ...
         sprintf('the five counts PARTITION the %d units (%d/%d/%d/%d/%d)', s11.n, ...
                 s11.nDone, s11.nRunning, s11.nStale, s11.nRetired, s11.nTodo));
% a live FINAL attempt is running, not retired: take a unit nobody has
% touched, so the count below is the whole of its history
id11 = s11.todo(find(s11.attempts(ismember(s11.units, s11.todo)) == 0, 1));
for k = 1:2
    r = claimSpecific(q, 'G', id11, 3);  work_queue('release', q, r.claim, r.token);
end
r3 = claimSpecific(q, 'G', id11, 3);
s11b = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
ok = chk(ok, ismember(id11, s11b.running) && ~ismember(id11, s11b.retired), ...
         sprintf('unit %d on its last attempt is RUNNING, not retired (was counted as both)', id11));
work_queue('release', q, r3.claim, r3.token);

% ---- 12. re-check under ownership ----------------------------------------
% a unit finished by someone else between the pre-check and the mkdir is
% given straight back: simulate by finishing it, then claiming
id12 = s11b.todo(1);  touchFile(outF(id12));
r12 = work_queue('claim', q, 'H', struct('maxAtt', 3));
ok = chk(ok, r12.id ~= id12 && ~isfolder(fullfile(q, sprintf('%d.claim', id12))), ...
         'a unit that became done is not claimed and no claim is left behind');
if ~isnan(r12.id), work_queue('release', q, r12.claim, r12.token); end

% ---- 4, 5. the five heartbeat states, plus UNKNOWN ------------------------
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
fid = fopen(fullfile(hb, 'W6.hb'), 'w');  fclose(fid);          % truncated record
r6 = campaign_heartbeat('read', hb, 'W6', struct('staleSec', 60));
ok = chk(ok, strcmp(r6.state, 'unknown'), 'an empty or garbled heartbeat is UNKNOWN, not running');

% ---- 6. a report cannot fail the work -----------------------------------
okr = safe_report(@() char(string(NaN)), 'the fine-sheet summary');
ok = chk(ok, ~okr, 'safe_report catches the print that reported a saved sheet as failed');
okr2 = safe_report(@() fprintf(''), 'a clean report');
ok = chk(ok, okr2, 'and returns true when the report is fine');

% ---- 7. fixed-width missing values ---------------------------------------
w = 8;
a1 = fmt_num(NaN, w, 4);  a2v = fmt_num([], w, 4);  a3 = fmt_num(17.7976, w, 4);
ok = chk(ok, numel(a1) == w && numel(a2v) == w && numel(a3) == w, ...
         sprintf('NaN, empty and a value all render at width %d ("%s" "%s" "%s")', w, a1, a2v, a3));
ok = chk(ok, all(a1 == '-'), 'the missing case is dashes, not an error');

if ok, fprintf('TEST_WORK_QUEUE: ALL PASS\n');
else,  fprintf('TEST_WORK_QUEUE: FAIL\n');  error('test_work_queue:fail', 'test_work_queue FAILED'); end
end

function r = claimSpecific(q, tag, id, maxAtt)
% CLAIMSPECIFIC  Claim until the queue hands out unit id, releasing others.
% INPUTS: q; tag; id; maxAtt.  OUTPUTS: r claim struct.
held = {};
while true
    r = work_queue('claim', q, tag, struct('maxAtt', maxAtt));
    assert(~isnan(r.id), 'claimSpecific: unit %d never came up', id);
    if r.id == id, break, end
    held{end+1} = r; %#ok<AGROW>
end
for k = 1:numel(held), work_queue('release', q, held{k}.claim, held{k}.token); end
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
