function ok = test_work_queue()
%% Purpose:
%
%   Tests the campaign primitives, one check per guarantee, each named for
%   the incident it exists for. THROWS on failure. These run in ONE MATLAB
%   process; the guarantees that only mean anything across PROCESSES (a
%   live owner cannot be stolen from, a killed owner's unit is taken over,
%   a final-attempt kill retires) are tested by test_campaign_processes,
%   which launches real workers.
%
%     1. two claims never get the same unit; each carries its own token;
%     2. DONE is the published artifact, published only through the queue,
%        only by the owner, only after validation;
%     3. a unit whose owner released without publishing goes back;
%     4. attempts persist across 'open'; 'open' is read-only and refuses a
%        different unit set (the first version wiped claims and attempts
%        on every init, re-arming the livelock);
%     5. maxAtt failed attempts RETIRE a unit, across re-opens; 'reset'
%        clears only unheld units; an unreadable attempt record blocks;
%     6. status gives ONE state per unit and the counts partition;
%        'finished' needs no held lock;
%     7. the heartbeat tells never / running / stalled / failed / done
%        apart and calls an empty, truncated or foreign record UNKNOWN;
%     8. a report cannot fail the work (safe_report), and fmt_num renders
%        missing and over-wide values at the same width;
%     9. publish_atomic replaces in one rename and refuses a directory;
%    10. unit_lock: a held lock is reported held by a probe, released on
%        release.
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
q = fullfile(tempdir, ['wqtest_' char(java.util.UUID.randomUUID())]);
mkdir(q);
cleanup = onCleanup(@() rmdir(q, 's'));
outF = @(id) fullfile(q, sprintf('unit%02d.done', id));
work_queue('init', q, 1:6, outF);

% ---- 1. claiming ---------------------------------------------------------
a = work_queue('claim', q, 'A');  b = work_queue('claim', q, 'B');
ok = chk(ok, a.id ~= b.id && ~isnan(a.id) && ~isnan(b.id), ...
         sprintf('two claims get different units (%d and %d)', a.id, b.id));
ok = chk(ok, ~isempty(a.token) && ~strcmp(a.token, b.token) && a.lock.held, 'each claim carries its own token and holds its lock');
ok = chk(ok, startsWith(a.tmpOut, outF(a.id)) && endsWith(a.tmpOut, '.part'), 'the claim names an attempt-specific temporary output');

% ---- 2. done is the PUBLISHED artifact ----------------------------------
P0 = work_queue('publish', q, a, a.tmpOut);
ok = chk(ok, ~P0.ok && contains(P0.msg, 'does not exist'), 'publishing nothing is refused');
touchFile(a.tmpOut, 'bad');
P1 = work_queue('publish', q, a, a.tmpOut, @(f) deal(false, 'validator says no'));
ok = chk(ok, ~P1.ok && ~isfile(outF(a.id)) && isfile(a.tmpOut), 'a result that fails validation is NOT published');
P2 = work_queue('publish', q, a, a.tmpOut, @(f) deal(true, ''));
ok = chk(ok, P2.ok && isfile(outF(a.id)) && ~isfile(a.tmpOut), 'a validated result is published in one move');
work_queue('release', q, a);
ok = chk(ok, ~a.lock.held || ~unit_lock('probe', fullfile(q, sprintf('%d.lock', a.id))).held, 'release frees the lock');
a2 = work_queue('claim', q, 'A');
ok = chk(ok, a2.id ~= a.id, sprintf('a done unit is not handed out again (got %d, not %d)', a2.id, a.id));
s = work_queue('status', q);
ok = chk(ok, s.nDone == 1 && ismember(a.id, s.done), 'status reads DONE from the artifact, not a flag');
% A RELEASED CLAIM HAS NO AUTHORITY, whatever its struct says. Astra's
% pass-3 reproduction: A releases but keeps its struct; B claims; A's
% old struct must not beat, publish or release over B.
stale = a2;                                          % A2 holds unit a2.id right now
work_queue('release', q, stale);                     % ... and lets it go, keeping the struct
ok = chk(ok, stale.lock.held, '(the released struct still SAYS held -- that is the trap)');
bOwn = claimSpecific(q, 'B2', stale.id, 3);
ok = chk(ok, bOwn.id == stale.id, sprintf('B2 now owns unit %d', bOwn.id));
touchFile(stale.tmpOut, 'stale result');
P3 = work_queue('publish', q, stale, stale.tmpOut);
ok = chk(ok, ~P3.ok && contains(P3.msg, 'lock') && ~isfile(outF(stale.id)), 'the released claim cannot PUBLISH over the new owner');
ok = chk(ok, ~work_queue('beat', q, stale), 'nor beat');
ok = chk(ok, ~work_queue('release', q, stale) && isfile(fullfile(q, sprintf('%d.owner', stale.id))), 'nor release: B2''s owner record survives');
delete(stale.tmpOut);
work_queue('release', q, bOwn);

% ---- 3. release without publish returns the unit -------------------------
work_queue('release', q, b);
b2 = work_queue('claim', q, 'C');
ok = chk(ok, b2.id == b.id, sprintf('a unit released unpublished is handed out again (%d)', b2.id));
work_queue('release', q, b2);  work_queue('release', q, a2);

% ---- 4. open is read-only; attempts persist -------------------------------
attBefore = work_queue('status', q).attempts;
work_queue('init', q, 1:6, outF);                   % what a re-run does
s4 = work_queue('status', q);
ok = chk(ok, isequal(s4.attempts, attBefore) && any(attBefore > 0), ...
         're-opening keeps the attempt counts (wiping them re-armed the livelock)');
try
    work_queue('open', q, 1:8, outF);  opened = true;
catch
    opened = false;
end
ok = chk(ok, ~opened, 'open refuses a different unit set (a campaign''s units are fixed at creation)');
try
    work_queue('create', q, 1:6, outF);  created = true;
catch
    created = false;
end
ok = chk(ok, ~created, 'create refuses to overwrite an existing queue');

% ---- 5. retirement, reset, unreadable record --------------------------------
id5 = b.id;   % claimed twice already (B, C)
r = claimSpecific(q, 'F', id5, 3);  work_queue('release', q, r);    % attempt 3
r = work_queue('claim', q, 'F', struct('maxAtt', 3));
ok = chk(ok, r.id ~= id5, sprintf('after 3 attempts unit %d is RETIRED and not handed out', id5));
if ~isnan(r.id), work_queue('release', q, r); end
work_queue('init', q, 1:6, outF);
s5 = work_queue('status', q, struct('maxAtt', 3));
ok = chk(ok, ismember(id5, s5.retired), 'and it STAYS retired across a re-open');
h = work_queue('claim', q, 'H');                     % hold one unit during reset
attH = work_queue('status', q).attempts(work_queue('status', q).units == h.id);
reset_ = work_queue('reset', q, struct('units', [id5, h.id]));
ok = chk(ok, isequal(reset_, id5), sprintf('reset clears the retired unit (%d) and refuses the held one (%d)', id5, h.id));
attH2 = work_queue('status', q).attempts(work_queue('status', q).units == h.id);
ok = chk(ok, attH2 == attH && attH > 0, 'and the held unit''s attempt record is untouched (reset holds the lock while it acts)');
try
    work_queue('reset', q, struct('units', 99));  badReset = false;
catch
    badReset = true;
end
ok = chk(ok, badReset, 'reset refuses a unit id that is not in the queue');
work_queue('release', q, h);
s5b = work_queue('status', q, struct('maxAtt', 3));
ok = chk(ok, ismember(id5, s5b.todo), 'a reset unit is claimable again');
corrupt = fullfile(q, '6.att');  touchFile(corrupt, 'garbage');
s5c = work_queue('status', q, struct('maxAtt', 3));
ok = chk(ok, ismember(6, s5c.retired), 'an unreadable attempt record BLOCKS the unit rather than re-opening it');
delete(corrupt);

% ---- 6. one state per unit ----------------------------------------------------
g = work_queue('claim', q, 'G');
s6 = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
ok = chk(ok, s6.nDone + s6.nRunning + s6.nAbandoned + s6.nRetired + s6.nTodo == s6.n, ...
         sprintf('the five counts PARTITION the %d units (%d/%d/%d/%d/%d)', s6.n, ...
                 s6.nDone, s6.nRunning, s6.nAbandoned, s6.nRetired, s6.nTodo));
ok = chk(ok, ismember(g.id, s6.running) && s6.held(s6.units == g.id) && strcmp(s6.owner{s6.units == g.id}, 'G'), ...
         sprintf('a held unit is RUNNING with its owner named (%d by G)', g.id));
ok = chk(ok, ~s6.finished, 'a campaign with a held lock is not finished');
% an owner record WITHOUT a lock = a dead owner: abandoned, claimable
ageFile(fullfile(q, sprintf('%d.owner', g.id)));
unit_lock('release', g.lock.file, g.lock.token);      % simulate death: lock gone, record stays
s6b = work_queue('status', q, struct('maxAtt', 3, 'staleSec', 60));
ok = chk(ok, ismember(g.id, s6b.abandoned), sprintf('a unit whose owner died holding it is ABANDONED (%d)', g.id));
g2 = work_queue('claim', q, 'I');
ok = chk(ok, g2.id == g.id, 'and is handed out again');
work_queue('release', q, g2);

% ---- 7. heartbeat states ------------------------------------------------------
hb = fullfile(q, 'hb');
campaign_heartbeat('beat', hb, 'W1', 'working');
campaign_heartbeat('done', hb, 'W2', 'finished');
campaign_heartbeat('fail', hb, 'W3', 'threw');
campaign_heartbeat('beat', hb, 'W4', 'old');   ageFile(fullfile(hb, 'W4.hb'));
r7 = campaign_heartbeat('read', hb, {'W1','W2','W3','W4','W5'}, struct('staleSec', 60));
ok = chk(ok, isequal({r7.state}, {'running','done','failed','stalled','never'}), ...
         sprintf('five states told apart: %s', strjoin({r7.state}, ' ')));
ok = chk(ok, contains(r7(5).msg, 'no check-in'), 'a worker with no check-in is not reported as alive');
touchFile(fullfile(hb, 'W6.hb'), '');
touchFile(fullfile(hb, 'W7.hb'), 'beat');
touchFile(fullfile(hb, 'W8.hb'), 'beat|13-Sep-2026 12:00:00|working');   % 3-field, old format
r8 = campaign_heartbeat('read', hb, {'W6','W7','W8'}, struct('staleSec', 60));
ok = chk(ok, all(strcmp({r8.state}, 'unknown')), 'an empty, truncated or old-format heartbeat is UNKNOWN, not running');

% ---- 8. reporting cannot fail the work; fixed width ------------------------------
okr = safe_report(@() error('deliberate:report', 'the summary threw'), 'a summary');
ok = chk(ok, ~okr, 'safe_report catches a throwing report and returns false');
okr2 = safe_report(@() fprintf(''), 'a clean report');
ok = chk(ok, okr2, 'and returns true when the report is fine');
w = 8;
f1 = fmt_num(NaN, w, 4);  f2 = fmt_num([], w, 4);  f3 = fmt_num(17.7976, w, 4);  f4 = fmt_num(123456789.123, w, 4);
ok = chk(ok, numel(f1) == w && numel(f2) == w && numel(f3) == w && numel(f4) == w, ...
         sprintf('NaN, empty, a value and an OVER-WIDE value all render at width %d ("%s" "%s" "%s" "%s")', w, f1, f2, f3, f4));

% ---- 9. publish_atomic ----------------------------------------------------------
dst = fullfile(q, 'pub.txt');  touchFile(dst, 'old');
tmp = fullfile(q, 'pub.txt.tmp');  touchFile(tmp, 'new');
publish_atomic(tmp, dst);
ok = chk(ok, strcmp(fileread(dst), 'new') && ~isfile(tmp), 'publish_atomic replaces an existing file');
mkdir(fullfile(q, 'adir'));  touchFile(tmp, 'x');
try
    publish_atomic(tmp, fullfile(q, 'adir'));  nested = true;
catch
    nested = false;
end
ok = chk(ok, ~nested && ~isfile(fullfile(q, 'adir', 'pub.txt.tmp')), 'publish_atomic refuses a directory destination (movefile would nest into it)');

% ---- 10. unit_lock: authority lives in the registry ------------------------------
lf = fullfile(q, 'x.lock');
L = unit_lock('try', lf);
pr = unit_lock('probe', lf);
ok = chk(ok, L.held && pr.held && pr.mine, 'a held lock probes as held, and as MINE, without a second channel');
L2 = unit_lock('try', lf);
ok = chk(ok, ~L2.held && unit_lock('holds', lf, L.token), 'a second try from the same process is refused and leaves the first intact');
alias = fullfile(q, '.', 'x.lock');
ok = chk(ok, unit_lock('probe', alias).mine, 'an alias of the path maps to the same lock (identity, not string)');
n0 = unit_lock('count');
clear functions
ok = chk(ok, unit_lock('count') == n0 && unit_lock('holds', lf, L.token), 'the registry survives `clear functions` (it is mlocked)');
ok = chk(ok, unit_lock('release', lf, L.token) && ~unit_lock('release', lf, L.token), 'release is idempotent: true once, then a no-op');
ok = chk(ok, ~unit_lock('holds', lf, L.token) && ~unit_lock('probe', lf).held, 'and afterwards the token holds nothing and the file is free');

% ---- 11. rib_validate is unit-aware ---------------------------------------------
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'DRO_tulip', 'indirect'));
rf = fullfile(q, 'rib.mat');
pt = struct('ok', true, 'sD', 1/24, 'z', ones(8,1));
R = struct('j', 4, 'sA', 0.2, 'pts', [pt, setfield(pt, 'sD', 2/24)], 'stop', 'complete', 'nSolve', 2);
problem = struct('thrustN', 0.07, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1, 'sD', 0);
save(rf, 'R', 'problem');
spec = struct('nPts', 2, 'col', 4, 'sA', 0.2, 'nD', 24, 'problem', problem);
[okV, ~, inf_] = rib_validate(rf, spec);
ok = chk(ok, okV && inf_.complete && inf_.nPts == 2, 'a good rib for its unit validates as complete');
ok = chk(ok, ~rib_validate(rf, setfield(spec, 'col', 5)), 'a rib for ANOTHER column is refused');
R.pts(2).ok = false;  save(rf, 'R', 'problem');
ok = chk(ok, ~rib_validate(rf, spec), 'an uncertified point is refused');
R.pts(2).ok = true;  problem.NpTulip = 8;  save(rf, 'R', 'problem');  problem.NpTulip = 7;
ok = chk(ok, ~rib_validate(rf, spec), 'a rib certified at another operating point is refused');
R.pts(2).sD = 0.03;  save(rf, 'R', 'problem');
ok = chk(ok, ~rib_validate(rf, spec), 'a point off the departure lattice is refused');
% a worker validates its TEMPORARY result, whose name ends in .part: load
% without '-mat' treats that as ASCII and would have refused every real unit
R.pts(2).sD = 2/24;  rp = fullfile(q, 'rib.mat.abc.part');  save(rp, 'R', 'problem');
ok = chk(ok, rib_validate(rp, spec), 'a valid rib saved under a .part name validates (load -mat)');
fid = fopen(rf, 'w'); fprintf(fid, 'not a mat file'); fclose(fid);
[okV, msgV] = rib_validate(rf, spec);
ok = chk(ok, ~okV && ~isempty(msgV), 'a corrupt file is a refusal with a reason, not an error');

% ---- 12. walk_checkpoint ------------------------------------------------------
cf = fullfile(q, 'col.mat.ckpt');
ident = struct('sA', 0.2, 'nD', 24, 'dirn', -1, 'targets', -(1:23)/24, 'sD0', 0, 'problem', struct('thrustN', 0.07));
[C0, why0] = walk_checkpoint('load', cf, ident);
ok = chk(ok, isempty(C0) && strcmp(why0, 'no checkpoint'), 'no checkpoint yet: load says so');
walk_checkpoint('save', cf, struct('identity', ident, 'k', 5, 'sD', 0.7917, 'z', (1:8)', 'Y', ones(14, 25), 'pts', struct('sD', {0.9583, 0.9167}), 'nSolve', 12));
[C1, ~] = walk_checkpoint('load', cf, ident);
ok = chk(ok, ~isempty(C1) && C1.k == 5 && numel(C1.pts) == 2 && C1.nSolve == 12, 'a matching checkpoint loads with its state');
[C2, why2] = walk_checkpoint('load', cf, setfield(ident, 'sA', 0.25));
ok = chk(ok, isempty(C2) && contains(why2, 'another walk'), 'a checkpoint for another column is refused by name');
[C3, ~] = walk_checkpoint('load', cf, setfield(ident, 'targets', -(1:11)/12));
ok = chk(ok, isempty(C3), 'a checkpoint with other targets is refused');
walk_checkpoint('clear', cf);
ok = chk(ok, ~isfile(cf), 'clear removes it');

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
    held{end+1} = r;
end
for k = 1:numel(held), work_queue('release', q, held{k}); end
end

function touchFile(f, txt)
% TOUCHFILE  Write a small file.  INPUTS: f; txt.  OUTPUTS: none.
fid = fopen(f, 'w');  fprintf(fid, '%s', txt);  fclose(fid);
end

function ageFile(f)
% AGEFILE  Backdate a file so it reads as stale.  INPUTS: f.  OUTPUTS: none.
system(sprintf('touch -t 202001010000 "%s"', f));
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
