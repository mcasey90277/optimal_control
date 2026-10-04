function ok = test_fill_holes_pool()
% TEST_FILL_HOLES_POOL  The hole filler keeps its fence alive (FINDINGS 92):
% 47 of the 24 x 48 filler's 58 failures were "The parallel pool has shut
% down" -- its direct solves run in-process, the pool idled past its
% 30-minute timeout, and every fenced certification after that threw.
% Checks: a fence pool never idles out; a dead pool is reopened; no pool is
% created when none was wanted (no PCT); the main loop revives the pool
% before each cell and hands the live one to the solver.
%
% Opens real one-worker pools: run in a batch process, not in a shared
% session (capped_pool would adopt that session's pool).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
H = fill_holes_direct('localfunctions');
delete(gcp('nocreate'));

% ---- a fence pool never idles out --------------------------------------
p = capped_pool(1);
ok = chk(ok, ~isempty(p) && isvalid(p) && isinf(p.IdleTimeout), 'capped_pool: a new pool has IdleTimeout = Inf');

% ---- a live pool is kept, and made idle-proof ----------------------------
p.IdleTimeout = 30;
q = H.livePool(p, true);
ok = chk(ok, q == p && isinf(q.IdleTimeout), 'livePool: a live pool is kept (same object), IdleTimeout set to Inf');

% ---- a dead pool is reopened ---------------------------------------------
delete(p);
q = H.livePool(p, true);
ok = chk(ok, ~isempty(q) && isvalid(q) && q.Connected && isinf(q.IdleTimeout), 'livePool: a deleted pool is reopened, idle-proof');
f = parfeval(q, @plus, 1, 2, 3);
ok = chk(ok, fetchOutputs(f) == 5, 'livePool: the reopened pool runs a call');
delete(q);

% ---- no fence wanted, none created ---------------------------------------
q = H.livePool([], false);
ok = chk(ok, isempty(q) && isempty(gcp('nocreate')), 'livePool: no pool is created when none was wanted');

% ---- wiring: the loop revives before each cell -----------------------------
src = fileread(fullfile(here, 'fill_holes_direct.m'));
loop = extractBetween(src, 'while kc < size(cells, 1)', 'function s = join_note');
ok = chk(ok, ~isempty(loop) && contains(loop{1}, 'pool = livePool(pool, wantFence)') ...
         && contains(loop{1}, 'solveOpts.pool = pool'), 'wiring: the cell loop revives the pool and passes it on');

% ---- the filler keeps its candidates (optimality-status spec 5) -----------
C1 = struct('ok', false, 'status', 3, 'z', ones(8, 1));  C0 = struct('ok', false, 'status', -1, 'z', nan(8, 1));
o = H.keepCandidate(struct([]), C1);  o = H.keepCandidate(o, C0);
ok = chk(ok, numel(o) == 1 && o(1).status == 3, 'keepCandidate keeps a status >= 1 result, drops one below the floor');
o = H.keepCandidate(o, C1);
ok = chk(ok, numel(o) == 1, 'keepCandidate keeps a repeated result once');
ok = chk(ok, contains(loop{1}, 'others = keepCandidate(others, C)') && contains(src, '''others'''), ...
         'wiring: the seed loop keeps candidates and the cell record stores them');

if ok, fprintf('test_fill_holes_pool: ALL PASS\n'); else, fprintf('test_fill_holes_pool: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
