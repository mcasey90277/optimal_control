function ok = test_run_capped()
%% Purpose:
%
%   Tests run_capped -- the parfeval hard per-call timeout -- on a real
%   pool (capped_pool). Checks:
%     1. FINISHES IN TIME: ok = true and all nout outputs are fcn's, in
%        order, with arguments passed through.
%     2. WORKER ERROR: ok = false and every output is [].
%     3. TIMEOUT: a call that outlives its cap returns ok = false and []
%        outputs within the cap plus a small margin -- not after the call.
%     4. INDISTINGUISHABLE: a timeout and an error return the same thing,
%        as the README records ("the caller cannot tell which").
%     5. THE POOL SURVIVES a cancelled call: the next call succeeds.
%   With no pool available (no Parallel Computing Toolbox, or its licence
%   held elsewhere) the test SKIPS and returns true, as the fence itself
%   degrades.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed (or
%                                                   skipped: no pool)
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
pool = capped_pool();
if isempty(pool), fprintf('  SKIP  no parallel pool available\n'); return, end

%% 1. Finishes in time:
[okA, s, p] = run_capped(pool, @(a, b) deal(a + b, a * b), 2, 60, 3, 4);
ok = chk(ok, okA && isequal(s, 7) && isequal(p, 12), ...
         sprintf('finishes: ok, outputs in order (%g, %g), arguments passed', s, p));

%% 2. Worker error:
[okE, e1, e2] = run_capped(pool, @(x) error('test:boom', 'boom %d', x), 2, 60, 1);
ok = chk(ok, ~okE && isempty(e1) && isempty(e2), 'a worker error: ok = false, outputs []');

%% 3. Timeout:
capSec = 2;
t0 = tic;
[okT, t1] = run_capped(pool, @(x) sleepThen(x), 1, capSec, 60);
wall = toc(t0);
ok = chk(ok, ~okT && isempty(t1), 'a call past its cap: ok = false, output []');
ok = chk(ok, wall < capSec + 30, sprintf('returned after %.1f s for a 60 s call capped at %g s', wall, capSec));

%% 4. Indistinguishable:
ok = chk(ok, isequal(okT, okE) && isequal(t1, e1), 'a timeout and an error look the same to the caller');

%% 5. The pool survives:
[okS, v] = run_capped(pool, @(x) 2*x, 1, 60, 21);
ok = chk(ok, okS && v == 42, 'the next call after a cancellation succeeds');

if ok, fprintf('TEST_RUN_CAPPED: ALL PASS\n');
else,  fprintf('TEST_RUN_CAPPED: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function y = sleepThen(x)
%% Purpose:
%
%   Sleep x seconds, then return x (a call that outlives its cap).
%
pause(x);
y = x;
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
