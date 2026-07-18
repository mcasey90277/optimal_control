% RUN_TASK11_TESTS  Complete no-solve test suite for Task 11 deliverable
% verification. Runs every test_*.m in earth_elliptic_to_geo that does not
% require a fresh, unbounded NLP solve (per Task-11 brief's named list),
% in one MATLAB launch, and reports PASS/FAIL per file plus a final summary.
%
% INPUTS:  none
% OUTPUTS: none (prints PASS/FAIL per test + summary; errors if any FAIL)
hereTop = fileparts(mfilename('fullpath'));
cd(hereTop);

testList = {'test_mee_rhs', 'test_mee_seed', 'test_mee_solver_smoke', ...
            'test_warmstart_mee', 'test_mintime_mee_guard', 'test_run_ladder', ...
            'test_psr_mee', 'test_verify_pmp_mee'};

nTests  = numel(testList);
passArr = false(1, nTests);
msgArr  = cell(1, nTests);
secArr  = zeros(1, nTests);

for kt = 1:nTests
    fprintf('\n================= RUNNING %s =================\n', testList{kt});
    [passArr(kt), msgArr{kt}, secArr(kt)] = run_one_test(hereTop, testList{kt});
end

fprintf('\n================= SUMMARY =================\n');
for kt = 1:nTests
    status = 'FAIL'; if passArr(kt), status = 'PASS'; end
    fprintf('%-6s %-28s %6.2f s  %s\n', status, testList{kt}, secArr(kt), msgArr{kt});
end
fprintf('%d/%d PASS\n', sum(passArr), nTests);
if ~all(passArr)
    error('run_task11_tests:someFailed', 'Not all no-solve tests passed.');
end

function [ok, msg, secs] = run_one_test(hereTop, testName)
% RUN_ONE_TEST  Executes one test_*.m script in its own function workspace
% (isolated from the caller's loop bookkeeping) and captures pass/fail.
t0 = tic;
try
    run(fullfile(hereTop, [testName '.m']));
    ok  = true;
    msg = 'OK';
catch ME
    ok  = false;
    msg = ME.message;
end
secs = toc(t0);
end
