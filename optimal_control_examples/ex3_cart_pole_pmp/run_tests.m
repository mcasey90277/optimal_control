function ok = run_tests()
%% Purpose:
%
%   THE runner for this folder: run every test and FAIL THE PROCESS if any
%   check failed. Each test returns a logical and prints its own lines, which
%   is a fine local interface but a poor automation one -- a failing check
%   still let `matlab -batch "test_x"` exit 0, so a broken demo could sit
%   green in a script (GPT-6 Astra's review, 2026-09-17). This file closes
%   that: it errors on any failure, and `exit(~run_tests())` gives a shell a
%   real exit code.
%
%   Order is cheapest-first, with the physical oracle FIRST: if the dynamics
%   are wrong, every downstream number is meaningless and the run should say
%   so in its first line.
%
%   Usage:
%     run_tests                                            % interactive
%     matlab -batch "cd('<this folder>'); exit(~run_tests())"   % automation
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (the
%                                                   function also THROWS on
%                                                   failure, so a caller that
%                                                   ignores the output still
%                                                   fails)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(here, fullfile(here, 'tests'), fullfile(root, 'oclib'), ...
        fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'));

names = {'test_cartpole_physics'      % the independent oracle: geometry only
         'test_cartpole_field'        % dynamics vs the ex2 helpers
         'test_cartpole_state_jac'    % the generated Jacobian vs two methods
         'test_cartpole_pmp_rhs'      % the PMP field
         'test_direct_ref'            % the committed fixture
         'test_cartpole_pmp_prop'     % propagator, STM, collapse contract
         'test_cartpole_pmp'          % the solve, end to end
         'test_minenergy_study'};     % the study script, gates and verdict (slowest)

res = false(1, numel(names));
secs = zeros(1, numel(names));
for k = 1:numel(names)
    fprintf('\n==== %s\n', names{k});
    t0 = tic;
    try
        %% A non-scalar verdict (e.g. []) must not reach res(k) = v: with res
        %% logical and k a scalar index, res(k) = [] is MATLAB's element-
        %% DELETION syntax, not an assignment -- it would shrink res instead
        %% of recording a failure, and because test_minenergy_study is LAST
        %% in names (line 50), the shrink would surface as an uncaught
        %% out-of-bounds error in the summary loop below, not the named
        %% run_tests:failed (ported from run_all_tests.m, same defect).
        v = feval(names{k});
        if ~isscalar(v)
            error('run_tests:badVerdict', ...
                  '%s did not return a scalar verdict', names{k});
        end
        res(k) = logical(v);
    catch err
        fprintf('  THREW  %s (%s)\n', err.message, err.identifier);
        res(k) = false;
    end
    secs(k) = toc(t0);
end

fprintf('\n==== SUMMARY ====\n');
for k = 1:numel(names)
    fprintf('  %-28s %-4s %5.0f s\n', names{k}, tern(res(k), 'PASS', 'FAIL'), secs(k));
end
ok = all(res);
fprintf('  %d of %d passed, %.0f s total\n', nnz(res), numel(res), sum(secs));

if ~ok
    error('run_tests:failed', '%d of %d tests FAILED: %s', ...
          nnz(~res), numel(res), strjoin(names(~res), ', '));
end
end

% ------------------------------------------------------------------------
function s = tern(c, a, b)
%% Purpose:
%
%   a if c, else b.
%
if c, s = a; else, s = b; end
end
