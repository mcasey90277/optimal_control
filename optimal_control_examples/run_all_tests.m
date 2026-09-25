function ok = run_all_tests()
%% Purpose:
%
%   Every suite under optimal_control_examples, with ONE exit code. A test
%   that prints FAIL and returns false still exits 0 on its own, which is
%   fine at the prompt and useless in automation; this throws.
%
%   Folder order is dependency order: the shared plant first -- if the
%   physics is wrong, nothing downstream means anything -- then each example.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every suite passed (the
%                                                   function also THROWS on
%                                                   failure, so a -batch run
%                                                   returns a real exit code)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'oclib'));

%% The suites, in dependency order. Each entry is a folder and the functions
%% in it that return a logical verdict (an example's own run_tests is just
%% another such function):
suites = { 'cartpole_common',   {'test_cartpole_params', 'test_cartpole_physics', ...
                                 'test_cartpole_field',  'test_cartpole_state_jac'}
           'ex3_cart_pole_pmp', {'run_tests'}
           'conjugate_points',  {'run_conjugate_tests'} };

nT      = sum(cellfun(@numel, suites(:,2)));
names   = cell(1, nT);
res     = false(1, nT);
secs    = zeros(1, nT);
n       = 0;
for k = 1:size(suites, 1)
    folder = suites{k,1};
    addpath(fullfile(here, folder), fullfile(here, folder, 'tests'));
    for m = 1:numel(suites{k,2})
        n = n + 1;
        names{n} = [folder '/' suites{k,2}{m}];
        fprintf('\n######## %s\n', names{n});
        t0 = tic;
        try
            %% Name-collision guard: today there is exactly one run_tests in
            %% the tree, so bare-name feval dispatch is unambiguous. It will
            %% not stay that way -- once ex4/ex5 each ship their own
            %% run_tests.m, feval resolution rests on path order, which
            %% MATLAB's CURRENT-FOLDER precedence overrides regardless of
            %% addpath order. The ex3 README teaches readers to cd into the
            %% example folder first, so once a second run_tests exists,
            %% running from inside ex4 (say) could silently re-resolve to
            %% ex3's run_tests and report a false PASS in ex4's row. Nothing
            %% is broken by this collision today; this assertion exists so
            %% that when it does happen, it fails loudly here instead of
            %% passing quietly:
            fn = suites{k,2}{m};
            cand1 = fullfile(here, folder, [fn '.m']);
            cand2 = fullfile(here, folder, 'tests', [fn '.m']);
            if isfile(cand1)
                expected = cand1;
            elseif isfile(cand2)
                expected = cand2;
            else
                error('run_all_tests:missingFile', ...
                      '%s: no source file found at %s or %s', fn, cand1, cand2);
            end
            actual = which(fn);
            assert(strcmp(actual, expected), 'run_all_tests:nameCollision', ...
                   ['%s resolved to\n  %s\ninstead of the expected\n  %s\n' ...
                    'A same-named function elsewhere on the path (or in the ' ...
                    'current folder) is shadowing this suite.'], ...
                   fn, actual, expected);

            %% A non-scalar verdict (e.g. []) must not reach res(n) = v: with
            %% res logical and n a scalar index, res(n) = [] is MATLAB's
            %% element-DELETION syntax, not an assignment -- it would shrink
            %% res instead of recording a failure, silently or (for the last
            %% suite) with an uncaught out-of-bounds error downstream:
            v = feval(fn);
            if ~isscalar(v)
                error('run_all_tests:badVerdict', ...
                      '%s did not return a scalar verdict', names{n});
            end
            res(n) = logical(v);
        catch err
            fprintf('  THREW  %s\n', err.message);
            res(n) = false;
        end
        secs(n) = toc(t0);
    end
end

fprintf('\n==== ALL SUITES ====\n');
for k = 1:nT
    if res(k), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-44s %-4s %6.1f s\n', names{k}, tag, secs(k));
end
ok = all(res);
fprintf('  %d of %d passed, %.1f s total\n', nnz(res), nT, sum(secs));
if ~ok
    error('run_all_tests:failed', '%d suite(s) FAILED: %s', ...
          nnz(~res), strjoin(names(~res), ', '));
end
end
