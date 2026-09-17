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
           'ex3_cart_pole_pmp', {'run_tests'} };

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
            res(n) = logical(feval(suites{k,2}{m}));
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
