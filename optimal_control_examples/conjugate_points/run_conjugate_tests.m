function ok = run_conjugate_tests()
%% Purpose:
%
%   THE runner for conjugate_points/: every test, one verdict, and a THROW
%   on failure so `matlab -batch` returns a real exit code. Named apart from
%   ex3's run_tests so the two never shadow each other on the path
%   (run_all_tests guards that collision).
%
%   Order: closed-form oracles first (if the oscillator is wrong nothing
%   else means anything), then the independent geometric oracle (Lindelof),
%   then the cross-instrument statements, the deliberate-failure
%   experiments, and the GUI last.
%
%   Usage:
%     run_conjugate_tests
%     matlab -batch "cd('<this folder>'); run_conjugate_tests"
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 every test passed (also
%                                                   throws on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'tests'));
names = {'test_cov_oscillator'      % closed forms: sin, cot, k pi, eigenvalues
         'test_cov_catenary'        % closed forms + Lindelof's tangents
         'test_cov_instruments'     % Morse, r'(p) = h(b), Delta J, delta -> 0
         'test_cov_mutation'        % each check fails on the input it targets
         'test_cov_gui'};           % the explorer says what the instruments say
res = false(1, numel(names));
secs = zeros(1, numel(names));
for k = 1:numel(names)
    fprintf('\n==== %s\n', names{k});
    t0 = tic;
    try
        v = feval(names{k});
        if ~isscalar(v)
            error('run_conjugate_tests:badVerdict', '%s returned a non-scalar verdict', names{k});
        end
        res(k) = logical(v);
    catch err
        fprintf('  THREW  %s\n', err.message);
    end
    secs(k) = toc(t0);
end
fprintf('\n==== conjugate_points ====\n');
for k = 1:numel(names)
    if res(k), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-26s %s  %6.1f s\n', names{k}, tag, secs(k));
end
ok = all(res);
if ~ok
    error('run_conjugate_tests:failed', '%d test(s) FAILED: %s', nnz(~res), ...
          strjoin(names(~res), ', '));
end
end
