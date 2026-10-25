function ok = run_tests()
%% Purpose:
%
%   THE runner for lagrangianHamiltonianMechanics: every mechlib test and
%   every system test, library first, one logical verdict, and a THROW on
%   any failure so `matlab -batch "exit(~run_tests())"` returns a real
%   exit code. New systems append their tests to names below.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (also
%                                                   THROWS on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'mechlib'), fullfile(here, 'mechlib', 'tests'), ...
        fullfile(here, 'nonholonomic', 'ChaplyginSleigh'), ...
        fullfile(here, 'nonholonomic', 'ChaplyginSleigh', 'tests'));

names = {'test_lagrange_dalembert', ...      % the builder on exact problems
         'test_conservation_check', ...      % energy drift + constraint residual verdict
         'test_animate_planar_body', ...  % 1280 x 720 MP4 frames, refusals
         'test_sleigh_oracle', ...       % hand-derived sleigh vs builder, exact solution, energy, constraint
         'test_sleigh_pipeline'};        % sleigh body, plots, movie and the script's parse

res  = false(1, numel(names));
secs = zeros(1, numel(names));
for kk = 1:numel(names)
    fprintf('\n==== %s\n', names{kk});
    t0 = tic;
    try
        v = feval(names{kk});
        if ~isscalar(v)
            error('run_tests:badVerdict', '%s did not return a scalar verdict', names{kk});
        end
        res(kk) = logical(v);
    catch err
        fprintf('  THREW  %s (%s)\n', err.message, err.identifier);
        res(kk) = false;
    end
    secs(kk) = toc(t0);
end

fprintf('\n==== SUMMARY (lagrangianHamiltonianMechanics) ====\n');
for kk = 1:numel(names)
    if res(kk), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-28s %-4s %5.0f s\n', names{kk}, verdict, secs(kk));
end
ok = all(res);
fprintf('  %d of %d passed, %.0f s total\n', nnz(res), numel(res), sum(secs));
if ~ok
    error('run_tests:failed', '%d of %d tests FAILED: %s', nnz(~res), numel(res), strjoin(names(~res), ', '));
end
end
