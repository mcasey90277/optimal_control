function ok = run_tests()
%% Purpose:
%
%   THE runner for optimal_control_examples/drone: every test in the
%   folder, cheapest and most fundamental first (the physics oracle before
%   anything that trusts the field), one logical verdict, and a THROW on
%   any failure so `matlab -batch "exit(~run_tests())"` returns a real exit
%   code (the ex3 runner's contract). Later stages append their tests to
%   names below, after the library's.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(here, 'drone_dynamics'), fullfile(here, 'drone_dynamics', 'tests'), ...
        fullfile(root, 'oclib'));

names = {'test_drone_params'          % the constants
         'test_euler_zyx'             % attitude kinematics vs longhand + successive rotations
         'test_drone_physics'         % THE ORACLE: mechanics, not a mirror
         'test_drone_field'           % form and complex-step safety
         'test_drone_state_jac'       % generated Jacobian vs complex step + FD, regeneration
         'test_layouts'               % the paper's J_LB column, pair_list
         'test_clearance_terms'       % gaps, unit vectors, hermite_resample
         'test_barrier_terms'};       % the barrier and its derivatives

res = false(1, numel(names));
secs = zeros(1, numel(names));
for k = 1:numel(names)
    fprintf('\n==== %s\n', names{k});
    t0 = tic;
    try
        v = feval(names{k});
        if ~isscalar(v)
            error('run_tests:badVerdict', '%s did not return a scalar verdict', names{k});
        end
        res(k) = logical(v);
    catch err
        fprintf('  THREW  %s (%s)\n', err.message, err.identifier);
        res(k) = false;
    end
    secs(k) = toc(t0);
end

fprintf('\n==== SUMMARY (drone) ====\n');
for k = 1:numel(names)
    if res(k), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-28s %-4s %5.0f s\n', names{k}, verdict, secs(k));
end
ok = all(res);
fprintf('  %d of %d passed, %.0f s total\n', nnz(res), numel(res), sum(secs));
if ~ok
    error('run_tests:failed', '%d of %d tests FAILED: %s', nnz(~res), numel(res), strjoin(names(~res), ', '));
end
end
