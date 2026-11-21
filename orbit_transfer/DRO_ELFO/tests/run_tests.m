function ok = run_tests()
%% Purpose:
%
%   Every test of the DRO_ELFO user layer, with ONE logical verdict and a
%   throw on failure so a scripted run returns a real exit code. It runs
%   the generic plot test (costate_common/tests/test_plot_pair_transfer)
%   the user-script test, the fleet front-door test and the comparison-script test (synthetic sheets, seconds); the first two
%   solve once through tempdir caches (cold about 3.5 min each, seconds
%   when cached), the front-door test solves a 1 x 1 two-rung sheet
%   uncached (about 4 min), so a run belongs in the background. The
%   Task 7 tests (picker, clearance sweep, repair tool, chain) cut their
%   real fixture sheets from results/compare (seconds to a minute; a fresh
%   clone without that folder solves the fixtures instead, slower). The
%   library-function test reads the pilot catalog (seconds; SKIP without
%   it).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (the
%                                                   function also THROWS
%                                                   on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  test_run_dro_elfo_torus added                    (c) 11/19/2026
%  M. Casey  Task 7 tests (pick, clearance, repair, chain)    (c) 11/19/2026
%  M. Casey  test_dro_elfo_library (library functions)        (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here    = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common', 'tests'));

names = {'test_plot_pair_transfer', 'test_dro_elfo_transfer', 'test_run_dro_elfo_torus', ...
         'test_compare_dro_targets', 'test_dro_elfo_pick', 'test_body_clearance_sweep', ...
         'test_repair_dro_elfo_cells', 'test_build_dro_elfo_catalog', 'test_dro_elfo_library'};
res   = false(size(names));
secs  = zeros(size(names));
for kk = 1:numel(names)
    fprintf('\n######## %s\n', names{kk});
    t0 = tic;
    try
        res(kk) = feval(names{kk});
    catch err
        fprintf('  ERROR %s\n', err.message);
        res(kk) = false;
    end
    secs(kk) = toc(t0);
end

fprintf('\n%-28s %-6s %8s\n', 'test', 'result', 'seconds');
for kk = 1:numel(names)
    if res(kk), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-28s %-6s %8.1f\n', names{kk}, s, secs(kk));
end
ok = all(res);
fprintf('\nDRO_ELFO: %d of %d passed\n', nnz(res), numel(res));
if ~ok
    error('DRO_ELFO:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
