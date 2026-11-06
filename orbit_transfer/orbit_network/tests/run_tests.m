function ok = run_tests()
%% Purpose:
%
%   Every test in this folder, in dependency order (loader, node table,
%   lib helpers, then the scripts end to end), with ONE logical verdict
%   and a throw on failure so a -batch run returns a real exit code.
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
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here, root, fullfile(root, 'lib'), fullfile(fileparts(root), 'costate_common'));

names = {'test_load_catalog_sheets', 'test_orbit_node_table', 'test_lib', 'test_pipeline'};
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
fprintf('\norbit_network: %d of %d passed\n', nnz(res), numel(res));
if ~ok
    error('orbit_network:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
