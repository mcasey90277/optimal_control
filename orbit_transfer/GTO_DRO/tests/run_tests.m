function ok = run_tests()
%% Purpose:
%
%   Every test in this folder, in dependency order, with ONE logical
%   verdict and a throw on failure so a scripted run returns a real exit
%   code. The study-script test solves one transfer walking the ladder
%   15 -> 12 N and takes ~2.5 minutes; the others are seconds.
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
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);

names = {'test_earth_clearance', 'test_gto_dro_entry', 'test_gto_dro_study'};
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
fprintf('\nGTO_DRO: %d of %d passed\n', nnz(res), numel(res));
if ~ok
    error('GTO_DRO:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
