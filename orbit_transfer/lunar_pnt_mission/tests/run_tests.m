function ok = run_tests()
%% Purpose:
%
%   Every test in this folder, in dependency order, with ONE logical
%   verdict and a throw on failure so a scripted run returns a real exit
%   code.
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
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
addpath(here);

names = {'test_mission_inputs', 'test_schedule_rendezvous', 'test_polish_rendezvous', ...
         'test_assemble_histories'};
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
fprintf('\nlunar_pnt_mission: %d of %d passed\n', nnz(res), numel(res));
if ~ok
    error('lunar_pnt_mission:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
