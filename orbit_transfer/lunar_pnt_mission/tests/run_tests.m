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
%  M. Casey  test_export_mission_transfers added              (c) 11/10/2026
%  M. Casey  test_leg1_inputs added                           (c) 11/11/2026
%  M. Casey  test_polish_gto_dro added                        (c) 11/11/2026
%  M. Casey  test_schedule_gto_dro added                      (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
addpath(here);

names = {'test_mission_inputs', 'test_leg1_inputs', 'test_polish_gto_dro', 'test_schedule_gto_dro', 'test_schedule_rendezvous', 'test_polish_rendezvous', ...
         'test_assemble_histories', 'test_mission_gdop', 'test_validate_pole_series', ...
         'test_pole_series_stats', ...
         'test_lunar_pnt_mission_study', 'test_mission_movie', 'test_export_mission_transfers'};
res   = false(size(names));
skp   = false(size(names));                % a test that returns a second output true was SKIPPED
secs  = zeros(size(names));
for kk = 1:numel(names)
    fprintf('\n######## %s\n', names{kk});
    t0 = tic;
    try
        if nargout(names{kk}) >= 2
            [res(kk), skp(kk)] = feval(names{kk});
        else
            res(kk) = feval(names{kk});
        end
    catch err
        fprintf('  ERROR %s\n', err.message);
        res(kk) = false;
    end
    secs(kk) = toc(t0);
end

fprintf('\n%-30s %-6s %8s\n', 'test', 'result', 'seconds');
for kk = 1:numel(names)
    if skp(kk), s = 'SKIP'; elseif res(kk), s = 'PASS'; else, s = 'FAIL'; end
    fprintf('%-30s %-6s %8.1f\n', names{kk}, s, secs(kk));
end
ok = all(res);
fprintf('\nlunar_pnt_mission: %d of %d passed', nnz(res & ~skp), numel(res));
if any(skp), fprintf(', %d SKIPPED (%s)', nnz(skp), strjoin(names(skp), ', ')); end
fprintf('\n');
if ~ok
    error('lunar_pnt_mission:run_tests:fail', '%d test(s) failed', nnz(~res));
end
end
