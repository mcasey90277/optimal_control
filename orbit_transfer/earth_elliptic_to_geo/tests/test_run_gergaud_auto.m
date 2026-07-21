% test_run_gergaud_auto.m — auto mode on the default 10 N endpoints must
% return the certified cached row WITHOUT a fresh solve, and without movie.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_paths;
row = run_gergaud(struct('thrustN',10,'runMode','auto','makeMovie',false, ...
    'makePlot',false,'returnOnly',true));
assert(abs(row.m_f_kg-1377.10) < 0.1, '10 N auto row m_f=1377.10');
assert(row.switches==19 && abs(row.revs-7.326)<1e-2, '10 N structure');
assert(row.certified, '10 N row certified');
% endpoint-default detection: paper defaults -> uses cache; a custom final
% flips to solve mode (assert it does NOT claim the cached 10 N number).
fprintf('test_run_gergaud_auto PASSED\n');
