%% JOB_TABLE1  Batch driver for drone_table1 (the full table, incl. stretch).
%
%  Launched by run_table1.sh in a clean MATLAB: adds its own paths, puts
%  CasADi on the path, logs to a file beside the results and writes a
%  JOB DONE marker. Edit the Ms line for the stretch rows (64, 128 take a
%  long time: thousands of pair constraints).
%
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------
here    = fileparts(mfilename('fullpath'));
route1  = fileparts(here);
addpath(route1, fullfile(fileparts(route1), 'drone_dynamics'));
logFile = fullfile(route1, 'results', 'table1_batch.log');
if ~isfolder(fullfile(route1, 'results')), mkdir(fullfile(route1, 'results')); end
try
    casadi_path;
    drone_table1(struct('Ms', [4 8 16 32], 'logFile', logFile));
catch err
    fid = fopen(logFile, 'a');
    fprintf(fid, 'JOB FAILED: %s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
    fclose(fid);
end
fid = fopen(logFile, 'a');  fprintf(fid, 'JOB DONE\n');  fclose(fid);
