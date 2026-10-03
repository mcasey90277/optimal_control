%% DEPARTURE_RIB96_JOB  (2026-10-03)  Measure the DEPARTURE axis the way the
% arrival axis was measured (FINDINGS 86, 90): one rib at 96 departure
% phases on one arrival column, so the departure-phase resolution an
% interpolated costate guess needs can be read off (INTERPOLATION_STATUS
% next step 1). Along departure at 1/24 the costates change by ~32-36% per
% step (FINDINGS 84-85); the 24 x 48 scored 1 of 60 because of it.
%
% The column: arrival phase 0.2837 (spine96 column 28 = record column 7),
% the smoothest departure ring of the record -- one family (anchor), 23 of
% 24 departure edges time-consistent, worst costate step 43%. A rib walked
% at 1/96 from sD = 0 passes the record's 24 departure phases at every 4th
% point, which doubles as a reproduction check.
%
% Launch with nohup, never in the interactive session:
%   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/departure_rib96.out 2>&1 &
% Progress: <outDir>/heartbeat.txt gets one line per solve; the checkpoint
% is rewritten after every accepted point. It resumes from the checkpoint:
% run it again.
%
% INPUTS:  none (edit the block below)
% OUTPUTS: <verdictF> (one line); <outDir>/rib96_col28.mat (build_ribs R)
verdictF = fullfile(getenv('HOME'), 'DEPARTURE_RIB96_VERDICT.txt');
res      = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results';
sheetMat = fullfile(res, 'sheet96_resolution_test', 'arrival_sheet_70mN_nA96.mat');
col      = 28;                                  % spine96 column: sA 0.2837
nD       = 96;
outDir   = fullfile(res, 'departure_rib96_test');

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect';
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    if ~isfolder(outDir), mkdir(outDir); end
    hb = fullfile(outDir, 'heartbeat.txt');
    beat = @() writelines(string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')) + " solve", hb, 'WriteMode', 'append');
    out = fullfile(outDir, sprintf('rib96_col%02d.mat', col));
    R = build_ribs(sheetMat, struct('nD', nD, 'nPts', nD - 1, 'direction', -1, 'only', col, 'wallSec', 900, ...
                                    'out', out, 'checkpoint', fullfile(outDir, sprintf('rib96_col%02d_ckpt.mat', col)), ...
                                    'progress', beat));
    assert(~isempty(R), 'departure_rib96_job: no rib was walked (column %d has no certified crossing?)', col);
    msg = sprintf('RIB96 FINISHED: column %d (sA %.4f): %d of %d points certified, %d solves, stop: %s', ...
                  col, R(1).sA, numel(R(1).pts), nD - 1, R(1).nSolve, R(1).stop);
catch ME
    msg = sprintf('RIB96 FAILED: %s | %s', ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);
