%% REPRODUCE_LIBRARY_24X48_JOB  (24 x 48 variant, launched 2026-09-20)  The detached job that rebuilds the 70 mN DRO -> tulip
% 24 x 24 library (reproduce_library_70mN, go = true) and writes its VERDICT TO
% A FILE. Launch it with nohup, never in the interactive session:
%
%   nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/reproduce_job.out 2>&1 &
%
% It takes about a DAY (process/REPRODUCE_LIBRARY_GUIDE.md, section 4). Watch
% <outDir>/torus.log (one line per stage, silent for hours in between), this
% job's .out, and the verdict file below. It resumes: run it again.
%
% INPUTS:  none (edit the two paths below)
% OUTPUTS: <verdictF> (one line), <outDir>/reproduce_result.mat
verdictF = fullfile(getenv('HOME'), 'REPRODUCE_24x48_VERDICT.txt');
opts = struct('go', true, 'nD', 24, 'nA', 48, 'nWorkers', 8);   % arrival axis refined (FINDINGS 86, 88)

here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = fileparts(fileparts(mfilename('fullpath')));                   % DRO_tulip/indirect
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
cd(ind);
try
    out = reproduce_library_70mN(opts);
    save(fullfile(out.spec.outDir, 'reproduce_result.mat'), 'out');
    c = out.comparison;
    msg = sprintf(['REPRODUCE FINISHED: REPRODUCED = %d | campaign %s | same library %d, no worse %d: %d of %d cells agree; ' ...
                   'missing %d, extra %d, new faster %d, new slower %d, costates over %d, family %d, nonfinite %d'], ...
                  out.ok, out.state, c.ok, c.noWorse, c.nAgree, c.nRef, c.nOnlyRef, c.nOnlyNew, c.nNewFaster, c.nNewSlower, c.nZOver, c.nFamilyDiff, c.nNonfinite);
    if ~isempty(out.audit), msg = sprintf('%s | audit %d ok / %d bad', msg, out.audit.nOk, out.audit.nBad); else, msg = [msg ' | audit NOT ESTABLISHED']; end
catch ME
    msg = sprintf('REPRODUCE FAILED: %s | %s', ME.identifier, strrep(ME.message, newline, ' '));
    if ~isempty(ME.stack), msg = sprintf('%s | at %s line %d', msg, ME.stack(1).name, ME.stack(1).line); end
end
fid = fopen(verdictF, 'w');  fprintf(fid, '%s\n', msg);  fclose(fid);
fprintf('%s\n', msg);
