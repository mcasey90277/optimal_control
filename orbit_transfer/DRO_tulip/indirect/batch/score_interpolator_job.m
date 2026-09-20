%% SCORE_INTERPOLATOR_JOB  The detached job that scores a costate library as a
% source of guesses (score_interpolator): random off-grid queries, each blended,
% polished under a hard cap, and compared with the nearest entry alone. Launch
% with nohup, never in the interactive session:
%
%   SCORE_LIBRARY=record nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch \
%       "run('<this file>')" > ~/score_record.out 2>&1 &
%
% SCORE_LIBRARY is 'record', 'spine96', or a catalog file's path; SCORE_NQUERY
% [60]. The verdict is <catalog's folder>/interp_score/SCORE_VERDICT.txt, one
% line per query in score.log beside it. It resumes: run it again.
%
% INPUTS:  none (environment variables above)
% OUTPUTS: the verdict file, score_rows.mat, score.log
here = pwd; cd('/Users/msc/Desktop/proj7/external/pumpkynPie'); startup(); cd(here);
ind = fileparts(fileparts(mfilename('fullpath')));                   % DRO_tulip/indirect
addpath(ind, fullfile(fileparts(fileparts(ind)), 'costate_common'));  cd(ind);
which_ = getenv('SCORE_LIBRARY');  if isempty(which_), which_ = 'record'; end
nQuery = str2double(getenv('SCORE_NQUERY'));  if ~isfinite(nQuery), nQuery = 60; end
switch which_
    case 'record',  catMat = fullfile(ind, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
    case 'spine96', catMat = fullfile(ind, 'results', 'sheet96_resolution_test', 'spine96_catalog.mat');
    otherwise,      catMat = which_;
end
try
    score_interpolator(catMat, struct('nQuery', nQuery));
catch ME
    outDir = fullfile(fileparts(catMat), 'interp_score');  if ~isfolder(outDir), mkdir(outDir); end
    fid = fopen(fullfile(outDir, 'SCORE_VERDICT.txt'), 'w');
    fprintf(fid, 'SCORE FAILED: %s | %s\n', ME.identifier, strrep(ME.message, newline, ' '));
    for q = 1:numel(ME.stack), fprintf(fid, '   at %s line %d\n', ME.stack(q).name, ME.stack(q).line); end
    fclose(fid);
end
