function d = batch_log_dir()
% BATCH_LOG_DIR  Folder for detached batch-job logs and one-line verdict files.
%
%   The *_job.m scripts write their completion verdicts here, and their
%   launch lines redirect stdout here, instead of into the home folder
%   (where ~60 of these files had piled up by Oct 2026). Created on first
%   call; its contents are git-ignored.
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   d - absolute path to DRO_tulip/indirect/results/batch_logs [char]
d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'batch_logs');
if ~isfolder(d)
    mkdir(d);
end
end
