function setup_paths()
% SETUP_PATHS  Add the pumpkyn toolbox to the MATLAB path for this solver.
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   (none) - modifies the MATLAB path in-place
%
% pumpkyn is used only for problem setup (GTO state construction, tulip
% seed, ballistic propagation) and for the indirect warm-start option; the
% NLP transcription and solve are entirely local code.

pumpkynSrc = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
                      'pumpkyn', 'src');
if ~exist(fullfile(pumpkynSrc, '+pumpkyn'), 'dir')
    error('setup_paths:missing', 'pumpkyn not found at %s', pumpkynSrc);
end
addpath(pumpkynSrc);
end
