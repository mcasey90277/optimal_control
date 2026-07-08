function setup_paths()
% SETUP_PATHS  Add the pumpkyn toolbox to the MATLAB path for this tutorial.
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   (none) - modifies the MATLAB path in-place
%
% The pumpkyn clone lives in the proj7 tree; only its src/ folder (which
% holds the +pumpkyn package) is required here.

pumpkynSrc = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
                      'pumpkyn', 'src');
if ~exist(fullfile(pumpkynSrc, '+pumpkyn'), 'dir')
    error('setup_paths:missing', 'pumpkyn not found at %s', pumpkynSrc);
end
addpath(pumpkynSrc);
end
