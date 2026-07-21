function setup_paths()
% SETUP_PATHS  Paths for the ZTL (Zhang-style thrust ladder) campaign.
%
% Adds the old indirect campaign (lowThrust_GTO_tulip: lt_pmp_eom*,
% solve_*_indirect, shoot_residual_*) and the pumpkyn toolbox (endpoints).
%
% INPUTS:
%   (none)
% OUTPUTS:
%   (none) - modifies the MATLAB path in-place

here = fileparts(mfilename('fullpath'));
oldCampaign = fullfile(here, '..', '..', 'lowThrust_GTO_tulip');
assert(exist(fullfile(oldCampaign, 'lt_pmp_eom_minfuel.m'), 'file') == 2, ...
    'setup_paths:missing', 'old indirect campaign not found at %s', oldCampaign);
addpath(oldCampaign);

pumpkynSrc = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
                      'pumpkyn', 'src');
assert(exist(fullfile(pumpkynSrc, '+pumpkyn'), 'dir') == 7, ...
    'setup_paths:missing', 'pumpkyn not found at %s', pumpkynSrc);
addpath(pumpkynSrc);
end
