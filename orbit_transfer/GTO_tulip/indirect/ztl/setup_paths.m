function setup_paths()
% SETUP_PATHS  ztl paths: the base indirect campaign (now a sibling under
% indirect/) + the shared CR3BP lib (params + pumpkyn, asserted there).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
oldCampaign = fullfile(here, '..', 'lowThrust_GTO_tulip');
assert(exist(fullfile(oldCampaign, 'lt_pmp_eom_minfuel.m'), 'file') == 2, ...
    'setup_paths:missing', 'indirect campaign not found at %s', oldCampaign);
addpath(oldCampaign);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
