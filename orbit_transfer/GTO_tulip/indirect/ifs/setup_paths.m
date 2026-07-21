function setup_paths()
% SETUP_PATHS  IFS paths.
% Adds: the DIRECT sundman_minfuel engine + its refine/ (cross-method edge:
% prep_refine_seed etc.), ../ms_band, ../lowThrust_GTO_tulip (lt_pmp_eom*),
% and the shared CR3BP lib (params + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel', 'refine'));
addpath(fullfile(here, '..', 'ms_band'));
addpath(fullfile(here, '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
