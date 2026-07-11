function setup_paths()
% SETUP_PATHS  Add IFS dependencies to the MATLAB path.
%
% Adds sundman_minfuel (cr3bp_lt_params, gto_tulip_endpoints, prep_refine_seed),
% sundman_minfuel/refine, ms_band (sms_eom, sms_seed_duals, sms_problem),
% lowThrust_GTO_tulip, and the pumpkyn parent (endpoint construction).
%
% INPUTS:  none
% OUTPUTS: none (path side effect)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'sundman_minfuel'));
addpath(fullfile(here, '..', 'sundman_minfuel', 'refine'));
addpath(fullfile(here, '..', 'ms_band'));
addpath(fullfile(here, '..', '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));
end
