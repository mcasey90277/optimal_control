function setup_paths()
% SETUP_PATHS  Add ms_band dependencies to the MATLAB path.
%
% Adds: ../sundman_minfuel (cr3bp_lt_params, gto_tulip_endpoints, dual .mats),
% ../../lowThrust_GTO_tulip (lt_pmp_eom_minfuel + min-time indirect machinery),
% and the pumpkyn package parent (endpoint construction only).
%
% INPUTS:  none
% OUTPUTS: none (path side effect)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));  % pumpkyn parent — verified against ../sundman_minfuel/setup_paths.m
end
