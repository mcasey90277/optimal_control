function setup_paths()
% SETUP_PATHS  ms_band paths.
% Adds: ../lowThrust_GTO_tulip (lt_pmp_eom* indirect EOM), the DIRECT
% sundman_minfuel engine (cross-method edge: dual-.mat seed data under its
% results/ + solver helpers), and the shared CR3BP lib (params + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
