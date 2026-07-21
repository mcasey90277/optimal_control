function setup_paths()
% SETUP_PATHS  Direct GTO->ELFO campaign paths.
% Adds: self, the tulip direct Sundman engine (cross-problem edge: this
% campaign reuses casadi_minfuel_sundman / insertion_states / minfuel_at_tf,
% retargeted to ELFO), and the shared CR3BP lib (params, gto_elfo_endpoints,
% pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'GTO_tulip', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
