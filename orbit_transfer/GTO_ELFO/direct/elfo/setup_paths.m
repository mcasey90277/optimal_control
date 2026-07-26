function setup_paths()
% SETUP_PATHS  Direct GTO->ELFO campaign paths.
% Adds: self, and the shared CR3BP lib (cr3bp_lt_params, minfuel_config,
% gto_elfo_endpoints, insertion_states, cr3bp_ipopt_opts, pumpkyn).
%
% This campaign does NOT reuse tulip's solver. It has its own -- the forked
% casadi_energy_freetf / casadi_mintime_freetf, which carry the 9th state
% (cScale), free final time and the two-primary clock that ELFO needs.
%
% Until 2026-07-26 this file added all 34 files of tulip's direct Sundman
% engine to reach exactly ONE function, insertion_states, while its header
% claimed a reuse of casadi_minfuel_sundman / minfuel_at_tf that never
% happened -- those two names appeared only in comments. insertion_states now
% lives in cr3bp_common (where all of its own callees already were), so the
% tulip path is gone, and with it a shadowing surface: an ELFO session that
% had also put PSR on the path carried two different definitions of
% casadi_minfuel_sundman at once.
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
