function setup_paths()
% SETUP_PATHS  Add the GTO->ELFO direct min-fuel pipeline's paths.
%
% elfo/ is a self-contained deliverable directory for the GTO->ELFO transfer,
% the sibling of PSR/ (which is the GTO->tulip deliverable). Unlike PSR -- which
% VENDORS its machinery into PSR/lib -- elfo/ uses a SHARED-PATH model: the two
% shared engine files (cr3bp_lt_params, minfuel_config) stay single-source in
% sundman_minfuel/ and are added to the path here. Nothing is copied, so there
% is no vendoring-drift surface. Everything ELFO-specific (the two-primary
% free-tf solver, seed generators, drivers, endpoints, export/verify, movie)
% lives in elfo/ itself.
%
% Paths added:
%   elfo             - the ELFO pipeline (solver, drivers, run_elfo_minfuel)
%   sundman_minfuel  - shared engine (cr3bp_lt_params, minfuel_config) + seed
%                      data under sundman_minfuel/results (referenced in place)
%   pumpkyn/src      - external tulip-construction toolbox (third-party, in proj7)
%
% CasADi is added by casadi_energy_freetf itself (CASADI_PATH env var or
% ~/casadi-3.7.0).
%
% INPUTS:  none
% OUTPUTS: none (path side effect)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'sundman_minfuel'));
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));
end
