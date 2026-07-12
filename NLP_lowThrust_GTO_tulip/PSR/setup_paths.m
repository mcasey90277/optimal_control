function setup_paths()
% SETUP_PATHS  Add the PSR pipeline's paths.
%
% PSR is SELF-CONTAINED as of 2026-07-12: all 19 machinery files the pipeline
% needs are VENDORED into PSR/lib (copies of ms_band + sundman_minfuel machinery;
% see PSR/lib/README.md for the manifest and provenance). This function adds
% only:
%   PSR         - the entry drivers (run_psr, psr_export_data, psr_movie)
%   PSR/lib     - the vendored machinery (solver, dual map, verifier, refine)
%   pumpkyn/src - the external tulip-construction toolbox (shared, third-party;
%                 lives in proj7, NOT vendored)
%
% The one input DATA dependency that is referenced in place (not vendored) is
% the min-energy backbone / bang-bang seed library under
% sundman_minfuel/results -- the PSR/lib copy of minfuel_config points there.
%
% CasADi is added by casadi_minfuel_sundman itself (CASADI_PATH env var or
% ~/casadi-3.7.0).
%
% INPUTS:  none
% OUTPUTS: none (path side effect)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'lib'));
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));
end
