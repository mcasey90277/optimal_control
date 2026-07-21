function setup_paths()
% SETUP_PATHS  Add the min_time module's paths.
%
% Self-contained min-time (PMP, always-burn) CR3BP transfer code. Adds:
%   min_time     - the module (mintime_solve, mintime_params, drivers)
%   pumpkyn/src  - external toolbox: analytic-STM min-time propagator
%                  (pumpkyn.cr3bp.tfMinProp/tfMinEoM/minDeltaV) + orbit builders
%                  (getTulip, fromOrb, orb2eci, fromPCI). Shared, third-party.
%
% INPUTS:  none    OUTPUTS: none (path side effect)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));
end
