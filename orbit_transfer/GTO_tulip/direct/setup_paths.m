function setup_paths()
% SETUP_PATHS  Paths for the direct GTO->tulip campaign (one campaign, one setup).
%
% Adds, in order: this folder (the front door), lib (solver + algorithms),
% certify (verification + second-order), viz (movies), tests, the shared CR3BP
% problem definitions, and the indirect PMP verifier the certify stage uses.
%
% ONE setup for the whole campaign. Before the 2026-07-26 flatten there were
% two -- sundman_minfuel/setup_paths.m and PSR/setup_paths.m -- for what was
% actually one pipeline plus its library, and the PSR one deliberately reached
% nothing outside itself. That isolation is gone, on purpose: it was already
% fictional (an addpath added in 2026-07-15 put the upstream solver ahead of
% PSR's frozen copies, so PSR had been running the upstream code for eleven
% days while its README claimed otherwise), and it cost a silently-dead feature
% and 20 duplicated files to maintain.
%
% THE INDIRECT DEPENDENCY IS REAL, NOT DUPLICATION. The certify stage runs the
% first-order PMP verifier (verify_direct_pmp + the sms_* costate machinery)
% from indirect/ms_band. The direct pipeline genuinely uses the indirect
% verifier -- that is the point of the check, an independent instrument. It was
% previously hidden by vendoring nine copies into PSR/lib; naming it here is
% more honest and keeps exactly one copy.
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));

addpath(here);                                   % front door + batch drivers
addpath(fullfile(here, 'lib'));                  % solver, homotopy, seeds, refine
addpath(fullfile(here, 'certify'));              % PMP + second-order verification
addpath(fullfile(here, 'viz'));                  % movies and frames
addpath(fullfile(here, 'tests'));                % guardrail tests

% Shared CR3BP problem definitions: cr3bp_lt_params, minfuel_config,
% gto_tulip_endpoints, insertion_states, cr3bp_ipopt_opts (+ pumpkyn).
addpath(fullfile(here, '..', '..', 'cr3bp_common'));
setup_cr3bp_common();

% Indirect PMP verifier used by the certify stage (see the note above).
addpath(fullfile(here, '..', 'indirect', 'ms_band'));
end
