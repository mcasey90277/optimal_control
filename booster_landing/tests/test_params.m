% TEST_PARAMS  Physical sanity of the booster parameter set.
%
% Checks the three facts the whole campaign leans on:
%   1. m0 > mdry (there is landing propellant);
%   2. min-throttle thrust-to-weight at DRY mass > 1 (hoverslam is forced:
%      the booster cannot hover, so the optimizer must produce a terminal
%      max-thrust arc arriving at v=0 exactly at touchdown);
%   3. bounds ordered: 0 < Tmin < Tmax, tf_lo < tf_hi, gs in (0,90) deg.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P = booster_params();
assert(P.m0 > P.mdry, 'no landing propellant');
twMin = P.Tmin / (P.mdry * P.g0);
assert(twMin > 1, 'min-throttle T/W at dry mass = %.3f, expected > 1', twMin);
assert(P.Tmin > 0 && P.Tmin < P.Tmax, 'thrust bounds disordered');
assert(P.tf_lo < P.tf_hi, 'tf bracket disordered');
assert(P.gs_deg > 0 && P.gs_deg < 90, 'glideslope out of range');
assert(P.drag.on == false, 'Phase 1 default must be vacuum');
fprintf('test_params PASS (min-throttle dry T/W = %.3f)\n', twMin);
