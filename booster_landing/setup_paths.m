function setup_paths()
% SETUP_PATHS  Paths for the booster_landing campaign (one campaign, one setup).
%
% Adds: this folder (front door), lib (dynamics + solvers + tracking),
% certify (gates), viz (plots + movie), tests, and CasADi 3.7.
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'lib'));
addpath(fullfile(here, 'certify'));
addpath(fullfile(here, 'viz'));
addpath(fullfile(here, 'tests'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));
addpath(fullfile(fileparts(here), 'oclib'));   % oc.* shared OC library
end
