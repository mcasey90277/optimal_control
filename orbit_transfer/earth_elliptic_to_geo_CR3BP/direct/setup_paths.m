function setup_paths()
% SETUP_PATHS  Path bootstrap for the elliptic->GEO CR3BP campaign.
%
% Adds this folder plus the 2-body campaign's core/ (shared solver, spec D3)
% so lunar_params / casadi_lt_mee / homotopy_mee / mee_seed all resolve.
%
% INPUTS:  none
% OUTPUTS: none (modifies the MATLAB path)
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md sec 4.
here = fileparts(mfilename('fullpath'));
addpath(here);
% Delegate to the 2-body campaign's OWN setup_paths (adds core/, lib/ --
% optdef.m lives there -- coords/, etc.) via the cwd trick: both files are
% named setup_paths.m, and cwd precedence resolves the local one (review
% amendment C, 2026-07-22).
e2b = fullfile(here, '..', '..', 'earth_elliptic_to_geo', 'direct');
oldd = cd(e2b);  setup_paths;  cd(oldd);
end
