function setup_lunar_pnt_paths(extraChecks)
%% Purpose:
%
%   lunar_pnt_mission campaign paths: this dir + lib/ + the deliverable-9
%   kit folder (transfers_at, catalog_transfer, the 70 mN catalog) +
%   costate_common (get_family_orbit, phase_state, fly_transfer, ms_tfmin,
%   the verification instruments) + the top-level oclib/ + casadi under
%   $HOME + proj7/pipelines/gdop (gdop_grid_from_states, tulip_phasing,
%   the access scorers). pumpkyn/pumpkynPie come through pumpkynPie's own
%   startup only if pumpkyn.cr3bp.tfMin is not found (the cd around it is
%   undone by an onCleanup, even on error).
%
%   SHADOW POLICY (exact): the folders are RE-FRONTED with addpath -begin
%   in the order gdop, kit, costate_common, lib, this dir, so the final
%   path front is: this dir, lib, costate_common, kit, gdop. The kit is
%   fronted before costate_common so the kit's byte-identical
%   get_family_orbit copy cannot beat costate_common's. After re-fronting,
%   every function the mission uses is ASSERTED with which to resolve in
%   the folder it is meant to (including this campaign's own
%   mission_inputs); a mismatch throws setup_lunar_pnt_paths:shadow naming
%   the function and where it resolved. Re-fronting makes a path-order
%   shadow impossible, so the throw is the last line of defence (a same-name
%   file inside a re-fronted folder, a missing file); it is tested through
%   the extraChecks seam.
%
%   Cross-campaign effect: after this setup costate_catalog_pick resolves
%   to the kit copy, which differs from GTO_DRO/catalog's; this campaign
%   does not call it.
%
%% Inputs:
%
%  extraChecks  {M x 2} cell, optional      {functionName, requiredFolder}
%                                           pairs asserted in addition to
%                                           the built-in list (test seam)
%
%% Outputs:
%
%  none
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, extraChecks = cell(0, 2); end
proj7Root = fullfile(getenv('HOME'), 'Desktop', 'proj7');
here = fileparts(mfilename('fullpath'));
lib  = fullfile(here, 'lib');
ot   = fileparts(here);                        % orbit_transfer/
kit  = fullfile(ot, 'DRO_tulip', 'deliverables', ...
                'costate_library_70mN_user_deliverable9');
cc   = fullfile(ot, 'costate_common');
gd   = fullfile(proj7Root, 'pipelines', 'gdop');

addpath(fullfile(fileparts(ot), 'oclib'), fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = fullfile(proj7Root, 'external', 'pumpkynPie');
    od = cd(pp);
    back = onCleanup(@() cd(od));
    startup();
    clear back;
end
addpath(gd,  '-begin');
addpath(kit, '-begin');
addpath(cc,  '-begin');
addpath(lib, '-begin');
addpath(here, '-begin');

want = {'mission_inputs', lib; ...
        'transfers_at', kit; 'catalog_transfer', kit; ...
        'fly_transfer', cc; 'ms_tfmin', cc; 'get_family_orbit', cc; ...
        'phase_state', cc; 'gdop_grid_from_states', gd; ...
        'tulip_phasing', gd; 'max_min_access_bounds', gd; ...
        'max_min_access_block_bounds', gd; ...
        'schedule_rendezvous', lib; 'polish_rendezvous', lib};
want = [want; extraChecks];
for kk = 1:size(want, 1)
    got = fileparts(which(want{kk, 1}));
    if ~strcmp(got, want{kk, 2})
        error('setup_lunar_pnt_paths:shadow', ...
              '%s resolves to "%s", expected "%s"', ...
              want{kk, 1}, got, want{kk, 2});
    end
end
end
