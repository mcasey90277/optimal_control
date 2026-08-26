function setup_paths()
%% Purpose:
%
%   GTO_tulip/catalog module paths: this dir + costate_common (family
%   provider, preflight/harvest helpers) + the DRO_tulip engine tree
%   (direct transcription + its lib/certify helpers + the indirect
%   ms_tfmin-based thrust-ladder driver, unmodified) + the top-level
%   oclib/+oc package (see below) + casadi + pumpkyn/pumpkynPie (via
%   pumpkynPie's own startup, if not already on the path).
%   Mirrors HALO_tulip/run_halo_catalog.m's addpath set -- that driver is
%   the proven working reference; the brief's slimmer two-dir sketch left
%   casadi_mintime_dro and certify_dro_mintime unresolved.
%
%   ADDED 2026-08-26 (Task 3 pilot launch): certify_dro_mintime ->
%   dro_residual now calls oc.local_residual (repo-wide "oclib move 3",
%   landed on main after Task 2's smoke test), so this scaffold also needs
%   the top-level oclib/ dir (sibling of orbit_transfer/, NOT under it) on
%   the path -- MATLAB packages are resolved by adding the PARENT of a
%   '+pkg' folder. Without it every thrust_ladder_library call fails at
%   the certify step with "Unable to resolve the name 'oc.local_residual'"
%   -- confirmed this breaks HALO_tulip's own reference driver identically
%   (it never added oclib/ either); not fixed there, out of this task's
%   scope, but worth flagging.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
ot = fileparts(fileparts(here));
droDir = fullfile(ot, 'DRO_tulip');
addpath(here, fullfile(ot, 'costate_common'), ...
        fullfile(droDir, 'direct'), fullfile(droDir, 'direct', 'lib'), ...
        fullfile(droDir, 'direct', 'certify'), ...
        fullfile(droDir, 'indirect'), ...
        fullfile(fileparts(ot), 'oclib'), ...
        fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = '/Users/msc/Desktop/proj7/external/pumpkynPie';
    od = cd(pp); startup(); cd(od);
end
end
