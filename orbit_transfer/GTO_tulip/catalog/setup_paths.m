function setup_paths()
%% Purpose:
%
%   GTO_tulip/catalog module paths: this dir + costate_common (family
%   provider, preflight/harvest helpers) + the DRO_tulip engine tree
%   (direct transcription + its lib/certify helpers + the indirect
%   ms_tfmin-based thrust-ladder driver, unmodified) + casadi + pumpkyn/
%   pumpkynPie (via pumpkynPie's own startup, if not already on the path).
%   Mirrors HALO_tulip/run_halo_catalog.m's addpath set exactly -- that
%   driver is the proven working reference; the brief's slimmer two-dir
%   sketch left casadi_mintime_dro and certify_dro_mintime unresolved.
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
        fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = '/Users/msc/Desktop/proj7/external/pumpkynPie';
    od = cd(pp); startup(); cd(od);
end
end
