function setup_gto_dro_paths()
%% Purpose:
%
%   GTO_DRO campaign paths: this dir + lib/ + costate_common (family
%   provider get_family_orbit, the shared ladder engine
%   thrust_ladder_library / casadi_mintime_dro / certify_dro_mintime, the
%   verification instruments fly_transfer / pmp_pointwise_checks /
%   mintime_hypothesis_gates / lift_margin / conj_spectrum) + the top-level
%   oclib/ dir (parent of +oc; certify_dro_mintime -> dro_residual calls
%   oc.local_residual) + casadi + pumpkyn/pumpkynPie through pumpkynPie's
%   own startup if they are not already on the path. No DRO_tulip or
%   GTO_tulip folder is added: the engine left DRO_tulip on 2026-09-16 and
%   GTO_tulip/catalog's helpers are tulip-keyed.
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
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
ot   = fileparts(here);                        % orbit_transfer/
addpath(here, fullfile(here, 'lib'), fullfile(ot, 'costate_common'), ...
        fullfile(fileparts(ot), 'oclib'), ...
        fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = '/Users/msc/Desktop/proj7/external/pumpkynPie';
    od = cd(pp); startup(); cd(od);
end
end
