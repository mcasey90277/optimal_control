function cat_ = build_hh_catalog(direction)
%% Purpose:
%
%   Packages the L1 <-> L2 halo-to-halo sheets into the shareable compact
%   catalog via costate_common/build_costate_catalog_family -- the FIRST
%   catalog on the schema-v2 arrival-period axis (sheets keyed by
%   tau_dep x tau_arr; Np = NaN). Run any time; packages whatever sheets
%   exist (the validator enforces schema consistency).
%
%% Inputs:
%
%  direction                char                    'AtoB' L1 -> L2 sheets
%                                                   (hh_d*, default) |
%                                                   'BtoA' L2 -> L1 sheets
%                                                   (hhB_d*)
%
%% Outputs:
%
%  cat_                     struct                  The catalog, saved to
%                                                   direct/results/
%                                                   costate_catalog_halo_halo[_B].mat
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('direction','var'), direction = 'AtoB'; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'costate_common'));

if strcmpi(direction, 'BtoA')
    glob = 'hhB_d*.mat';
    name = 'costate_catalog_halo_halo_B';
    dirTxt = 'L2 -> L1';
else
    glob = 'hh_d*.mat';
    name = 'costate_catalog_halo_halo';
    dirTxt = 'L1 -> L2';
end

spec = struct( ...
  'glob', glob, ...
  'name', name, ...
  'description', sprintf(['Minimum-time low-thrust %s HALO-TO-HALO costate ', ...
    'CATALOG (Earth-Moon L1 and L2 southern halos), thrust rungs walked ', ...
    '15->1 N. FIRST schema-v2 catalog: sheets are keyed by (tau_dep x ', ...
    'tau_arr) -- both ARRIVAL and departure are halos, so the petal-count ', ...
    'field Np is NaN and costate_catalog_pick takes the ARRIVAL PERIOD ', ...
    '(ND) as its third argument. Each entry z8 = [lambda_r(3); ', ...
    'lambda_v(3); lambda_m; tf] is a converged PMP solution accepted ', ...
    'unchanged by pumpkyn.cr3bp.tfMin. Reconstruction recipes ride in ', ...
    'dep_family/dep_params and arr_family/arr_params.'], dirTxt), ...
  'provenance', sprintf(['Halo-to-halo catalog run (run_halo_halo_catalog, ', ...
    '%s stages): per sheet, 6x6 phasing torus anchored cold at 15 N and ', ...
    'walked down; three gates per entry (ms residual, flown arrival, tfMin ', ...
    'acceptance). Southern branches (pm = -1) both ends; the northern pair ', ...
    'is the exact z-mirror (CR3BP symmetry). Only two L1 members are ', ...
    'admissible under the 500 km / 100 Mm criteria. Casey/Koblick 2026-08.'], ...
    dirTxt), ...
  'depReconstruction', ['both ends: get_family_orbit(''halo'', ', ...
    'dep_params or arr_params) -- getHalo(tau, Lpt, pm)->cont_np->prop;  ', ...
    'state at phase f: interp1(t, rv, mod(f,1)*t(end), ''spline'')']);

catDir = fullfile(here, 'direct', 'results', 'catalog');
if strcmpi(direction, 'BtoA')
    outMat = fullfile(here, 'direct', 'results', 'costate_catalog_halo_halo_B.mat');
else
    outMat = fullfile(here, 'direct', 'results', 'costate_catalog_halo_halo.mat');
end
cat_ = build_costate_catalog_family(catDir, outMat, spec);
end
