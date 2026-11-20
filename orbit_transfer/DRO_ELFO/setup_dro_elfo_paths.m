function setup_dro_elfo_paths()
%% Purpose:
%
%   DRO_ELFO campaign paths: this dir, lib/ and catalog/ (each only when
%   the folder exists), costate_common (family provider get_family_orbit,
%   the shared ladder engine, pair_transfer, plot_pair_transfer), the
%   top-level oclib/ dir (parent of +oc), casadi, proj7/pipelines/gdop
%   (im_elfo_optimum) only when that function is not already on the path,
%   and pumpkyn / pumpkynPie through pumpkynPie's own startup when
%   pumpkyn.cr3bp.tfMin is missing. The current folder is restored after
%   the startup call.
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
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
ot   = fileparts(here);                        % orbit_transfer/
addpath(here);
for sub = {'lib', 'catalog'}
    if isfolder(fullfile(here, sub{1})), addpath(fullfile(here, sub{1})); end
end
addpath(fullfile(ot, 'costate_common'), fullfile(fileparts(ot), 'oclib'), ...
        fullfile(getenv('HOME'), 'casadi-3.7.0'));
if isempty(which('im_elfo_optimum'))
    addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'pipelines', 'gdop'));
end
if isempty(which('pumpkyn.cr3bp.tfMin'))
    pp = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkynPie');
    od = cd(pp); startup(); cd(od);
end
end
