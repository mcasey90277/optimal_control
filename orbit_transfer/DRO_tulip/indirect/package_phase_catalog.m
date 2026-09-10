function cat_ = package_phase_catalog(sheetMat, ribMats, opts)
%% Purpose:
%
%   Package a certified arrival sheet and its departure ribs into a
%   shareable costate CATALOG, in the same compact format as the six
%   min-time catalogs already shipped to Darin.
%
%   This regime needs its own catalog rather than a new rung on the shipped
%   DRO -> tulip one: a catalog carries ONE thruster, and this engine is
%   Isp 900 s against the shipped catalog's 1710 s. It needs no new schema --
%   the min-time schema has always keyed sheets by sD_frac x sA_frac with a
%   thrust axis, so a one-rung phase sheet is a layout conversion
%   (sheet_to_catalog_file) followed by the standard packager.
%
%   Only CERTIFIED entries are carried, so the catalog's availability grid
%   means what it says.
%
%% Inputs:
%
%  sheetMat                 char                    build_arrival_sheet output
%  ribMats                  cellstr | char | []     build_ribs outputs
%  opts                     struct (optional)
%   .thrustN [0.070] .ispS [900] .m0kg [150] .nD [12] .sD0 [0]
%   .tag ['70mN'] .outDir [results/] .name ['costate_catalog_dro_tulip_70mN']
%
%% Outputs:
%
%  cat_                     struct                  the catalog (also saved)
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, ribMats = {}; end
if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
outDir = d('outDir', fullfile(here, 'results'));
tag = d('tag', '70mN');
name = d('name', 'costate_catalog_dro_tulip_70mN');

L = load(sheetMat);  S = L.S;
if ischar(ribMats), ribMats = {ribMats}; end
ribs = {};
for k = 1:numel(ribMats)
    if ~isfile(ribMats{k}), fprintf('  (missing rib file %s -- skipped)\n', ribMats{k}); continue, end
    Rk = load(ribMats{k});
    % a rib must belong to the SHEET'S problem, not merely exist
    if isfield(Rk, 'problem') && isfield(S, 'problem')
        assertSameProblem(Rk.problem, S.problem, ribMats{k});
    else
        fprintf('  (rib %s carries no problem identity -- accepted on trust, rebuild it)\n', ribMats{k});
    end
    for m = 1:numel(Rk.R), ribs{end+1} = Rk.R(m); end %#ok<AGROW>
end
fprintf('packaging: %d arrival phases certified, %d ribs, %d rib points\n', ...
    nnz(isfinite(S.TF)), numel(ribs), sum(cellfun(@(r) numel(r.pts), ribs)));

% one sheet FILE in the packager's layout, in its own folder (the packager
% globs a directory)
% A FRESH staging directory. The packager globs this folder, so a stale
% sheet left from an earlier run would be consumed alongside the new one and
% shipped unaudited. (Astra chain review 2026-09-10.)
sheetDir = fullfile(outDir, ['phase_sheet_' tag]);
if isfolder(sheetDir), rmdir(sheetDir, 's'); end
mkdir(sheetDir);
sheetFile = fullfile(sheetDir, sprintf('dro_tulip_%s_tau1_Np7.mat', tag));
Q = sheet_to_catalog_file(S, ribs, sheetFile, opts);
fprintf('  sheet file: %d of %d grid points certified\n', nnz(Q.OK), numel(Q.OK));

cat_ = build_costate_catalog_family(sheetDir, fullfile(outDir, [name '.mat']), struct( ...
    'glob', sprintf('dro_tulip_%s_*.mat', tag), 'name', name, ...
    'description', sprintf(['DRO (tau = 1) -> 7-petal tulip minimum-time PMP costates at a ' ...
        'SECOND propulsion regime: %.0f mN, Isp %g s, m0 %g kg (the shipped DRO -> tulip ' ...
        'catalog is Isp 1710 s). One thrust rung, %d x %d departure x arrival phase grid. ' ...
        'Every entry passed the same gate stack as the shipped catalogs -- multiple-shooting ' ...
        'residual, flown arrival in position AND velocity, pumpkyn tfMin acceptance, the ' ...
        'free-time conjugate test, and the three BCT sufficiency-hypothesis gates.'], ...
        d('thrustN', 0.070)*1000, d('ispS', 900), d('m0kg', 150), d('nD', 12), numel(S.sA)), ...
    'provenance', ['Arrival axis by pseudo-arclength continuation (arclength_ms / ' ...
        'arclength_arrival), departure axis by the bisecting walker (rib_from_crossing); ' ...
        'assembled by sheet_from_arcs with the certified library as seeds; ' ...
        'FINDINGS 37-38, 2026-09-09.'], ...
    'depReconstruction', 'DRO of period tau_dep, pumpkyn get_family_orbit(''dro'', tau)'));
end

function assertSameProblem(a, b, src)
% ASSERTSAMEPROBLEM  A rib may only join a sheet certified at the same
% operating point.  INPUTS: a; b; src.
f = {'thrustN', 'ispS', 'm0kg', 'tauDRO', 'NpTulip', 'pmTulip', 'sD'};
for k = 1:numel(f)
    assert(isfield(a, f{k}) && isfield(b, f{k}) && ...
           abs(a.(f{k}) - b.(f{k})) <= 1e-12*max(abs(b.(f{k})), 1), ...
           'rib %s was certified at a different %s (%g vs %g)', src, f{k}, ...
           getfielddef(a, f{k}), getfielddef(b, f{k}));
end
end

function v = getfielddef(s, f)
% GETFIELDDEF  Field or NaN.  INPUTS: s; f.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = NaN; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
