function T = load_catalog_sheets(root)
%% Purpose:
%
%   The four min-time costate catalogs (DRO_tulip, HALO_tulip, DPO_tulip,
%   HALO_HALO) flattened into ONE sheet table: one row per (departure
%   orbit, arrival orbit) sheet, 56 rows, with the per-sheet phase grids
%   attached. The v1 DRO catalog is mapped onto the v2 field names here
%   (dep_family 'dro', dep_params struct('tau', tauDRO), ...) so that no
%   consumer ever branches on schema version. Catalog order is DRO_tulip,
%   HALO_tulip, DPO_tulip, HALO_HALO; sheet order is as stored. The
%   `accepted` mask is the certified witness: has_solution & conj_pass ==
%   1 & gate_min_lamv > 0 & gate_min_qmt > 0 & gate_dimS == 1 (conj_pass
%   stores -1 for "not run", which is not a pass). Every consumer counts
%   and minimizes over `accepted`; has_solution is reported beside it.
%
%   Refusals, by name: a missing catalog file (:missingCatalog, names the
%   path), rung vectors that differ between catalogs (:rungs), a sheet
%   whose arrays are not nD x nA x numel(rungs) (:shape).
%
%% Inputs:
%
%  root                     char                    The orbit_transfer
%                                                   folder; [] for the
%                                                   parent of this file
%                                                   (a bare call runs the
%                                                   demo)
%
%% Outputs:
%
%  T                        struct [56 x 1]         Fields: catalog,
%                                                   sheet_index,
%                                                   dep_family, dep_params,
%                                                   arr_family, arr_params,
%                                                   dep_key, arr_key,
%                                                   T_dep, T_arr, rungs_N
%                                                   [1 x 11], has_solution
%                                                   / conj_pass (== 1) /
%                                                   accepted (solved,
%                                                   conjugate-pass, three
%                                                   hypothesis gates) /
%                                                   tf_nd [6 x 6 x 11],
%                                                   sD_frac, sA_frac
%                                                   [1 x 6], constants
%                                                   (muStar, lStar_km,
%                                                   tStar_s)
%
%% Revision History:
%  M. Casey                                                   (c) 11/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: load, print the census of closed cells at 1 N and 0.5 N per sheet
   %(callers pass root or [] -- a bare call is this demo)
     T = load_catalog_sheets([]);
     r1 = find(T(1).rungs_N == 1);  r05 = find(T(1).rungs_N == 0.5);
     fprintf('%-11s %-18s %-14s %5s %5s\n', 'catalog', 'dep', 'arr', 'n1N', 'n05N');
     for kk = 1:numel(T)
         fprintf('%-11s %-18s %-14s %5d %5d\n', T(kk).catalog, T(kk).dep_key, T(kk).arr_key, ...
                 nnz(T(kk).has_solution(:, :, r1)), nnz(T(kk).has_solution(:, :, r05)));
     end
     return
end

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
end

catalogs = { 'dro_tulip',  fullfile(root, 'DRO_tulip',  'direct', 'results', 'costate_catalog_dro_tulip.mat')
             'halo_tulip', fullfile(root, 'HALO_tulip', 'direct', 'results', 'costate_catalog_halo_tulip.mat')
             'dpo_tulip',  fullfile(root, 'DPO_tulip',  'direct', 'results', 'costate_catalog_dpo_tulip.mat')
             'halo_halo',  fullfile(root, 'HALO_HALO',  'direct', 'results', 'costate_catalog_halo_halo.mat') };

nCat = size(catalogs, 1);
Cs   = cell(nCat, 1);
for cc = 1:nCat
    file = catalogs{cc, 2};
    if exist(file, 'file') ~= 2
        error('orbit_network:load_catalog_sheets:missingCatalog', 'catalog file not found: %s', file);
    end
    S      = load(file);
    fn     = fieldnames(S);
    Cs{cc} = S.(fn{1});
end
rungs = Cs{1}.rungs_N(:).';
for cc = 2:nCat
    if ~isequal(Cs{cc}.rungs_N(:).', rungs)
        error('orbit_network:load_catalog_sheets:rungs', '%s: rungs %s differ from %s', ...
              catalogs{cc, 1}, mat2str(Cs{cc}.rungs_N(:).'), mat2str(rungs));
    end
end
nR   = numel(rungs);
rows = cell(sum(cellfun(@(C) numel(C.sheets), Cs)), 1);
n    = 0;
for cc = 1:nCat
    C = Cs{cc};
    for ss = 1:numel(C.sheets)
        sh = C.sheets(ss);
        r.catalog     = catalogs{cc, 1};
        r.sheet_index = ss;
        if isfield(sh, 'dep_family')                       % schema v2
            r.dep_family = sh.dep_family;
            r.dep_params = sh.dep_params;
            r.arr_family = sh.arr_family;
            r.arr_params = sh.arr_params;
            r.T_dep      = sh.tau_dep;
            r.T_arr      = sh.tau_arr;
        else                                               % schema v1 (DRO catalog)
            r.dep_family = 'dro';
            r.dep_params = struct('tau', sh.tauDRO);
            r.arr_family = 'tulip';
            r.arr_params = struct('Np', sh.Np, 'pm', sh.pm);
            r.T_dep      = sh.tauDRO;
            r.T_arr      = sh.period_tulip_nd;
        end
        r.dep_key = node_key(r.dep_family, r.dep_params);
        r.arr_key = node_key(r.arr_family, r.arr_params);
        r.rungs_N = rungs;
        nD = numel(sh.sD_frac);  nA = numel(sh.sA_frac);
        for f = {'has_solution', 'tf_nd', 'conj_pass', 'gate_min_lamv', 'gate_min_qmt', 'gate_dimS'}
            v = sh.(f{1});
            if ~isequal(size(v), [nD nA nR])
                error('orbit_network:load_catalog_sheets:shape', '%s sheet %d: %s is %s, expected [%d %d %d]', ...
                      catalogs{cc, 1}, ss, f{1}, mat2str(size(v)), nD, nA, nR);
            end
            r.(f{1}) = v;
        end
        r.has_solution = logical(r.has_solution);
        r.conj_pass    = (r.conj_pass == 1);              % the catalog stores -1 for "not run": logical(-1) would be a pass
        r.accepted     = r.has_solution & r.conj_pass & (r.gate_min_lamv > 0) & (r.gate_min_qmt > 0) & (r.gate_dimS == 1);
        r = rmfield(r, {'gate_min_lamv', 'gate_min_qmt', 'gate_dimS'});
        r.sD_frac      = sh.sD_frac(:).';
        r.sA_frac      = sh.sA_frac(:).';
        r.constants    = struct('muStar', C.constants.muStar, 'lStar_km', C.constants.lStar_km, 'tStar_s', C.constants.tStar_s);
        n       = n + 1;
        rows{n} = r;
    end
end
T = vertcat(rows{:});
end

% ------------------------------------------------------------------------
function key = node_key(family, p)
%% Purpose:
%
%   Canonical node name from a family and its get_family_orbit parameter
%   struct: dro_tau1, dpo_tau2, halo_L2_tau2.2, halo_L1_tau1.8037,
%   tulip_Np7. Numbers print with %g so 1 and 1.0 agree.
%
switch family
    case 'dro',   key = sprintf('dro_tau%g', p.tau);
    case 'dpo',   key = sprintf('dpo_tau%g', p.tau);
    case 'halo',  key = sprintf('halo_L%d_tau%g', p.Lpt, p.tau);
    case 'tulip', key = sprintf('tulip_Np%d', p.Np);
    otherwise,    error('orbit_network:load_catalog_sheets:family', 'unknown family %s', family);
end
end
