function cat_ = build_costate_catalog_family(catDir, outMat, spec)
%% Purpose:
%
%   FAMILY-AGNOSTIC catalog packager: packages any campaign's thrust-ladder
%   sheets into ONE shareable costate CATALOG in the COMPACT format (data
%   minimization per D. Koblick): only canonical nondimensional quantities
%   are stored -- phase fractions, the z8 vectors (which already contain
%   t_f), thrust rungs, and per-sheet availability/flight-time lookup
%   grids. Everything else (days, Delta-V, masses) is DERIVABLE and the
%   formulas ride along in .derive.
%
%   Every entry is a converged minimum-time PMP solution that passed all
%   three gates: multiple-shooting residual, end-to-end PMP flight to the
%   target, and acceptance UNCHANGED by pumpkyn.cr3bp.tfMin.
%
%   This is the second-campaign generalization of DRO_tulip's
%   build_costate_catalog (migration rule: reused-twice code moves to
%   costate_common). The sheet schema is IDENTICAL to the DRO catalog's --
%   sheets(k).tauDRO holds the DEPARTURE-ORBIT PERIOD whatever the family
%   (legacy field name, kept so costate_catalog_pick and friends work on
%   every catalog unchanged) -- plus two additions: dep_family and
%   dep_params, which carry the full get_family_orbit reconstruction
%   recipe (a halo needs Lpt and pm, not just its period).
%
%% Inputs:
%
%  catDir                   char                    Folder of sheet files
%
%  outMat                   char                    Output .mat
%
%  spec                     struct
%   .glob                   char                    Sheet file pattern
%                                                   (e.g. 'halo_tau*_Np*.mat')
%   .name                   char                    Catalog/variable name
%                                                   (e.g.
%                                                   'costate_catalog_halo_tulip')
%   .description            char                    One-paragraph statement
%   .provenance             char                    How it was made
%   .depReconstruction      char                    Derive-string for the
%                                                   departure orbit
%
%% Outputs:
%
%  cat_                     struct                  The catalog (also saved
%                                                   to outMat under
%                                                   spec.name)
%
%% Revision History:
%  M. Casey                                                   (c) 08/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

files = dir(fullfile(catDir, spec.glob));
files = files(~contains({files.name}, 'pre_'));
assert(~isempty(files), 'no sheet files matching %s in %s', spec.glob, catDir);

cat_ = struct();
cat_.name = spec.name;
cat_.description = spec.description;
cat_.created = datestr(now, 'yyyy-mm-dd');
cat_.provenance = spec.provenance;
cat_.schema = catalog_schema('version');       % versioned since v2

nTot = 0;
sheets = struct([]);
nS = 0;
for kf = 1:numel(files)
    Q = load(fullfile(files(kf).folder, files(kf).name));
    if nnz(Q.OK) == 0, continue, end     % an empty sheet would crash lookups
    nS = nS + 1;
    ob = Q.meta;
    if nS == 1
        cat_.constants = struct('muStar', ob.muStar, 'lStar_km', ob.lStar, ...
                                'tStar_s', ob.tStar);
        g0 = 9.80665*ob.tStar^2/(1000*ob.lStar);
        cat_.thruster = struct('isp_s', ob.ispS, 'm0_kg', ob.m0kg, ...
            'c_nd', (ob.ispS/ob.tStar)*g0, ...
            'thrust_nd_formula', 'Tmax_nd = (T_N/m0_kg)*tStar_s^2/(lStar_km*1000)');
        cat_.rungs_N = Q.rungs(:)';
    end
    assert(isequal(Q.rungs(:)', cat_.rungs_N), 'sheet %s: rung mismatch', files(kf).name);
    assert(ob.ispS == cat_.thruster.isp_s && ob.m0kg == cat_.thruster.m0_kg ...
        && abs(ob.muStar - cat_.constants.muStar) < 1e-15, ...
        'sheet %s: thruster/constants differ from the catalog''s', files(kf).name);
    [nD, nA, nR] = size(Q.OK);
    z8 = zeros(8, nnz(Q.OK));
    idx = zeros(nD, nA, nR);
    n = 0;
    for iD = 1:nD
        for iA = 1:nA
            for kr = 1:nR
                if ~Q.OK(iD,iA,kr), continue, end
                n = n + 1;
                z8(:,n) = Q.Z8(:,iD,iA,kr);
                idx(iD,iA,kr) = n;
            end
        end
    end
    tf = Q.TF;  tf(~Q.OK) = NaN;
    % v2 canonical keys + full recipes BOTH ends; legacy aliases kept so
    % every v1 consumer works on v2 tulip catalogs (tauDRO = tau_dep; Np
    % is the tulip petal count, NaN for non-tulip arrivals):
    isTulipArr = strcmpi(ob.arrFamily, 'tulip');
    % KEYS are the REQUESTED family parameters (what a user asks for),
    % not the propagated periods (those are derivable via the recipes):
    sheets(nS,1).tau_dep         = ob.tauDRO;     % departure period key (ND)
    if isfield(ob.arrParams, 'tau')
        sheets(nS,1).tau_arr     = ob.arrParams.tau;   % requested (halo etc.)
    else
        sheets(nS,1).tau_arr     = ob.periodTulip;     % tulip: locked period
    end
    sheets(nS,1).tauDRO          = ob.tauDRO;     % LEGACY alias of tau_dep
    sheets(nS,1).dep_family      = ob.depFamily;
    sheets(nS,1).dep_params      = ob.depParams;  % full reconstruction recipe
    sheets(nS,1).arr_family      = ob.arrFamily;
    sheets(nS,1).arr_params      = ob.arrParams;
    if isTulipArr
        sheets(nS,1).Np          = ob.NpTulip;
        sheets(nS,1).pm          = ob.pmTulip;
    else
        sheets(nS,1).Np          = NaN;           % petal count: tulip-only
        sheets(nS,1).pm          = ob.arrParams.pm;
    end
    sheets(nS,1).period_tulip_nd = ob.periodTulip;% LEGACY alias of tau_arr
    sheets(nS,1).sD_frac         = Q.sD(:)';
    sheets(nS,1).sA_frac         = Q.sA(:)';
    sheets(nS,1).has_solution    = Q.OK;
    sheets(nS,1).tf_nd           = tf;
    sheets(nS,1).entry_index     = idx;
    sheets(nS,1).z8              = z8;
    nTot = nTot + n;
end
% order sheets by (departure period, arrival period) for human readability
[~, ord] = sortrows([[sheets.tau_dep]', [sheets.tau_arr]']);
cat_.sheets = sheets(ord);
cat_.n_entries = nTot;

cat_.derive = struct( ...
    'days_from_nd',   't_days = t_nd * constants.tStar_s / 86400', ...
    'dep_period_nd',  'sheets(k).tauDRO (LEGACY NAME: the departure-orbit period, any family)', ...
    'phase_days',     'frac * period_nd * tStar_s/86400 (period_nd: tauDRO or period_tulip_nd)', ...
    'thrust_nd',      'see thruster.thrust_nd_formula', ...
    'm_final',        'mf = 1 - Tmax_nd*tf_nd/c_nd   (all-burn minimum time)', ...
    'deltaV_kms',     'dV = c_nd*log(1/mf) * lStar_km/tStar_s', ...
    'orbit_reconstruction', spec.depReconstruction, ...
    'pm_plus_one',    'z-mirror the pm=-1 orbit and costates (flip all z components)');

cat_.usage = sprintf(['%% [tf_nd, z8, info] = costate_catalog_pick(cat, tauDep, Np, ...\n', ...
    '%%                       depDays, arrDays, thrustN);\n', ...
    '%% info reports exactly which sheet/pair/rung you were actually given.\n']);

problems = catalog_schema('validate', cat_);
assert(isempty(problems), 'catalog fails schema validation:\n%s', ...
       strjoin(problems, sprintf('\n')));
S = struct(spec.name, cat_);
save(outMat, '-struct', 'S');
fprintf('%s: %d sheets, %d entries -> %s\n', spec.name, numel(cat_.sheets), nTot, outMat);
end
