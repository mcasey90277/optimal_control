function cat_ = build_costate_catalog(catDir, outMat)
%% Purpose:
%
%   Packages the multi-orbit thrust-ladder sheets into ONE shareable costate
%   CATALOG: every (DRO period x tulip petal count) sheet of the coarse
%   sweep, in the COMPACT format (data minimization per D. Koblick): only
%   canonical nondimensional quantities are stored -- phase fractions, the
%   z8 vectors (which already contain t_f), thrust rungs, and per-sheet
%   availability/flight-time lookup grids. Everything else (days, Delta-V,
%   masses) is DERIVABLE and the formulas ride along in .derive.
%
%   Every entry is a converged minimum-time PMP solution that passed all
%   three gates: multiple-shooting residual, end-to-end PMP flight to the
%   target, and acceptance UNCHANGED by pumpkyn.cr3bp.tfMin.
%
%% Inputs:
%
%  catDir                   char                    Folder of sheet files
%                                                   ladder_tau*_Np*.mat
%
%  outMat                   char                    Output .mat holding the
%                                                   structure
%                                                   costate_catalog_dro_tulip
%
%% Outputs:
%
%  cat_                     struct
%   .constants              muStar, lStar_km, tStar_s      (stored ONCE)
%   .thruster               isp_s, m0_kg, c_nd, thrust_nd_formula
%   .rungs_N                [1 x nR] thrust rungs shared by all sheets
%   .sheets(k)              tauDRO (= DRO period, ND), Np, pm,
%                           period_tulip_nd, sD_frac, sA_frac,
%                           has_solution [nD x nA x nR], tf_nd [same],
%                           entry_index [same], z8 [8 x nEntries]
%   .derive                 formula strings for every derivable quantity
%   .n_entries              total across sheets
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

files = dir(fullfile(catDir, 'ladder_tau*_Np*.mat'));
files = files(~contains({files.name}, 'pre_'));
assert(~isempty(files), 'no sheet files in %s', catDir);

cat_ = struct();
cat_.name = 'costate_catalog_dro_tulip';
cat_.description = ['Minimum-time low-thrust DRO->Tulip costate CATALOG: ', ...
    'one phasing-torus sheet per (DRO period x tulip petal count), thrust ', ...
    'rungs walked 15->1 N. COMPACT format: canonical ND quantities only; ', ...
    'derived values reconstructed via .derive. Each entry z8 = [lambda_r(3); ', ...
    'lambda_v(3); lambda_m; tf] is a converged PMP solution accepted ', ...
    'unchanged by pumpkyn.cr3bp.tfMin.'];
cat_.created = datestr(now, 'yyyy-mm-dd');
cat_.provenance = ['Coarse catalog sweep (run_catalog_sweep): per sheet, ', ...
    '6x6 phasing torus anchored cold at 15 N and walked down; three gates ', ...
    'per entry (ms residual, flown arrival, tfMin acceptance). pm = -1 ', ...
    'branch only: pm = +1 tulips are exact z-mirrors (CR3BP symmetry), so ', ...
    'mirror the costates for that branch. Casey/Koblick 2026-08.'];

nTot = 0;
sheets = struct([]);
for kf = 1:numel(files)
    Q = load(fullfile(files(kf).folder, files(kf).name));
    ob = Q.meta;
    if kf == 1
        cat_.constants = struct('muStar', ob.muStar, 'lStar_km', ob.lStar, ...
                                'tStar_s', ob.tStar);
        g0 = 9.80665*ob.tStar^2/(1000*ob.lStar);
        cat_.thruster = struct('isp_s', ob.ispS, 'm0_kg', ob.m0kg, ...
            'c_nd', (ob.ispS/ob.tStar)*g0, ...
            'thrust_nd_formula', 'Tmax_nd = (T_N/m0_kg)*tStar_s^2/(lStar_km*1000)');
        cat_.rungs_N = Q.rungs(:)';
    end
    assert(isequal(Q.rungs(:)', cat_.rungs_N), 'sheet %s: rung mismatch', files(kf).name);
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
    sheets(kf,1).tauDRO          = ob.tauDRO;     % = DRO period (ND)
    sheets(kf,1).Np              = ob.NpTulip;
    sheets(kf,1).pm              = ob.pmTulip;
    sheets(kf,1).period_tulip_nd = ob.periodTulip;
    sheets(kf,1).sD_frac         = Q.sD(:)';
    sheets(kf,1).sA_frac         = Q.sA(:)';
    sheets(kf,1).has_solution    = Q.OK;
    sheets(kf,1).tf_nd           = tf;
    sheets(kf,1).entry_index     = idx;
    sheets(kf,1).z8              = z8;
    nTot = nTot + n;
end
% order sheets by (tauDRO, Np) for human readability
[~, ord] = sortrows([[sheets.tauDRO]', [sheets.Np]']);
cat_.sheets = sheets(ord);
cat_.n_entries = nTot;

cat_.derive = struct( ...
    'days_from_nd',   't_days = t_nd * constants.tStar_s / 86400', ...
    'dro_period_nd',  'sheets(k).tauDRO  (tau IS the period; getDRO selects by it)', ...
    'phase_days',     'frac * period_nd * tStar_s/86400 (period_nd: tauDRO or period_tulip_nd)', ...
    'thrust_nd',      'see thruster.thrust_nd_formula', ...
    'm_final',        'mf = 1 - Tmax_nd*tf_nd/c_nd   (all-burn minimum time)', ...
    'deltaV_kms',     'dV = c_nd*log(1/mf) * lStar_km/tStar_s', ...
    'orbit_reconstruction', ['DRO: getDRO(tauDRO)->cont_np->prop;  tulip: ', ...
        'getTulip(2*pi*(Np-2)/(Np-1), Np, pm)->cont_np->prop;  state at ', ...
        'phase f: interp1(t, rv, mod(f,1)*t(end), ''spline'')'], ...
    'pm_plus_one',    'z-mirror the pm=-1 orbit and costates (flip all z components)');

cat_.usage = sprintf(['%% [tf_nd, z8, info] = costate_catalog_pick(cat, tauDRO, Np, ...\n', ...
    '%%                       depDays, arrDays, thrustN);\n', ...
    '%% info reports exactly which sheet/pair/rung you were actually given.\n']);

costate_catalog_dro_tulip = cat_; %#ok<NASGU>
save(outMat, 'costate_catalog_dro_tulip');
fprintf('costate_catalog_dro_tulip: %d sheets, %d entries -> %s\n', ...
        numel(cat_.sheets), nTot, outMat);
end
