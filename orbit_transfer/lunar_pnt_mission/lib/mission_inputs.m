function M = mission_inputs()
%% Purpose:
%
%   THE VERIFIED INPUTS of the lunar PNT mission study: the 70 mN DRO ->
%   tulip catalog (single sheet, 24 x 48 cells), its regime constants, the
%   two orbits as periodic phase_state handles (ONE orbit definition for
%   DRO states, slot states and transfer endpoints), and the tulip
%   phasing. Every number the later tasks use is read from here, none is
%   retyped. Verifications are recorded in M.checks and enforced: any
%   failure throws mission_inputs:check.
%
%   Checks: one sheet; all 1152 cells has_solution with status 4;
%   sh.tau_dep / sh.tau_arr equal the get_family_orbit periods (1e-12);
%   tulip_phasing tau0 / Np / pm equal tauT / sh.Np / sh.pm; cnd from
%   the kit README formula equals cat_.thruster.c_nd (1e-12
%   relative); the regime literals (muStar, lStar, tStar, 0.07 N, 900 s,
%   150 kg) equal the plan's pinned values. THE ORACLE for Tnd, cnd and the
%   endpoint definition is the flown arrival: catalog cells (3,5) and
%   (17,40), flown by fly_transfer from depState(sD) to arrState(sA), must
%   land within 100 km and 10 m/s. (The catalog's derive field holds
%   formula text only, and the flown-mass vs 1 - Tnd*tf/cnd comparison is
%   an integrator consistency check, recorded but not an oracle.)
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  M.cat       struct                  Catalog struct
%  M.sh        struct                  Its single sheet
%  M.muStar    scalar                  CR3BP mass parameter [ND]
%  M.lStar     scalar                  Length unit [km]
%  M.tStar     scalar                  Time unit [s]
%  M.thrustN   scalar                  Thrust [N]
%  M.ispS      scalar                  Specific impulse [s]
%  M.m0kg      scalar                  Initial mass [kg]
%  M.Tnd       scalar                  ND acceleration at m = 1
%  M.cnd       scalar                  ND exhaust velocity
%  M.phys      struct                  Tnd cnd muStar lStar tStar m0kg
%                                      (the fly_transfer physics struct)
%  M.tauD      scalar                  DRO period [ND]
%  M.tauT      scalar                  Tulip period [ND]
%  M.sD        1 x 24                  Catalog departure phases [fraction]
%  M.sA        1 x 48                  Catalog arrival phases [fraction]
%  M.tf        24 x 48                 Catalog min-time table [ND]
%  M.depState  function handle         DRO phase fraction -> [6 x 1]
%  M.arrState  function handle         Tulip phase fraction -> [6 x 1]
%  M.dtau6     1 x 6                   Slot orbit-time offsets [ND]
%  M.droPhase0 1 x 6                   DRO phase of spacecraft k at t = 0
%  M.slotPhase0 1 x 6                  Tulip phase of slot j at t = 0
%  M.propKgPerNd scalar                Propellant per ND of burn [kg/ND]
%  M.checks    struct array            .name .value .pass
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ot  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
kit = fullfile(ot, 'DRO_tulip', 'deliverables', ...
               'costate_library_70mN_user_deliverable9');
L    = load(fullfile(kit, 'costate_catalog_dro_tulip_70mN.mat'));
cat_ = L.costate_catalog_dro_tulip_70mN;
sh   = cat_.sheets(1);

M.cat     = cat_;
M.sh      = sh;
M.muStar  = cat_.constants.muStar;
M.lStar   = cat_.constants.lStar_km;
M.tStar   = cat_.constants.tStar_s;
M.thrustN = cat_.rungs_N(1);
M.ispS    = cat_.thruster.isp_s;
M.m0kg    = cat_.thruster.m0_kg;
M.Tnd     = (M.thrustN/M.m0kg)*M.tStar^2/(M.lStar*1000);
M.cnd     = cat_.thruster.c_nd;
M.phys    = struct('Tnd', M.Tnd, 'cnd', M.cnd, 'muStar', M.muStar, ...
                   'lStar', M.lStar, 'tStar', M.tStar, 'm0kg', M.m0kg);

dp = sh.dep_params;  dp.muStar = M.muStar;
ap = sh.arr_params;  ap.muStar = M.muStar;
[tD, rvD] = get_family_orbit(sh.dep_family, dp);
[tT, rvT] = get_family_orbit(sh.arr_family, ap);
M.tauD = tD(end);
M.tauT = tT(end);
M.depState = phase_state(tD, rvD);
M.arrState = phase_state(tT, rvT);

M.sD = sh.sD_frac(:).';
M.sA = sh.sA_frac(:).';
M.tf = sh.tf_nd;

[M.dtau6, ~, npPh, tau0, pmPh] = tulip_phasing();
M.droPhase0   = (0:5)/6;
M.slotPhase0  = M.dtau6/M.tauT;
M.propKgPerNd = M.m0kg*M.Tnd/M.cnd;

% ---- verifications ----
g0nd = 9.80665*M.tStar^2/(1000*M.lStar);
cnd2 = (M.ispS/M.tStar)*g0nd;                 % kit README formula
ck = {};
ck(end+1, :) = {'one sheet', numel(cat_.sheets), isscalar(cat_.sheets)};
ck(end+1, :) = {'all cells has_solution', nnz(sh.has_solution), ...
                nnz(sh.has_solution) == 1152 && numel(sh.has_solution) == 1152};
ck(end+1, :) = {'all cells status 4', nnz(sh.status == 4), ...
                nnz(sh.status == 4) == 1152};
ck(end+1, :) = {'tau_dep vs DRO period', abs(sh.tau_dep - M.tauD), ...
                abs(sh.tau_dep - M.tauD) < 1e-12};
ck(end+1, :) = {'tau_arr vs tulip period', abs(sh.tau_arr - M.tauT), ...
                abs(sh.tau_arr - M.tauT) < 1e-12};
ck(end+1, :) = {'phasing tau0 vs tauT', abs(tau0 - M.tauT), ...
                abs(tau0 - M.tauT) < 1e-12};
ck(end+1, :) = {'phasing Np vs sheet', npPh - sh.Np, npPh == sh.Np};
ck(end+1, :) = {'phasing pm vs sheet', pmPh - sh.pm, pmPh == sh.pm};
ck(end+1, :) = {'cnd vs README formula', abs(M.cnd - cnd2)/abs(cnd2), ...
                abs(M.cnd - cnd2)/abs(cnd2) < 1e-12};
% THE oracle for Tnd, cnd and the endpoint definition: a certified catalog
% z8 flown from depState(sD) lands on arrState(sA) only if they are right
% (a 1 percent Tnd error misses by far more than 100 km). Two well-
% separated cells.
cells = [3 5; 17 40];
for kc = 1:size(cells, 1)
    iD = cells(kc, 1);  iA = cells(kc, 2);
    Fc = fly_transfer(sh.z8(:, sh.entry_index(iD, iA)), M.depState(M.sD(iD)), ...
                      M.arrState(M.sA(iA)), M.phys);
    ck(end+1, :) = {sprintf(['catalog cell (%d,%d) flies onto the tulip ' ...
                    '(Tnd, cnd, endpoints): %.3g km %.3g m/s'], iD, iA, ...
                    Fc.flyKm, Fc.flyVms), Fc.flyKm, ...
                    Fc.flyKm < 100 && Fc.flyVms < 10};
    if kc == 1, Fm = Fc; end
end
% integrator mass-law consistency only (NOT an oracle: both sides use the
% same Tnd and cnd)
mfLaw = 1 - M.Tnd*M.tf(3, 5)/M.cnd;
ck(end+1, :) = {'integrator mass law consistency (not an oracle)', ...
                abs(Fm.mf - mfLaw), abs(Fm.mf - mfLaw) < 1e-6};
% regime: the plan's pinned values (the one place these literals live)
ck(end+1, :) = {'muStar pinned', M.muStar - 0.012150585609624, ...
                abs(M.muStar - 0.012150585609624) < 1e-15};
ck(end+1, :) = {'lStar pinned [km]', M.lStar - 389703.264829278, ...
                abs(M.lStar - 389703.264829278) < 1e-6};
ck(end+1, :) = {'tStar pinned [s]', M.tStar - 382981.289129055, ...
                abs(M.tStar - 382981.289129055) < 1e-6};
ck(end+1, :) = {'thrust pinned [N]', M.thrustN - 0.07, abs(M.thrustN - 0.07) < 1e-15};
ck(end+1, :) = {'Isp pinned [s]', M.ispS - 900, M.ispS == 900};
ck(end+1, :) = {'m0 pinned [kg]', M.m0kg - 150, M.m0kg == 150};
M.checks = struct('name', ck(:, 1), 'value', ck(:, 2), 'pass', ck(:, 3));
M.checks = reshape(M.checks, 1, []);
if ~all([M.checks.pass])
    bad = M.checks(~[M.checks.pass]);
    error('mission_inputs:check', 'failed: %s', strjoin({bad.name}, '; '));
end
if nargout == 0
    for kk = 1:numel(M.checks)
        fprintf('  %-32s %12.4g  %s\n', M.checks(kk).name, ...
                M.checks(kk).value, string(M.checks(kk).pass));
    end
end
end
