function ok = test_mission_inputs()
%% Purpose:
%
%   Verifies mission_inputs (catalog regime, grids, periods, slot phases,
%   state handles, cell (3,5) flown through fly_transfer from the
%   phase_state endpoints, propellant rate, recorded checks), the
%   deliberate failure of gate I1's oracle (1.01 x Tnd misses the tulip by
%   thousands of km), and the path
%   setup (re-fronting over a shadowing DRO_tulip/indirect, and the
%   setup_lunar_pnt_paths:shadow throw through the extraChecks seam). The
%   session path is restored exactly afterwards.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  I1 deliberate failure (1.01 x Tnd)               (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
M = mission_inputs();
assert(isequal(size(M.tf), [24 48]) && all(isfinite(M.tf(:))), 'tf table');
assert(abs(M.thrustN - 0.07) < 1e-15 && M.ispS == 900 && M.m0kg == 150, 'regime');
assert(abs(M.tauD - 1) < 1e-12 && abs(M.tauT - (5/6)*2*pi) < 1e-12, 'periods');
assert(abs(M.sA(1) - 0.0129) < 1e-9 && abs(M.sD(2) - 1/24) < 1e-12, 'catalog grids');
assert(max(abs(M.slotPhase0 - [0.1486 0.2188 0.4808 0.5525 0.8149 0.8846])) < 1e-4, 'slot phases');
assert(isequal(size(M.depState(0.3)), [6 1]) && isequal(size(M.arrState(0.3)), [6 1]), 'state handles');
% the catalog's own endpoints: flying cell (3,5) from depState(sD(3)) lands on arrState(sA(5))
z8 = M.sh.z8(:, M.sh.entry_index(3, 5));
F  = fly_transfer(z8, M.depState(M.sD(3)), M.arrState(M.sA(5)), M.phys);
assert(F.flyKm < 100 && F.flyVms < 10, 'cell (3,5) flown miss %.3g km', F.flyKm);
assert(abs(M.propKgPerNd*M.tf(3,5) - F.propellantKg) < 1e-6, 'propellant rate');
assert(all([M.checks.pass]), 'input checks');
fprintf('  cell (3,5): flyKm %.4g  flyVms %.4g  propKg %.6f\n', F.flyKm, F.flyVms, F.propellantKg);
% deliberate failure of gate I1's oracle: the same flight with Tnd scaled
% by 1.01 (a copy of M.phys) must miss by thousands of km, so a 1 percent
% thrust-scaling error cannot pass the < 100 km check
physBad = M.phys;  physBad.Tnd = 1.01*M.phys.Tnd;
Fb = fly_transfer(z8, M.depState(M.sD(3)), M.arrState(M.sA(5)), physBad);
fprintf('  cell (3,5) with 1.01 x Tnd: flyKm %.4g  flyVms %.4g (check limit 100 km)\n', Fb.flyKm, Fb.flyVms);
assert(Fb.flyKm > 1000, 'a 1 percent Tnd error passes the flown-arrival oracle (%.3g km)', Fb.flyKm);

% shadow experiment: put the DRO_tulip/indirect copies at the FRONT; the
% setup must re-front the kit so the kit copy wins
p0  = path;
restore = onCleanup(@() path(p0));
ot  = fileparts(fileparts(here));
kit = fullfile(ot, 'DRO_tulip', 'deliverables', 'costate_library_70mN_user_deliverable9');
shadow = fullfile(ot, 'DRO_tulip', 'indirect');
addpath(shadow, '-begin');
assert(strcmp(fileparts(which('transfers_at')), shadow), 'shadow not in effect');
setup_lunar_pnt_paths();
assert(strcmp(fileparts(which('transfers_at')), kit), 'kit transfers_at not re-fronted');
assert(strcmp(fileparts(which('catalog_transfer')), kit), 'kit catalog_transfer not re-fronted');
assert(strcmp(fileparts(which('get_family_orbit')), fullfile(ot, 'costate_common')), ...
       'get_family_orbit not from costate_common');
assert(strcmp(fileparts(which('mission_inputs')), fullfile(fileparts(ot), 'orbit_transfer', 'lunar_pnt_mission', 'lib')), ...
       'mission_inputs not from this campaign');

% the throw: an extra check demanding a wrong folder must raise :shadow
thrown = '';
try
    setup_lunar_pnt_paths({'tulip_phasing', tempdir});
catch err
    thrown = err.identifier;
    fprintf('  shadow throw message: %s\n', err.message);
end
assert(strcmp(thrown, 'setup_lunar_pnt_paths:shadow'), 'shadow throw not raised (%s)', thrown);
clear restore;
assert(strcmp(path, p0), 'path not restored exactly');
ok = true;  fprintf('test_mission_inputs: PASS\n');
end
