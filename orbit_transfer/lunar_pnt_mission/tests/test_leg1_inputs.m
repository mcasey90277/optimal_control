function [ok, skipped] = test_leg1_inputs()
%% Purpose:
%
%   Verifies leg1_inputs (leg 1, GTO -> DRO): regime values read from the
%   package (5 N, 1710 s, 150 kg), DRO period 1 ND, GTO period against the
%   Kepler period of the default elements, the clean-cell table (34 on the
%   0 deg sheet of the pilot package; the count is skipped with a printed
%   note when the package identity is not the pilot's), the catalog's own
%   departure endpoint reproduced by gtoState, the SIGN of the orientation
%   drift measured independently with pumpkyn fromPCI at two epochs (and
%   its deliberate failure with the opposite sign), one clean entry
%   re-flown onto its slot (< 1 km), the path setup (earth_clearance
%   resolves to GTO_DRO/lib), and the recorded checks.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%  skipped                  logical                 true if the pilot-only
%                                                   34-count was skipped
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_lunar_pnt_paths.m'));
ot = fileparts(fileparts(here));
assert(strcmp(fileparts(which('earth_clearance')), fullfile(ot, 'GTO_DRO', 'lib')), ...
       'earth_clearance not from GTO_DRO/lib');
assert(strcmp(fileparts(which('leg1_inputs')), fullfile(fileparts(here), 'lib')), ...
       'leg1_inputs not from this campaign');

L = leg1_inputs(struct());
nS = numel(L.sheetOrient);

% regime, read from the package
assert(L.thrustN == 5 && L.ispS == 1710 && L.m0kg == 150, 'regime');
assert(abs(L.phys.Tnd - (5/150)*L.tStar^2/(L.lStar*1000)) < 1e-12*L.phys.Tnd, 'Tnd');
assert(abs(L.tauD - 1) < 1e-9, 'DRO period');
assert(isequal(size(L.tf, 1:3), [12 6 nS]) && isequal(size(L.clean, 1:3), [12 6 nS]) && ...
       isequal(size(L.z8, 1:4), [8 12 6 nS]), 'table shapes');
assert(all(isnan(L.tf(~L.clean))) && all(isfinite(L.tf(L.clean))), 'tf NaN exactly off the clean cells');
assert(all([L.checks.pass]), 'input checks');

% GTO period = Kepler period of the default elements
muE  = 6.67384e-20*(1 - L.muStar)*(5.9736E24 + 7.35E22);
smaKm = (6378 + 350 + 6378 + 35786)/2;
tKep = 2*pi*sqrt(smaKm^3/muE)/L.tStar;
assert(abs(L.tauG - tKep) < 1e-6, 'tauG %.9g vs Kepler %.9g', L.tauG, tKep);
fprintf('  tauG %.9g  Kepler %.9g  tauD %.12g  sgn %+d\n', L.tauG, tKep, L.tauD, L.sgn);

% clean-cell count on the pilot (0 deg sheet), skipped otherwise
skipped = false;
isPilot = L.identity.n_entries == 271 && isequal(L.identity.orientations(:).', 0);
if isPilot
    assert(nnz(L.clean(:, :, 1)) == 34, 'pilot clean count %d', nnz(L.clean(:, :, 1)));
    fprintf('  pilot package: %d clean cells at 5 N on the 0 deg sheet\n', nnz(L.clean(:, :, 1)));
else
    skipped = true;
    fprintf('  NOTE: package is not the pilot (n_entries %d); 34-count skipped; clean per sheet: %s\n', ...
            L.identity.n_entries, mat2str(squeeze(sum(sum(L.clean, 1), 2)).'));
end
assert(nnz(L.clean(:, :, 1)) > 0, 'no clean cell on the first sheet');

% gtoState equals the catalog's own departure endpoint
S  = load(L.identity.file);
fn = fieldnames(S);  cat_ = S.(fn{1});
shs = cat_.sheets;
k1 = find(arrayfun(@(s) abs(s.arr_params.tau - 1) < 1e-9, shs), 1);
dp = shs(k1).dep_params;  dp.muStar = L.muStar;
[tG, rvG] = get_family_orbit('gto', dp);
depCat = phase_state(tG, rvG);
for iD = [1 5 12]
    d = max(abs(L.gtoState(L.sheetOrient(1), L.sD(iD)) - depCat(L.sD(iD))));
    assert(d < 1e-12, 'gtoState vs catalog endpoint at sD(%d): %.3g', iD, d);
end
assert(isequal(size(L.gtoState(0, 0.3)), [6 1]) && isequal(size(L.slotState(0.3)), [6 1]), 'handles');

% the SIGN, measured here independently of leg1_inputs: a fixed inertial
% state converted with fromPCI at 0 and dt is the GTO locus at the drifted
% orientation, for sgn as stored and NOT for -sgn
o0 = L.sheetOrient(1);  s0 = 0.3;  dt = 0.2;
xi  = pumpkyn.cr3bp.toPCI(0, L.gtoState(o0, s0).', L.muStar, L.tStar, L.lStar, 1);
x1  = pumpkyn.cr3bp.fromPCI(dt*L.tStar, xi, L.muStar, L.tStar, L.lStar, 1).';
x0  = pumpkyn.cr3bp.fromPCI(0,          xi, L.muStar, L.tStar, L.lStar, 1).';
assert(max(abs(x0 - L.gtoState(o0, s0))) < 1e-9, 'fromPCI at t = 0 does not return the state');
dOk  = max(abs(x1 - L.gtoState(o0 + L.sgn*(180/pi)*dt,  s0)));
dBad = max(abs(x1 - L.gtoState(o0 - L.sgn*(180/pi)*dt, s0)));
fprintf('  sign test: residual %.3g with sgn, %.3g with -sgn\n', dOk, dBad);
assert(dOk < 1e-9, 'stored sgn fails the fromPCI comparison (%.3g)', dOk);
assert(dBad > 1e-3, 'the opposite sign also passes (%.3g): the test cannot tell the sign', dBad);

% one clean entry re-flown onto its slot
idx = find(L.clean(:, :, 1), 1);
[iD, iA] = ind2sub([12 6], idx);
F = fly_transfer(L.z8(:, iD, iA, 1), L.gtoState(L.sheetOrient(1), L.sD(iD)), ...
                 L.slotState(L.sA(iA)), L.phys);
fprintf('  clean cell (%d,%d): flyKm %.4g  flyVms %.4g  tf %.6f ND\n', iD, iA, F.flyKm, F.flyVms, L.tf(iD, iA, 1));
assert(F.flyKm < 1, 'clean cell flown miss %.3g km', F.flyKm);
assert(abs(F.tfNd - L.tf(iD, iA, 1)) < 1e-12, 'tf table vs z8(8)');

% deliberate failure: the flight with the wrong slot phase must miss
Fb = fly_transfer(L.z8(:, iD, iA, 1), L.gtoState(L.sheetOrient(1), L.sD(iD)), ...
                  L.slotState(L.sA(iA) + 0.05), L.phys);
assert(Fb.flyKm > 1000, 'a wrong slot phase passes the miss gate (%.3g km)', Fb.flyKm);

ok = true;  fprintf('test_leg1_inputs: PASS\n');
end
