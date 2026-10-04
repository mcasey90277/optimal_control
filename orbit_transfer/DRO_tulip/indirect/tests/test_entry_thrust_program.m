function ok = test_entry_thrust_program()
% TEST_ENTRY_THRUST_PROGRAM  Flying an entry segment by segment from its
% junctions reaches the arrival state the ORBIT gives at its phase (oracle:
% B.stateA, not the catalog), with small junction defects; the reported
% thrust direction is the one the equations of motion apply (oracle: EoM
% with thrust minus EoM without); K is read from the array (K = 12 works).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
% Fixture: a PHASE-CATALOG entry (sheet 1, sD_frac = 0, sA_frac = 0.0754), whose
% keyed phases its audit flies to within 1 km (the anchor's stored sA is rounded).
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));
c = L.costate_catalog_dro_tulip_70mN;
sh = c.sheets(1);  iD = find(sh.sD_frac == 0, 1);  iA = 2;
assert(abs(sh.sA_frac(iA) - 0.0754) < 1e-12, 'fixture: column 2 is not sA_frac 0.0754');
z8 = sh.z8(:, sh.entry_index(iD, iA));
[B, ~] = arclength_arrival('setup', catalog_setup_request(c, 0));
rv0 = B.stateD(0);  rv0 = rv0(1:6);  rvf = B.stateA(sh.sA_frac(iA));  rvf = rvf(1:6);
phys = struct('Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu, 'lStar', 389703.264829278, 'tStar', 382981.289129055, 'm0kg', 150);
% Defect tol 5e-2 km: the pchip cut of one flight measures 1.1e-2 km (K=24), 6.9e-3 km (K=12).
for K = [24 12]
    seed = seed_from_z8(z8, rv0, K, B.Tnd, B.cnd, B.mu);
    P = entry_thrust_program(z8, seed.Y(:, 1:K), rv0, rvf, phys);
    ok = chk(ok, P.K == K && P.flyKm < 1 && P.flyVms < 0.1, sprintf('K = %d: arrives at the orbit''s state (%.3f km, %.4f m/s)', K, P.flyKm, P.flyVms));
    ok = chk(ok, max(P.defectKm) < 5e-2 && P.startErr < 1e-12, sprintf('K = %d: junction defects %.1e km, start %.1e', K, max(P.defectKm), P.startErr));
end
k = round(numel(P.t)/2);
dyT = pumpkyn.cr3bp.tfMinEoM(P.t(k), P.Y(k, :).', B.Tnd, B.cnd, B.mu);
dy0 = pumpkyn.cr3bp.tfMinEoM(P.t(k), P.Y(k, :).', 0, B.cnd, B.mu);
a = dyT(4:6) - dy0(4:6);  a = a/sqrt(sum(a.^2));
ok = chk(ok, sqrt(sum((a(:) - P.u(k, :).').^2)) < 1e-9, 'thrust direction = the EoM''s applied acceleration direction');
ok = chk(ok, abs(P.massKg(end) - 150*(1 - B.Tnd*z8(8)/B.cnd)) < 1e-6, 'final mass = the all-burn mass law');
if ok, fprintf('test_entry_thrust_program: ALL PASS\n'); else, fprintf('test_entry_thrust_program: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
