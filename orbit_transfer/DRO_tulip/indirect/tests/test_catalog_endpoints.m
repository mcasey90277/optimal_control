function ok = test_catalog_endpoints()
% TEST_CATALOG_ENDPOINTS  The shipped endpoint rule is the certified one: on
% the library of record its departure and arrival closures equal
% arclength_arrival's (via catalog_setup_request) bitwise at every grid
% phase and at the wrap, and Tnd / cnd agree to 1e-15.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
[B, orb] = catalog_endpoints(c);
[Bref, ~] = arclength_arrival('setup', catalog_setup_request(c, 0));
sh = c.sheets(1);
eqD = true;  eqA = true;
for s = sh.sD_frac, r = B.stateD(s);  q = Bref.stateD(s);  eqD = eqD && isequal(r(1:6), q(1:6)); end
for s = sh.sA_frac, r = B.stateA(s);  q = Bref.stateA(s);  eqA = eqA && isequal(r(1:6), q(1:6)); end
ok = chk(ok, eqD && eqA, 'departure and arrival states bitwise-equal to the certified rule at every grid phase');
r = B.stateD(1 - 1e-12);  q = B.stateD(0);
gap = sqrt(sum((r(1:6) - q(1:6)).^2));
ok = chk(ok, gap < 1e-8, sprintf('REVIEW FOCUS 2: the wrap s = 1 - 1e-12 gives the s = 0 state (gap %.2e ND = the orbit table''s closure error)', gap));
ok = chk(ok, abs(B.Tnd - Bref.Tnd) <= 1e-15*Bref.Tnd && abs(B.cnd - Bref.cnd) <= 1e-15*Bref.cnd && B.mu == Bref.mu, 'Tnd, cnd, mu agree');
ok = chk(ok, isequal(size(B.stateA(0.3)), [6 1]) && all(isfield(B.problem, {'lStar', 'tStar', 'muStar', 'Tnd', 'cnd', 'thrustN', 'NpTulip', 'tauDRO'})), 'closures return 6 x 1; B.problem has the consumers'' fields');
ok = chk(ok, size(orb.rvD, 2) >= 6 && orb.tD(1) == 0 && size(orb.rvT, 2) >= 6, 'orbit samples returned for plotting');
B2 = catalog_endpoints(c, 1, 2*c.rungs_N(1));
ok = chk(ok, abs(B2.Tnd - 2*B.Tnd) <= 1e-15*B2.Tnd && B2.problem.thrustN == 2*c.rungs_N(1) && B.problem.thrustN == c.rungs_N(1), 'optional thrustN: Tnd doubles with the thrust, B.problem.thrustN carried, default unchanged');
if ok, fprintf('test_catalog_endpoints: ALL PASS\n'); else, fprintf('test_catalog_endpoints: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
