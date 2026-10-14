function ok = test_catalog_transfer()
% TEST_CATALOG_TRANSFER  catalog_transfer flies one library entry with
% pumpkyn's propagator and hands back its history: (1) a transfers_at row
% flown from its junctions arrives on the arrival orbit at the row's phase
% (oracle from the orbit table), runs 0 .. z8(8), passes through every
% stored junction, carries t_f / Delta-V equal to the row's, a unit thrust
% direction equal to -lam_v/|lam_v| from the state columns, and a mass that
% only decreases; (2) a forced z8 flight ends where the junction flight does
% and is the method recorded; (3) the struct('iD', 'iA') form gives the same
% flight as the cell's primary row; (4) an extremes row (dep_frac / arr_frac)
% is accepted, and an extremes row's OWN thrust and sheet are used when opts
% does not name them (a row re-labelled 0.05 N flies at 0.05 N; on a
% two-sheet copy of the catalog a row naming sheet 2 gets sheet 2's orbits,
% a row naming no sheet is refused, the struct(iD, iA) form is refused);
% (5) opts.plot returns a figure whose line data is the flight; (6)
% 'junctions' on a row without them refuses with its id, malformed junctions
% warn and fall back to z8; (7) the repo dependency closure is exactly
% {get_family_orbit, status_key, transfers_at}. Checks marked IDENTITY are
% regression guards that compare the code's own formulas, not oracles.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
repoRoot = fileparts(fileparts(fileparts(here)));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
[prim, alt2] = sample_status_pair(c);
sh = c.sheets(1);  lStar = c.constants.lStar_km;  tStar = c.constants.tStar_s;
pA = sh.arr_params;  pA.muStar = c.constants.muStar;
[tT, rvT] = get_family_orbit(sh.arr_family, pA);

%% (1) junction flight of a transfers_at row
F = catalog_transfer(c, alt2);
rvf = interp1(tT/tT(end), rvT, mod(alt2.sA, 1), 'spline');
missKm = sqrt(sum((F.Y(end, 1:3) - rvf(1:3)).^2))*lStar;
ok = chk(ok, strcmp(F.method, 'junctions'), 'a row with junctions is flown from them (auto)');
ok = chk(ok, missKm < 1, sprintf('arrives on the arrival orbit at sA = %.4f (miss %.3f km)', alt2.sA, missKm));
ok = chk(ok, abs(F.flyKm - missKm) < 1e-9, 'IDENTITY: F.flyKm is that miss');
ok = chk(ok, F.t(1) == 0 && issorted(F.t) && abs(F.t(end) - alt2.z8(8)) < 1e-12, 'time base runs 0 .. z8(8), sorted');
J = alt2.junctions;  K = size(J, 2);  h = alt2.z8(8)/K;  worst = 0;
for q = 2:K
    Yq = interp1(F.t, F.Y, (q - 1)*h);
    worst = max(worst, sqrt(sum((Yq(1:3) - J(1:3, q).').^2))*lStar);
end
ok = chk(ok, worst < 1e-3, sprintf('passes through its %d stored junctions (worst %.2e km)', K, worst));
ok = chk(ok, abs(F.tfDays - alt2.tfDays) < 1e-9 && abs(F.dvKms - alt2.dvKms) < 1e-6 && abs(F.propellantKg - alt2.propellantKg) < 1e-6, ...
    sprintf('t_f, Delta-V and propellant equal the row''s (%.3f d, %.4f km/s, %.3f kg)', F.tfDays, F.dvKms, F.propellantKg));
lv = F.Y(:, 11:13);  uRef = -lv ./ sqrt(sum(lv.^2, 2));
ok = chk(ok, max(abs(sqrt(sum(F.u.^2, 2)) - 1)) < 1e-12 && max(abs(F.u(:) - uRef(:))) < 1e-12, 'IDENTITY: u is the unit vector -lam_v/|lam_v| from the state columns');
ok = chk(ok, all(diff(F.massKg) <= 0) && abs(F.massKg(1) - c.thruster.m0_kg) < 1e-9 && abs(F.finalMassKg - F.massKg(end)) < 1e-12, ...
    'mass starts at m0 and only decreases; finalMassKg is the last sample');
ok = chk(ok, F.status == 2 && strcmp(F.statusName, alt2.statusName), 'status and its name are carried from the row');
ok = chk(ok, max(abs(F.rKm(:) - reshape(F.Y(:, 1:3)*lStar, [], 1))) < 1e-9 && ...
    max(abs(F.vKms(:) - reshape(F.Y(:, 4:6)*lStar/tStar, [], 1))) < 1e-12 && ...
    max(abs(F.tDays - F.t*tStar/86400)) < 1e-12, 'IDENTITY: unit conversions ride along (km, km/s, days)');

%% (2) forced z8 flight
Fz = catalog_transfer(c, prim, struct('method', 'z8'));
Fj = catalog_transfer(c, prim);
dEnd = sqrt(sum((Fz.Y(end, 1:3) - Fj.Y(end, 1:3)).^2))*lStar;
ok = chk(ok, strcmp(Fz.method, 'z8') && dEnd < 1, sprintf('a forced z8 flight ends within %.3f km of the junction flight', dEnd));

%% (3) the cell form
T = transfers_at(c, sh.sD_frac(1), sh.sA_frac(1));
ip = find(strcmp({T.kind}, 'primary'), 1);
Fc = catalog_transfer(c, struct('iD', 1, 'iA', 1));
Fr = catalog_transfer(c, T(ip));
ok = chk(ok, isequal(Fc.Y, Fr.Y) && isequal(Fc.t, Fr.t), 'struct(iD, iA) flies the cell''s primary row');

%% (4) an extremes row, and its own thrust and sheet
[eMin, ~] = costate_catalog_extremes(c, struct(), false, 'time');
Fe = catalog_transfer(c, eMin);
ok = chk(ok, abs(Fe.sD - eMin.dep_frac) < 1e-15 && abs(Fe.sA - eMin.arr_frac) < 1e-15 && Fe.flyKm < 1, ...
    sprintf('an extremes row (dep_frac / arr_frac) is accepted with no opts and arrives (miss %.3f km)', Fe.flyKm));
eLow = eMin;  eLow.thrustN = 0.05;                 % re-labelled: the flight must use 0.05 N, not rungs_N(1)
Fl = catalog_transfer(c, eLow, struct('method', 'z8'));   % one flight from m0: the rocket equation holds exactly
TndLow = (0.05/c.thruster.m0_kg)*tStar^2/(lStar*1000);
mfLow  = 1 - TndLow*eLow.z8(8)/c.thruster.c_nd;
ok = chk(ok, abs(Fl.finalMassKg - mfLow*c.thruster.m0_kg) < 1e-6 && abs(Fl.finalMassKg - Fe.finalMassKg) > 0.1, ...
    sprintf('an extremes row''s own thrustN is honored (final mass %.3f kg at 0.05 N vs %.3f kg at the rung)', Fl.finalMassKg, Fe.finalMassKg));
c2 = c;  c2.sheets(2) = c.sheets(1);  c2.sheets(2).tauDRO = 1.1;  c2.sheets(2).dep_params.tau = 1.1;  c2.sheets(2).tau_dep = 1.1;
eS2 = eMin;  eS2.tauDRO = 1.1;
F2 = catalog_transfer(c2, eS2);
ok = chk(ok, abs(F2.orbits.tD(end) - 1.1) < 1e-6, sprintf('on a two-sheet catalog a row naming tauDRO 1.1 gets sheet 2''s orbit (period %.4f ND)', F2.orbits.tD(end)));
refusedS = false;
try, catalog_transfer(c2, prim); catch ME, refusedS = strcmp(ME.identifier, 'catalog_transfer:sheet'); end
ok = chk(ok, refusedS, 'on a two-sheet catalog a row naming no sheet is refused (catalog_transfer:sheet)');
refusedC = false;
try, catalog_transfer(c2, struct('iD', 1, 'iA', 1)); catch ME, refusedC = strcmp(ME.identifier, 'catalog_transfer:multiSheet'); end
ok = chk(ok, refusedC, 'on a two-sheet catalog the struct(iD, iA) form is refused (catalog_transfer:multiSheet)');
F2s = catalog_transfer(c2, prim, struct('sheet', 1));
ok = chk(ok, isequal(F2s.Y, Fj.Y), 'opts.sheet = 1 on the two-sheet catalog reproduces the one-sheet flight');

%% (5) the plot
[Fp, hFig] = catalog_transfer(c, prim, struct('plot', true, 'visible', false));
isFig = ~isempty(hFig) && isgraphics(hFig, 'figure');
lineOk = false;
if isFig
    ln = findobj(hFig, 'Type', 'line', 'LineStyle', '-', 'Color', 'b');
    lineOk = ~isempty(ln) && isequal(ln(1).XData(:), Fp.Y(:, 1)) && isequal(ln(1).ZData(:), Fp.Y(:, 3));
    close(hFig);
end
ok = chk(ok, isFig && lineOk, 'opts.plot returns a figure whose transfer line is the flight');
[~, hNone] = catalog_transfer(c, prim);
ok = chk(ok, isempty(hNone), 'no figure without opts.plot');

%% (6) refusal
bare = struct('z8', prim.z8, 'sD', prim.sD, 'sA', prim.sA);
refused = false;
try, catalog_transfer(c, bare, struct('method', 'junctions')); catch ME, refused = strcmp(ME.identifier, 'catalog_transfer:noJunctions'); end
ok = chk(ok, refused, 'method junctions on a bare z8 entry refuses with catalog_transfer:noJunctions');
Fb = catalog_transfer(c, bare);
ok = chk(ok, strcmp(Fb.method, 'z8') && isnan(Fb.status) && Fb.flyKm < 1, 'a bare z8 entry flies from z8 (auto), status not annotated');
badJ = prim;  badJ.junctions = prim.junctions(1:13, :);
lastwarn('');  Fw = catalog_transfer(c, badJ);  [~, wid] = lastwarn;
ok = chk(ok, strcmp(Fw.method, 'z8') && strcmp(wid, 'catalog_transfer:junctionsShape'), 'malformed junctions warn (catalog_transfer:junctionsShape) and fall back to z8');

%% (7) closure
req = matlab.codetools.requiredFilesAndProducts(fullfile(here, 'catalog_transfer.m'));
inRepo = req(startsWith(req, repoRoot));
[~, names] = cellfun(@fileparts, inRepo, 'UniformOutput', false);
names = setdiff(names, 'catalog_transfer');
want = {'get_family_orbit', 'status_key', 'transfers_at'};
ok = chk(ok, isequal(sort(names(:)), sort(want(:))), sprintf('repo dependencies are exactly {%s} (got {%s})', strjoin(want, ', '), strjoin(names, ', ')));
if ok, fprintf('test_catalog_transfer: ALL PASS\n'); else, fprintf('test_catalog_transfer: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
