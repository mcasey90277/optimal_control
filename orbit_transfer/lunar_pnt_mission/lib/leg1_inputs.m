function L = leg1_inputs(opts)
%% Purpose:
%
%   THE VERIFIED INPUTS of leg 1 of the lunar PNT mission (six 5 N stages,
%   one inertial GTO -> the six DRO slots): the GTO -> DRO costate catalog
%   package, its regime constants, the clean-cell tables of every sheet
%   with DRO period 1, the GTO and DRO as periodic phase_state handles (ONE
%   orbit definition for transfer endpoints, departure states and slot
%   states), and the SIGN that relates elapsed leg-1 time to the rotating-
%   frame orientation of the inertial GTO. Every number later tasks use is
%   read from here; none is retyped (the regime literals are pinned only
%   as checks).
%
%   SHEET KEY TRAP: in this package sheets(k).tauDRO and .tau_dep hold the
%   GTO orientation IN DEGREES (the GTO sheet key rule), not the DRO
%   period. The DRO period is sheets(k).arr_params.tau (ND) and the
%   orientation is sheets(k).dep_params.orientDeg; sheets are selected by
%   arr_params.tau == 1 and carry their orientation in sheetOrient.
%
%   ORIENTATION DRIFT: an inertial GTO fixed in space has, in the rotating
%   frame at leg-1 time t [ND], the orientation orient0 + sgn*(180/pi)*t
%   deg. sgn is MEASURED at load time (never assumed): a fixed inertial
%   state is converted with pumpkyn.cr3bp.fromPCI at t = 0 and at dt and
%   compared with gtoState at the drifted orientation for sgn = +1 and
%   -1; exactly one must match (1e-9 ND) and the other must miss.
%
%   A cell is CLEAN iff has_solution && conj_pass == 1 at the rung whose
%   thrust is 5 N (the last of rungs_N). Failing any recorded check throws
%   leg1_inputs:check.
%
%   Self-demo (nargin == 0): prints the checks and the clean counts.
%
%% Inputs:
%
%  opts                     struct (optional)       .catalogMat  path of the
%                                                   package [GTO_DRO/results/
%                                                   costate_catalog_gto_dro.mat]
%                                                   .thrustN  rung thrust
%                                                   [5]
%
%% Outputs:
%
%  L.muStar, L.lStar, L.tStar  scalars             Package constants [ND, km, s]
%  L.thrustN, L.ispS, L.m0kg   scalars             Regime [N, s, kg]
%  L.phys       struct                  Tnd cnd muStar lStar tStar m0kg
%                                       (fly_transfer physics, nd_propulsion)
%  L.tauG       scalar                  GTO Kepler period [ND]
%  L.tauD       scalar                  DRO period [ND]
%  L.sgn        scalar                  +1 or -1 (see above)
%  L.sheetOrient 1 x nS                 Sheet GTO orientations [deg]
%  L.sD         1 x 12                  Departure (GTO) phases [fraction]
%  L.sA         1 x 6                   Arrival (DRO) phases [fraction]
%  L.tf         12 x 6 x nS             Min-time at 5 N [ND], NaN if not clean
%  L.clean      12 x 6 x nS logical     has_solution & conj_pass == 1
%  L.z8         8 x 12 x 6 x nS         [lambda0(7); tf], NaN if not clean
%  L.gtoState   function handle         (orientDeg, phase) -> [6 x 1]
%  L.slotState  function handle         DRO phase fraction -> [6 x 1]
%  L.identity   struct                  .file .sha256 .created .n_entries
%                                       .orientations .rung
%  L.checks     struct array            .name .value .pass
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    L = leg1_inputs(struct());
    for kk = 1:numel(L.checks)
        fprintf('  %-52s %12.4g  %s\n', L.checks(kk).name, L.checks(kk).value, ...
                string(L.checks(kk).pass));
    end
    fprintf('  clean cells per sheet: %s (orientations %s)\n', ...
            mat2str(squeeze(sum(sum(L.clean, 1), 2)).'), mat2str(L.sheetOrient));
    return;
end
if nargin < 1, opts = struct(); end
ot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
catFile = fdef(opts, 'catalogMat', fullfile(ot, 'GTO_DRO', 'results', ...
                                            'costate_catalog_gto_dro.mat'));
thrustN = fdef(opts, 'thrustN', 5);

S  = load(catFile);
fn = fieldnames(S);
if numel(fn) ~= 1
    error('leg1_inputs:package', 'expected ONE variable in %s, found %d', catFile, numel(fn));
end
cat_ = S.(fn{1});

L.muStar  = cat_.constants.muStar;
L.lStar   = cat_.constants.lStar_km;
L.tStar   = cat_.constants.tStar_s;
rung = find(abs(cat_.rungs_N - thrustN) < 1e-12);
if numel(rung) ~= 1
    error('leg1_inputs:rung', 'thrust %g N is not one rung of [%s]', thrustN, num2str(cat_.rungs_N));
end
L.thrustN = cat_.rungs_N(rung);
L.ispS    = cat_.thruster.isp_s;
L.m0kg    = cat_.thruster.m0_kg;
pr = nd_propulsion(L.thrustN, L.ispS, L.m0kg, L.lStar, L.tStar);
L.phys = struct('Tnd', pr.Tnd, 'cnd', pr.cnd, 'muStar', L.muStar, ...
                'lStar', L.lStar, 'tStar', L.tStar, 'm0kg', L.m0kg);

%% Sheets with DRO period 1 (sheets(k).tauDRO is the ORIENTATION, see above):
shs  = cat_.sheets(:).';
tauA = arrayfun(@(s) s.arr_params.tau, shs);
keep = find(abs(tauA - 1) < 1e-9);
if isempty(keep)
    error('leg1_inputs:sheets', 'no sheet with DRO period 1 in %s', catFile);
end
orient = arrayfun(@(s) s.dep_params.orientDeg, shs(keep));
[orient, ord] = sort(orient);
keep = keep(ord);
shs  = shs(keep);
nS   = numel(shs);
L.sheetOrient = orient(:).';
L.sD = shs(1).sD_frac(:).';
L.sA = shs(1).sA_frac(:).';
L.tf    = NaN(12, 6, nS);
L.clean = false(12, 6, nS);
L.z8    = NaN(8, 12, 6, nS);
gridsOk = true;
for ks = 1:nS
    sh = shs(ks);
    gridsOk = gridsOk && isequal(sh.sD_frac(:).', L.sD) && isequal(sh.sA_frac(:).', L.sA);
    cl = sh.has_solution(:, :, rung) & (sh.conj_pass(:, :, rung) == 1);
    L.clean(:, :, ks) = cl;
    tfk = sh.tf_nd(:, :, rung);
    t2  = NaN(12, 6);   t2(cl) = tfk(cl);
    L.tf(:, :, ks) = t2;
    for ia = 1:6
        for id = 1:12
            if cl(id, ia)
                L.z8(:, id, ia, ks) = sh.z8(:, sh.entry_index(id, ia, rung));
            end
        end
    end
end

%% Orbits: ONE definition each, as phase_state handles
depP = shs(1).dep_params;  depP.muStar = L.muStar;
orbCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
p0 = depP;  p0.orientDeg = L.sheetOrient(1);
[tG0, rvG0, infoG] = get_family_orbit('gto', p0);
L.tauG = tG0(end);
orbCache(okey(L.sheetOrient(1))) = phase_state(tG0, rvG0);
L.gtoState = @(orientDeg, phase) gto_at(orbCache, depP, orientDeg, phase);
arrP = shs(1).arr_params;  arrP.muStar = L.muStar;
[tD, rvD] = get_family_orbit('dro', arrP);
L.tauD = tD(end);
droH = phase_state(tD, rvD);
L.slotState = @(phase) droH(phase);

%% The sign, measured with fromPCI at two epochs:
o0 = L.sheetOrient(1);  s0 = 0.3;  dtM = 0.2;
x0 = L.gtoState(o0, s0);
xi = pumpkyn.cr3bp.toPCI(0, x0.', L.muStar, L.tStar, L.lStar, 1);
x1 = pumpkyn.cr3bp.fromPCI(dtM*L.tStar, xi, L.muStar, L.tStar, L.lStar, 1).';
rP = max(abs(x1 - L.gtoState(o0 + (180/pi)*dtM, s0)));
rM = max(abs(x1 - L.gtoState(o0 - (180/pi)*dtM, s0)));
if rP < rM, L.sgn = +1; else, L.sgn = -1; end
rSgn = min(rP, rM);  rOther = max(rP, rM);

%% Identity
md  = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(catFile, 'r');
raw = fread(fid, Inf, '*uint8');
fclose(fid);
md.update(raw);
L.identity = struct('file', catFile, ...
    'sha256', lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, [])), ...
    'created', char(string(cat_.created)), 'n_entries', cat_.n_entries, ...
    'orientations', L.sheetOrient, 'rung', rung);

%% Checks
ck = {};
ck(end+1, :) = {'thrust rung (5 N)', L.thrustN, L.thrustN == 5};
ck(end+1, :) = {'Isp pinned [s]', L.ispS - 1710, L.ispS == 1710};
ck(end+1, :) = {'m0 pinned [kg]', L.m0kg - 150, L.m0kg == 150};
ck(end+1, :) = {'muStar pinned', L.muStar - 0.012150585609624, ...
                abs(L.muStar - 0.012150585609624) < 1e-15};
ck(end+1, :) = {'lStar pinned [km]', L.lStar - 389703.264829278, ...
                abs(L.lStar - 389703.264829278) < 1e-6};
ck(end+1, :) = {'tStar pinned [s]', L.tStar - 382981.289129055, ...
                abs(L.tStar - 382981.289129055) < 1e-6};
ck(end+1, :) = {'cnd vs package c_nd', abs(L.phys.cnd - cat_.thruster.c_nd), ...
                abs(L.phys.cnd - cat_.thruster.c_nd) < 1e-9*cat_.thruster.c_nd};
ck(end+1, :) = {'DRO period is 1 ND', abs(L.tauD - 1), abs(L.tauD - 1) < 1e-9};
muE  = 6.67384e-20*(1 - L.muStar)*(5.9736E24 + 7.35E22);
tKep = 2*pi*sqrt(infoG.sma_km^3/muE)/L.tStar;
ck(end+1, :) = {'GTO period vs Kepler', abs(L.tauG - tKep), abs(L.tauG - tKep) < 1e-6};
ck(end+1, :) = {'sheets share the phase grids', double(gridsOk), gridsOk};
ck(end+1, :) = {'sign matches fromPCI (1e-9) ...', rSgn, rSgn < 1e-9};
ck(end+1, :) = {'... and the other sign does not (> 1e-3)', rOther, rOther > 1e-3};
% the catalog's own endpoints fly: first clean cell of the first sheet
idx = find(L.clean(:, :, 1), 1);
if isempty(idx)
    ck(end+1, :) = {'a clean cell flies onto its slot (no clean cell)', NaN, false};
else
    [id, ia] = ind2sub([12 6], idx);
    Fc = fly_transfer(L.z8(:, id, ia, 1), L.gtoState(L.sheetOrient(1), L.sD(id)), ...
                      L.slotState(L.sA(ia)), L.phys);
    ck(end+1, :) = {sprintf('clean cell (%d,%d) flies onto its slot: %.3g km %.3g m/s', ...
                    id, ia, Fc.flyKm, Fc.flyVms), Fc.flyKm, Fc.flyKm < 1 && Fc.flyVms < 1};
end
L.checks = struct('name', ck(:, 1), 'value', ck(:, 2), 'pass', ck(:, 3));
L.checks = reshape(L.checks, 1, []);
if ~all([L.checks.pass])
    bad = L.checks(~[L.checks.pass]);
    error('leg1_inputs:check', 'failed: %s', strjoin({bad.name}, '; '));
end
end

% ------------------------------------------------------------------------
function x = gto_at(cache, depP, orientDeg, phase)
% GTO_AT  State on the GTO locus of orientation orientDeg at a phase
%   fraction; loci are built once per orientation through
%   get_family_orbit and cached.  INPUTS: cache (containers.Map handle);
%   depP (struct, family parameters); orientDeg (scalar [deg]); phase
%   (scalar or array [fraction]).  OUTPUTS: x [6 x 1] (or [6 x N]).
key = okey(orientDeg);
if ~isKey(cache, key)
    p = depP;  p.orientDeg = orientDeg;
    [t, rv] = get_family_orbit('gto', p);
    cache(key) = phase_state(t, rv);
end
h = cache(key);
x = h(phase);
end

% ------------------------------------------------------------------------
function k = okey(orientDeg)
% OKEY  Cache key of an orientation.  INPUTS: orientDeg (scalar [deg]).
%   OUTPUTS: k (char).
k = sprintf('%.12g', orientDeg);
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
