%% GTO_DRO_STUDY  One minimum-time low-thrust transfer: GTO -> DRO
%   Edit the knobs, press Run. Every step is done here rather than behind a
%   front door so the machinery is visible; the SOLVE goes through the one
%   shared engine (lib/gto_dro_entry -> costate_common/thrust_ladder_library:
%   direct solve -> fly -> ms_tfmin -> pumpkyn tfMin acceptance) and is never
%   reimplemented. Five instruments stay in the library because a second copy
%   here would be a second, unverified implementation: the pointwise PMP
%   checks (N2, N4-N6), the hypothesis gates (S1-S3, X2), lift_margin (S3),
%   conj_spectrum (S4) and fly_transfer (the ONE flight every number below
%   reads).
%
%   The GTO "orbit" is get_family_orbit's algebraic Kepler locus frozen at
%   epoch zero in the rotating frame: a SET of departure states keyed by
%   orientation and time-since-perigee fraction, not a trajectory. The DRO
%   is pumpkynPie's getDRO seed refined by cont_np and propagated one period.
%
%     0  knobs                             (one source; study_opts overrides)
%     1  paths and constants
%     2  the DEPARTURE locus (GTO)         V1 geometry
%     3  the ARRIVAL orbit (DRO)           V2 closure
%     4  endpoints and the energy bound
%     5  solve down the ladder             N0 ladder, V3 ratio screen
%     6  the flight                        N3 arrival
%     7  clearances                        N7a Earth, N7b Moon
%     8  necessary conditions              N1, N2, N4, N5, N6
%     9  sufficiency hypotheses            S1-S4, X2
%    10  acceptance by pumpkyn             X1
%    11  figures
%    12  save and the verdict table
%   IDs: N necessary (fail one and this is not an extremal), S sufficiency
%   hypotheses (Bonnard-Caillau-Trelat; results, never thrown), V validity
%   of the inputs and instruments, X cross-checks against pumpkyn.
%
%   study_opts: a struct in the workspace before running overrides any knob
%   (e.g. study_opts.thrustN = 7; study_opts.out.visible = 'off'). Used by
%   tests/test_gto_dro_study and by anyone scripting several runs.
%
%  M. Casey                                                   (c) 11/09/2026
%  Copyright Coorbital Inc.

clearvars -except study_opts; clc
if ~exist('study_opts', 'var'), study_opts = struct(); end
here = fileparts(mfilename('fullpath'));

%% ========================================================================
%  0. KNOBS -- every adjustable parameter, with units; one source.
%% ========================================================================
thrustN = 15;                        % returned rung [N]; the ladder walks 15 N down to it
m0kg    = 150;                       % initial mass [kg]
ispS    = 1710;                      % specific impulse [s]
gto     = struct( ...
    'orientDeg', 0, ...              % Earth->Moon line to perigee [deg], rotating sense
    'sma_km',    24446, ...          % semi-major axis [km] (350 x 35786 km GTO)
    'ecc',       0.724781150290436, ...
    'depFrac',   0.10);              % departure phase, TIME fraction since perigee (0.5 = apogee). 0.10 by measurement: the perigee cell (0) at 15 N hits the iteration cap with ~5% mass left (diag 2026-11-09); 0.10 converges in ~20 s
dro     = struct( ...
    'tau',       1, ...              % DRO period [ND]; 0.5 1 2 3 are the catalog periods (1 = 4.43 d)
    'arrFrac',   0);                 % arrival phase: time fraction along the DRO from its +x crossing
mesh    = struct('N', 400);          % Hermite-Simpson intervals per rung
recipe  = struct( ...
    'thrLock',   true, ...           % pin the throttle to 1 (classical min-time)
    'tf0',       0.30, ...           % top-rung t_f guess [ND] (tuned for tulips; printed against the result)
    'maxCpuSec', 600, ...            % IPOPT wall cap per rung [s]
    'K',         24, ...             % multiple-shooting segments for the conjugate re-solve
    'floorEarthKm', 200);            % Earth altitude floor imposed in the direct solve [km] (engine option added 2026-11-09; the post-solve check N7a uses tol.earthKm)
tol     = struct( ...
    'gtoKm',     1, ...              % V1 perigee/apogee radius vs a(1 -+ e) [km]
    'closure',   1e-7, ...           % V2 DRO periodicity [ND]
    'ratioLo',   1/3, 'ratioHi', 3, ...   % V3 t_f ratio between consecutive rungs (steps only)
    'msNormR',   1e-8, ...           % N1 multiple-shooting residual, inf-norm
    'H',         1e-6, ...           % N2 |H| along the arc
    'km',        100, 'ms', 10, ...  % N3 flown arrival position [km] and velocity [m/s]
    'lamm',      1e-6, ...           % N4 |lambda_m(t_f)|
    'adj',       1e-7, ...           % N5 adjoint equations, relative
    'gap',       1e-12, 'throttle', 1e-10, ...   % N6 minimum-principle gap, applied throttle
    'earthKm',   200, ...            % N7a Earth altitude floor [km] (below the 350 km perigee by design)
    'moonKm',    500, ...            % N7b lunar altitude floor [km] (the engine's own floorKm)
    'lift',      1e-6, 'liftMargin', 10, ...     % S3 lift backward error, Eckart-Young margin
    'accDz',     1e-6, ...           % X1 pumpkyn tfMin must reproduce z8 to this
    'field',     1e-10);             % X2 pumpkyn field vs independent field, relative
out     = struct( ...
    'dir',       '', ...             % results folder ('' = results/study_o<orient>_tau<tau>_<thrust>N)
    'visible',   'on', ...           % figure visibility
    'conjTest',  true, ...           % run the conjugate re-solve (S4); false = NOT APPLICABLE
    'reuse',     '', ...             % path of a saved gto_dro_study.mat to gate instead of solving
    'solverType', 'fsolve');         % X1: pumpkyn tfMin solver ('fsolve' | 'lsqnonlin')
selfCheck = struct('strict', false); % true: any FAIL or UNRESOLVED gate throws (tests); false: code gates only

% study_opts overrides (struct knobs merge field by field)
for f = fieldnames(study_opts).'
    name = f{1};
    if exist(name, 'var') && isstruct(study_opts.(name)) && isstruct(eval(name))
        cur = eval(name);
        for g = fieldnames(study_opts.(name)).'
            cur.(g{1}) = study_opts.(name).(g{1});
        end
        eval([name ' = cur;']);
    else
        eval([name ' = study_opts.(name);']);
    end
end
if isempty(out.dir)
    out.dir = fullfile(here, 'results', sprintf('study_o%03d_tau%g_%gN', gto.orientDeg, dro.tau, thrustN));
end
if ~isfolder(out.dir), mkdir(out.dir); end
knobs = struct('thrustN', thrustN, 'm0kg', m0kg, 'ispS', ispS, 'gto', gto, 'dro', dro, ...
               'mesh', mesh, 'recipe', recipe, 'tol', tol, 'out', out);
verdicts = cell(0, 4);                                       % {id, status, value, threshold}
tStudy = tic;

%% ========================================================================
%  1. PATHS AND CONSTANTS
%% ========================================================================
fprintf('1. PATHS AND CONSTANTS\n');
run(fullfile(here, 'setup_gto_dro_paths.m'));
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
rEarthKm = 6378.137;  rMoonKm = 1737.4;
day    = @(tnd) tnd*tStar/86400;
ndp    = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
Tnd    = ndp.Tnd;   cnd = ndp.cnd;
fprintf('   mu* = %.12f, L* = %.3f km, T* = %.3f s (1 TU = %.4f d)\n', muStar, lStar, tStar, day(1));
fprintf('   thrust %.3g N on %g kg at Isp %g s -> T/m0 = %.3e ND, c = %.4f ND (%.3f km/s)\n', ...
        thrustN, m0kg, ispS, Tnd, cnd, cnd*lStar/tStar);

%% ========================================================================
%  2. THE DEPARTURE LOCUS (GTO)
%% ========================================================================
fprintf('\n2. DEPARTURE: GTO locus, orientation %g deg\n', gto.orientDeg);
[tD, rvD, infoD] = get_family_orbit('gto', struct('orientDeg', gto.orientDeg, 'sma_km', gto.sma_km, 'ecc', gto.ecc));
rDkm   = lStar*sqrt(sum((rvD(:, 1:3) - [-muStar, 0, 0]).^2, 2));
rPeri  = min(rDkm);  rApo = max(rDkm);
TkepD  = day(tD(end));
fprintf('   %d locus points; perigee altitude %.1f km, apogee altitude %.1f km, Kepler period %.4f d\n', ...
        size(rvD, 1), rPeri - rEarthKm, rApo - rEarthKm, TkepD);
fprintf('   (a SET of departure states frozen at epoch zero in the rotating frame, not a trajectory)\n');
v1val = max(abs(rPeri - gto.sma_km*(1 - gto.ecc)), abs(rApo - gto.sma_km*(1 + gto.ecc)));
verdicts = addv(verdicts, 'V1 GTO geometry', pass3(v1val < tol.gtoKm), v1val, tol.gtoKm);

%% ========================================================================
%  3. THE ARRIVAL ORBIT (DRO)
%% ========================================================================
fprintf('\n3. ARRIVAL: DRO, period tau = %g ND\n', dro.tau);
[tA, rvA, infoA] = get_family_orbit('dro', struct('tau', dro.tau));
closure = max(abs(rvA(end, :) - rvA(1, :)));
CJ      = pumpkyn.cr3bp.jacobi(rvA(1, :), muStar);
dAkm    = lStar*sqrt(sum((rvA(:, 1:3) - [1 - muStar, 0, 0]).^2, 2));
fprintf('   period %.4f d, Jacobi C = %.5f, perilune altitude %.0f km, apolune altitude %.0f km, closure %.1e ND\n', ...
        day(infoA.periodND), CJ, min(dAkm) - rMoonKm, max(dAkm) - rMoonKm, closure);
verdicts = addv(verdicts, 'V2 DRO closure', pass3(closure < tol.closure), closure, tol.closure);

%% ========================================================================
%  4. ENDPOINTS AND THE ENERGY BOUND
%% ========================================================================
fprintf('\n4. ENDPOINTS (costate_common/phase_state, spline rule)\n');
depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
rv0 = depAt(gto.depFrac).';     rvf = arrAt(dro.arrFrac).';
rv0km = pumpkyn.cr3bp.toPCI(0, rv0, muStar, tStar, lStar, 1);     % Earth-centred, km and km/s
rvfkm = pumpkyn.cr3bp.toPCI(0, rvf, muStar, tStar, lStar, 2);     % Moon-centred
sepKm = lStar*sqrt(sum((rvf(1:3) - rv0(1:3)).^2));
dVlb  = pumpkyn.cr3bp.minDeltaV(rv0, rvf, muStar)*lStar/tStar;
fprintf('   depart  (phase %.3f) ND [%s]\n', gto.depFrac, sprintf('%+.6f ', rv0));
fprintf('           Earth-centred: r = %.1f km, v = %.4f km/s\n', sqrt(sum(rv0km(1:3).^2)), sqrt(sum(rv0km(4:6).^2)));
fprintf('   arrive  (phase %.3f) ND [%s]\n', dro.arrFrac, sprintf('%+.6f ', rvf));
fprintf('           Moon-centred:  r = %.1f km, v = %.4f km/s\n', sqrt(sum(rvfkm(1:3).^2)), sqrt(sum(rvfkm(4:6).^2)));
fprintf('   straight-line separation %.0f km; pumpkyn minDeltaV energy bound %.3f km/s\n', sepKm, dVlb);

%% ========================================================================
%  5. SOLVE DOWN THE LADDER (lib/gto_dro_entry -> the shared engine)
%% ========================================================================
fprintf('\n5. SOLVE: 15 N down to %g N\n', thrustN);
entryOpts = struct('m0kg', m0kg, 'ispS', ispS, 'thrLock', recipe.thrLock, 'tf0', recipe.tf0, ...
    'maxCpuSec', recipe.maxCpuSec, 'N', mesh.N, 'K', recipe.K, 'conjTest', out.conjTest, ...
    'floorEarthKm', recipe.floorEarthKm, ...
    'gtoParams', struct('sma_km', gto.sma_km, 'ecc', gto.ecc), ...
    'outMat', fullfile(out.dir, 'ladder_scratch.mat'), 'keepScratch', true);
if isempty(out.reuse)
    tSolve = tic;
    [entry, gates, P] = gto_dro_entry(gto.orientDeg, dro.tau, gto.depFrac, dro.arrFrac, thrustN, entryOpts);
    solveWall = toc(tSolve);
    fprintf('   solved in %.0f s wall\n', solveWall);
else
    S = load(out.reuse, 'entry', 'gates', 'P', 'knobs');
    same = S.knobs.thrustN == thrustN && S.knobs.m0kg == m0kg && S.knobs.ispS == ispS && ...
           S.knobs.gto.orientDeg == gto.orientDeg && S.knobs.gto.depFrac == gto.depFrac && ...
           S.knobs.gto.sma_km == gto.sma_km && S.knobs.gto.ecc == gto.ecc && ...
           S.knobs.dro.tau == dro.tau && S.knobs.dro.arrFrac == dro.arrFrac;
    assert(same, 'gto_dro_study:reuse', ...
           '%s was solved with different knobs (thrust, mass, Isp, GTO or DRO parameters); refusing to gate it as this case', out.reuse);
    entry = S.entry;  gates = S.gates;  P = S.P;  solveWall = NaN;
    fprintf('   REUSED %s (no solve)\n', out.reuse);
end
fprintf('   %6s %9s %8s %8s %8s %9s %9s %9s %6s\n', 'N', 't_f ND', 't_f d', 'm_f kg', 'dV km/s', 'flown km', 'ms |R|', 'acc |dz|', 'ok');
for kr = 1:numel(P.rungs)
    tfk = P.TF(1, 1, kr);  mfk = m0kg*(1 - nd_propulsion(P.rungs(kr), ispS, m0kg, lStar, tStar).Tnd*tfk/cnd);
    dVk = cnd*log(m0kg/mfk)*lStar/tStar;
    fprintf('   %6.3g %9.4f %8.3f %8.2f %8.3f %9.3f %9.1e %9.1e %6d\n', P.rungs(kr), tfk, day(tfk), mfk, dVk, ...
            P.FLYKM(1, 1, kr), P.RES(1, 1, kr), P.ACCDZ(1, 1, kr), P.OK(1, 1, kr));
end
fprintf('   top-rung t_f / tf0 = %.2f (tf0 = %.2f ND was tuned for tulips)\n', entry.tfPerRung(1)/recipe.tf0, recipe.tf0);
verdicts = addv(verdicts, 'N0 ladder', pass3(all(entry.okPerRung)), nnz(entry.okPerRung), numel(entry.okPerRung));
if numel(P.rungs) > 1
    rs = gates.ratioScreen;
    fprintf('   step ratios t_f(k+1)/t_f(k): %s\n', sprintf('%.3f ', rs));
    verdicts = addv(verdicts, 'V3 ratio screen', pass3(all(rs > tol.ratioLo & rs < tol.ratioHi)), max(abs(log(rs))), log(tol.ratioHi));
else
    verdicts = addv(verdicts, 'V3 ratio screen', 'NOT APPLICABLE', NaN, NaN);
end
assert(entry.ok, 'gto_dro_study:noRoot', 'the bottom rung did not produce an accepted entry; nothing below can be gated');

%% ========================================================================
%  6. THE FLIGHT -- pumpkyn tfMinProp from the accepted z8 (ONE flight; every
%     number in sections 7-11 reads it). The engine's own G1b flown miss of
%     the DIRECT solution at this rung is printed beside it (P.FLYKM): the
%     direct control is not saved by the engine, so it is read, not re-flown.
%% ========================================================================
fprintf('\n6. FLIGHT of the accepted costates\n');
phys   = struct('Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'm0kg', m0kg);
flight = fly_transfer(entry.z8, rv0(1:6).', rvf(1:6).', phys);
fprintf('   t_f = %.4f ND = %.3f d; m_f = %.2f kg (%.2f kg propellant); dV = %.3f km/s (bound %.3f)\n', ...
        entry.z8(8), flight.tfDays, flight.finalMassKg, flight.propellantKg, flight.dvKms, dVlb);
fprintf('   arrival miss %.4f km, %.4f m/s (z8 flight); engine G1b flown miss of the direct solution %.4f km\n', ...
        flight.flyKm, flight.flyVms, gates.flownKm);
verdicts = addv(verdicts, 'N3 arrival', pass3(flight.flyKm < tol.km && flight.flyVms < tol.ms), flight.flyKm, tol.km);

%% ========================================================================
%  7. CLEARANCES on the flown arc (lib/earth_clearance, one instrument twice)
%% ========================================================================
fprintf('\n7. CLEARANCES\n');
[hEarth, tE, pE, oE] = earth_clearance(flight.t, flight.Y, tol.earthKm, struct('muStar', muStar, 'lStar', lStar));
[hMoon,  tM, pM, oM] = earth_clearance(flight.t, flight.Y, tol.moonKm, ...
                           struct('muStar', muStar, 'lStar', lStar, 'center', [1 - muStar, 0, 0], 'radiusKm', rMoonKm));
fprintf('   N7a Earth altitude min %9.1f km at t/t_f = %.3f  (floor %4.0f)  %s\n', hEarth, tE/flight.t(end), tol.earthKm, pass3(pE));
fprintf('   N7b Moon  altitude min %9.1f km at t/t_f = %.3f  (floor %4.0f)  %s\n', hMoon,  tM/flight.t(end), tol.moonKm,  pass3(pM));
verdicts = addv(verdicts, 'N7a Earth floor', pass3(pE), hEarth, tol.earthKm);
verdicts = addv(verdicts, 'N7b Moon floor',  pass3(pM), hMoon,  tol.moonKm);

%% ======================= sections 8-12 follow (Task 4) ===================
save(fullfile(out.dir, 'gto_dro_study.mat'), 'entry', 'gates', 'P', 'knobs', 'verdicts', 'flight');
fprintf('\n   (provisional) saved %s after %.0f s\n', fullfile(out.dir, 'gto_dro_study.mat'), toc(tStudy));

%% ------------------------------------------------------------------------
function v = fieldd(s, f, d_)
% FIELDD  s.(f) if present else d_.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = d_; end
end

function s = pass3(c)
% PASS3  PASS / FAIL from a logical.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function V = addv(V, id, status, value, threshold)
% ADDV  Append one verdict row {id, status, value, threshold}.  INPUTS: V; id; status; value; threshold.  OUTPUTS: V.
V(end+1, :) = {id, status, value, threshold};
end
