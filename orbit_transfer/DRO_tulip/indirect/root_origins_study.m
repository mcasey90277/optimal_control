%% ROOT_ORIGINS_STUDY  Where a root comes from: COLD LOTTERY -> LADDER -> BASIN HUNT.
%
%   THE PREQUEL TO anchor_study.m. That script starts from a certified root
%   of a neighbouring transfer and makes a new one. This script starts from
%   nothing, and answers the question underneath the whole library: if every
%   seed is a previous root, where did the FIRST root come from?
%
%   HONEST PROVENANCE FIRST. The very first DRO -> tulip root on this
%   machine was not built here. It arrived with pumpkynPie as a converged
%   indirect solution of exactly this cell -- Darin's walk-down -- and it is
%   the "reference" every early FINDINGS table is measured against. This
%   script is the route you take when no such gift exists, which is the
%   situation for every other cell, every other thrust and every other
%   orbit pair in the program.
%
%   THREE MECHANISMS, measured here rather than asserted:
%
%     A. THE COLD LOTTERY (section 2). Solve the operating point cold, from
%        nothing but the endpoints and a time-of-flight guess, at several
%        mesh densities. Every mesh converges. They do NOT agree. Mesh
%        density acts as a de facto random seed for which basin you land in,
%        and no output of the solver tells you which one you got. This is
%        the measurement that makes continuation required equipment rather
%        than a convenience.
%
%     B. THE LADDER (sections 3-5). Start where the problem is EASY -- high
%        thrust, near-impulsive, sub-revolution, converges cold in seconds
%        -- and walk the ENGINE down to where you actually want it, each
%        rung warm-started from the rung above. The walk has two halves and
%        the seam between them is the whole subject of anchor_study.m:
%          3. a DIRECT ladder, 15 N -> 0.5 N, each rung a fresh collocation
%             NLP seeded by the previous rung's discrete trajectory;
%          4. the HANDOFF: harvest the bottom rung's constraint multipliers
%             into costates, shoot them, and you have an indirect root;
%          5. an INDIRECT ladder, 0.5 N -> the operating point, each rung a
%             multiple-shooting solve seeded by the previous rung's junction
%             states. Below about 0.5 N the direct NLP loses the thread and
%             the shooting solver does not; that is why the seam is there.
%
%     C. THE BASIN HUNT (section 7). At the operating point you have
%        reached, ask the direct solver for a transfer to a DIFFERENT
%        arrival phase, warm-started from the root you hold. The solver is
%        BRANCH-BLIND -- it lands in whatever basin its guess is near. That
%        is a hazard if you assume it returns the family you started from,
%        and a discovery tool if you do not. It is how four of the five
%        families in the 70 mN library were found (FINDINGS 61, 63, 68).
%
%   THE SEQUENCE OF STUDY SCRIPTS:
%     root_origins_study   (this)  where a root comes from when there is none
%     anchor_study                 how a root becomes a SEED, and a new root
%     transfer_study               is that root a minimum? the gate stack
%     run_phase_torus              a whole library over a phase grid
%
%   RUNTIME AND WHAT THE DEFAULTS ACTUALLY DO. About 50 minutes, measured
%   2026-09-15 on this cell under R2026a. Every rung prints as it lands and
%   the run saves after each one, so it can be watched and picked apart.
%   The defaults reach 0.12 N and then WALL: rungs at 0.11 and 0.10 N refuse,
%   with revolutions climbing 0.81 -> 1.31 across the walk. That is a winding
%   wall and it is a property of THIS CELL -- the campaign's own chain reached
%   0.09 N and stalled at 0.067 N, walking from the FASTEST 0.5 N entry of the
%   sheet rather than from the anchor phase. The run still certifies the root
%   it reaches, and section 6 says plainly what was and was not demonstrated.
%   Fewer meshes in section 0 and a shorter ladder cut this to a few minutes.
%
%  M. Casey                                                   (c) 09/15/2026
%  Copyright Coorbital Inc.

%% paths
clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
D = fullfile(fileparts(here), 'direct');
addpath(D, fullfile(D, 'lib'), fullfile(D, 'certify'));
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));

%% CR3BP constants
muStar = 0.012150585609624;                 % Earth-Moon mass ratio
lStar  = 389703.264829278;                  % km
tStar  = 382981.289129055;                  % s
day    = @(tnd) tnd*tStar/86400;

%% ========================================================================
%  0. USER INPUTS
%     Everything the run depends on is here: the two orbits, the cell, the
%     engine we want to arrive at, the two ladders that get us there, where
%     to hunt, and the budgets and tolerances every printed gate uses.
%% ========================================================================

%% the problem
dep   = struct('family', 'dro',   'tau', 1.0);            % departure orbit
arr   = struct('family', 'tulip', 'Np', 7, 'pm', -1);     % arrival orbit
phase = struct('sD', 0.0, 'sA', 0.0754);                  % the cell we solve

%% the OPERATING POINT we want to reach, and the campaign's answer for it
op = struct( ...
    'thrustN', 0.070, ...    % N   -- the abstract's Hall thruster
    'ispS',    900,   ...    % s
    'm0kg',    150,   ...    % kg
    'tfDays',  17.7976);     % d   -- the campaign's certified root at this cell,
                             %        reported against, NEVER used to steer

%% A. the cold lottery: the same cold solve at several mesh densities
cold = struct( ...
    'run',     true, ...
    'meshes',  [400 800 1600], ...   % one cold solve per mesh; 1600 is the slow one
    'tf0',     4.0,  ...     % ND time-of-flight guess (the solver's own default)
    'floor',   false);       % the recorded lottery ran UNCONSTRAINED, so a
                             % Moon-diving basin can show itself. true adds the
                             % clearance path constraint to the cold solves.

%% B. the DIRECT ladder: each row is an engine setting, solved by collocation
%     warm-started from the row above. Isp stays high up here on purpose: at
%     15 N with a 900 s engine the vehicle would burn most of its mass.
%     THE RUNG RATIOS MATTER. These are the shipped catalog's, about 0.75
%     per step. A coarser set jumps basins: measured on this cell
%     2026-09-15, a 1 -> 0.5 N step (ratio 0.5) returned a converged,
%     defect-1e-14, clearance-safe extremal at 54.7 d and dV 46.8 km/s --
%     16x its own guess, 94% of the mass burned, a different winding number
%     entirely. Section 3's screen is what catches that.
dladder = struct( ...
    'rungs', [15 1710; 12 1710; 10 1710; 7 1710; 5 1710; 3 1710; ...
              2 1710; 1.5 1710; 1 1710; 0.75 1710; 0.5 1710], ...        % [N, s]
    'tfExp', 0.6, ...        % next rung's guess: t_f ~ T^-0.6 in this band
    'tfBand', 3, ...         % refuse a rung outside [1/band, band] x its guess
    'N',     800);           % collocation intervals

%% B2. the INDIRECT ladder: below ~0.5 N the direct NLP loses the thread, so
%      the walk continues by multiple shooting on banked junction states.
%      This is the real route: thrust down at high Isp, then Isp down at
%      fixed thrust, then the last thrust steps. The chain stalled at
%      0.067 N on a basin/winding wall, so the bottom rungs are where this
%      stops being a demonstration and becomes a search.
%
%      WHERE THE Isp STAGE GOES IS A DESIGN CHOICE, and it was measured
%      here on 2026-09-15. Walking thrust all the way down at Isp 1710 s
%      walls at 0.105 N on this cell (revolutions climbing 0.81 -> 1.31: a
%      winding wall, not a sensitivity one). Taking the Isp stage EARLIER,
%      while the walk is still converging, is the remedy probe_abstract_case
%      used: dropping Isp raises the depletion rate, the vehicle gets
%      lighter sooner, and late-arc acceleration goes UP -- which helps the
%      descent that follows rather than hurting it.
iladder = struct( ...
    'run',   true, ...
    'rungs', [0.375 1710; 0.30 1710; 0.24 1710; 0.20 1710; 0.16 1710; ...
              0.14  1710; 0.12 1710; ...
              0.12  1400; 0.12 1150; 0.12  900; ...
              0.11   900; 0.10  900; 0.09  900; 0.08  900; 0.07  900], ... % [N, s]
    'tfExps', [1.0 1.3 1.7 2.2 3.0], ...   % per-rung t_f guess sweep: past a winding
    'attemptSec', 600, ...                 % per-rung budget across all exponents
    'maxMiss', 2);                         % wall the right answer has a different
                                           % number of revolutions, and only a
                                           % different exponent reaches its basin.
                                           % maxMiss: consecutive refusals = a wall

%% C. the basin hunt at whatever operating point we reached
%     A hunt probe is a lottery ticket, not a required result, so it gets a
%     short budget: a probe that has not converged in a few minutes is
%     telling you there is no nearby basin, and the next phase is cheaper to
%     try than this one is to finish.
hunt = struct( ...
    'run',    true, ...
    'cpuSec', 300, ...     % per-probe IPOPT budget
    'sA',     [0.0754 + 6/24, 0.0754 + 12/24, 0.0754 + 17/24]);   % arrival phases to probe

%% budgets
num = struct('maxIter', 3000, 'maxCpuSec', 900, 'coldCpuSec', 1800, ...
             'clearKm', 1900, 'K', 24, 'shootSec', 240, 'certSec', 900, 'gateKm', 100);

%% tolerances -- one source, used by every printed line AND the verdict
tol = struct( ...
    'closure', 1e-7,  ...    % orbit periodicity
    'defect',  1e-8,  ...    % collocation defect on an accepted rung
    'termErr', 1e-7,  ...    % terminal boundary residual
    'unit',    1e-8,  ...    % |u| = 1
    'R',       1e-8,  ...    % multiple-shooting residual, inf-norm
    'basin',   0.02,  ...    % |t_f - t_f_ref|/t_f above which we call it ANOTHER basin
    'same',    1e-2);        % |t_f - the campaign's|/t_f to call it the same root

outMat = fullfile(here, 'results', 'root_origins_study.mat');
g0  = 9.80665*tStar^2/(1000*lStar);
ndT = @(TN, m) (TN/m)*tStar^2/(lStar*1000);       % N -> ND thrust acceleration at m = 1
ndC = @(isp)   (isp/tStar)*g0;                    % s -> ND exhaust speed
rMoonKm = 1737.4;  floorKm = num.clearKm - rMoonKm;

fprintf('ROOT ORIGINS -- how a transfer is found when no solution exists yet\n');
fprintf('   problem  : DRO tau %.2f -> %d-petal tulip (pm %+d) at sD %.4f, sA %.4f\n', ...
        dep.tau, arr.Np, arr.pm, phase.sD, phase.sA);
fprintf('   target   : %.0f mN, Isp %g s, %g kg   (the campaign certifies %.4f d here)\n', ...
        op.thrustN*1000, op.ispS, op.m0kg, op.tfDays);
fprintf('   A lottery: %s\n', tern(cold.run, sprintf('meshes %s, t_f guess %.2f ND, clearance floor %s', ...
        mat2str(cold.meshes), cold.tf0, tern(cold.floor, 'ON', 'OFF')), 'off'));
fprintf('   B direct : %d rungs, %g N -> %g N at Isp %g s, N = %d, guess band %gx\n', ...
        size(dladder.rungs, 1), dladder.rungs(1, 1), dladder.rungs(end, 1), dladder.rungs(1, 2), ...
        dladder.N, dladder.tfBand);
fprintf('   B2 indir : %s\n', tern(iladder.run, sprintf('%d rungs, %g N -> %g N, Isp %g -> %g s, %g s per rung', ...
        size(iladder.rungs, 1), iladder.rungs(1, 1), iladder.rungs(end, 1), ...
        iladder.rungs(1, 2), iladder.rungs(end, 2), iladder.attemptSec), 'off'));
fprintf('   C hunt   : %s\n', tern(hunt.run, sprintf('arrival phases %s, %g s per probe', ...
        mat2str(hunt.sA, 4), hunt.cpuSec), 'off'));

%% ========================================================================
%  1. THE TWO ORBITS AND THE TWO ENDPOINTS
%     Both orbits are generated from their own parameters, never loaded, so
%     the run is reproducible from section 0 alone.
%% ========================================================================
[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
closD = norm(rvD(end, :) - rvD(1, :));   closA = norm(rvT(end, :) - rvT(1, :));
[stateD, seamD] = phase_state(tD, rvD);
[stateA, seamA] = phase_state(tT, rvT);
rv0 = stateD(phase.sD);   rv0 = rv0(1:6);
rvf = stateA(phase.sA);
orbOK = closD < tol.closure && closA < tol.closure && max(seamD.deriv, seamA.deriv) < 1e-6;
fprintf('\n1. ORBITS: DRO period %.3f d, tulip period %.3f d\n', day(infoD.periodND), day(infoT.periodND));
fprintf('   closures %.1e / %.1e, interpolant seam derivative %.1e            %s\n', ...
        closD, closA, max(seamD.deriv, seamA.deriv), pass(orbOK));
fprintf('   departure r = [%+.5f %+.5f %+.5f], arrival r = [%+.5f %+.5f %+.5f]\n', rv0(1:3), rvf(1:3));
assert(orbOK, 'the generated orbits are not periodic to tolerance');

TndOp = ndT(op.thrustN, op.m0kg);   cndOp = ndC(op.ispS);
S = struct('phase', phase, 'op', op, 'lottery', [], 'direct', [], 'handoff', [], ...
           'indirect', [], 'reached', [], 'hunt', []);

%% ========================================================================
%  2. THE COLD LOTTERY -- the measurement that makes continuation necessary.
%
%     X0 = U0 = [] hands the NLP nothing but the endpoints, a time-of-flight
%     guess and the mesh. It builds its own crude guess (states linear
%     between the endpoints, a mass ramp, thrust along the chord) and knows
%     nothing about costates or about what the answer looks like.
%
%     Nothing about the PROBLEM changes between these runs. Only the
%     discretisation does. If they disagree on t_f then a cold solve is
%     choosing a basin for reasons that have nothing to do with the physics,
%     and "it converged and the defects are tiny" is not evidence that you
%     have the transfer you wanted.
%
%     The recorded result on this exact cell (2026-08-03, FINDINGS):
%        N = 400  -> 4.3807 ND, +9.1%,  periselene 4819 km
%        N = 800  -> 4.4507 ND, +10.8%, periselene -343 km, INSIDE THE MOON
%        N = 1600 -> 4.7833 ND, +19.1%, 7.8 m worst-interval error, PHYSICAL
%     The N = 1600 row is the decisive one: a genuinely accurate, physical,
%     safe minimum-time extremal that is 19% slower than the reference. The
%     accuracy machinery works cold. BASIN SELECTION DOES NOT.
%% ========================================================================
dCold = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', false, ...
               'returnModel', true, 'maxCpuSec', num.coldCpuSec);
if cold.floor, dCold.minAltKm = floorKm; end
L = struct('N', {}, 'tf', {}, 'defect', {}, 'perisKm', {}, 'rel', {}, 'ok', {}, 'wall', {});
if cold.run
    fprintf('\n2. THE COLD LOTTERY at %.0f mN / Isp %g s -- one solve per mesh, same problem\n', ...
            op.thrustN*1000, op.ispS);
    fprintf('      %6s  %10s  %9s  %8s  %10s  %11s  %s\n', ...
            'N', 't_f (ND)', 't_f (d)', 'vs ref', 'defect', 'peris (km)', 'status');
    for N = cold.meshes
        tM = tic;
        try
            oM = solveDirect(rv0, rvf, TndOp, cndOp, muStar, N, [], [], cold.tf0, dCold);
            okM = oM.success && oM.maxDefect < tol.defect;
            relM = (oM.tf*tStar/86400 - op.tfDays)/op.tfDays;
            L(end+1) = struct('N', N, 'tf', oM.tf, 'defect', oM.maxDefect, ...
                              'perisKm', oM.altMinKm + rMoonKm, 'rel', relM, 'ok', okM, 'wall', toc(tM));
            fprintf('      %6d  %10.4f  %9.3f  %+7.1f%%  %10.1e  %11.0f  %s (%.0f s)\n', ...
                    N, oM.tf, day(oM.tf), 100*relM, oM.maxDefect, oM.altMinKm + rMoonKm, ...
                    tern(okM, 'converged', oM.ipoptStatus), toc(tM));
        catch ME
            L(end+1) = struct('N', N, 'tf', NaN, 'defect', NaN, 'perisKm', NaN, ...
                              'rel', NaN, 'ok', false, 'wall', toc(tM));
            fprintf('      %6d  THREW: %s\n', N, firstline(ME.message));
        end
        S.lottery = L;  saveq(outMat, S);
    end
    good = L([L.ok]);
    if numel(good) >= 2
        spread = (max([good.tf]) - min([good.tf]))/min([good.tf]);
        nBasin = numel(uniquetol([good.tf], tol.basin));
        below  = nnz([good.perisKm] < num.clearKm);
        fprintf('   %d of %d meshes converged. t_f spread %.1f%%; %d distinct basin(s) at a %.0f%% threshold;\n', ...
                numel(good), numel(L), 100*spread, nBasin, 100*tol.basin);
        fprintf('   %d of them came inside the %g km lunar clearance.\n', below, num.clearKm);
        fprintf('   every one of these is a converged minimum-time extremal. None of them is\n');
        fprintf('   the campaign''s root, and nothing the solver printed says which is which.\n');
        if nBasin == 1
            fprintf('   (These meshes happened to agree. That is a result about THESE meshes,\n');
            fprintf('    not a property of cold solving -- see the recorded table above.)\n');
        end
    else
        fprintf('   fewer than two meshes converged: the lottery cannot be read from this run\n');
    end
else
    fprintf('\n2. THE COLD LOTTERY: skipped (cold.run = false)\n');
end

%% ========================================================================
%  3. THE DIRECT LADDER -- from where the problem is easy toward where we
%     want it. At 15 N the transfer is near-impulsive and sub-revolution and
%     the NLP converges cold in seconds. Each rung below it is a fresh
%     collocation solve warm-started from the rung above, with a
%     time-of-flight guess scaled by T^-exp. Nothing here is a search: every
%     rung starts next to its answer. Compare these wall times with the
%     lottery's, and notice that the answers now form a curve instead of a
%     scatter.
%% ========================================================================
dLad = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', false, ...
              'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', num.maxCpuSec);
fprintf(['\n3. THE DIRECT LADDER (N = %d, t_f guess scaled by T^-%.1f, clearance floor ON,\n' ...
         '   a rung refused if its t_f leaves [1/%g, %g] x its own guess)\n'], ...
        dladder.N, dladder.tfExp, dladder.tfBand, dladder.tfBand);
fprintf('      %8s %7s  %10s  %9s  %10s  %11s  %8s  %s\n', ...
        'T (N)', 'Isp (s)', 't_f (ND)', 't_f (d)', 'defect', 'peris (km)', 'dV km/s', 'status');
R = struct('thrustN', {}, 'ispS', {}, 'tf', {}, 'defect', {}, 'perisKm', {}, 'dv', {}, 'wall', {});
seedX = [];  seedU = [];  seedTf = cold.tf0;  Tprev = [];  oPrev = [];
for kr = 1:size(dladder.rungs, 1)
    TN = dladder.rungs(kr, 1);  isp = dladder.rungs(kr, 2);
    Tnd = ndT(TN, op.m0kg);  cnd = ndC(isp);
    if ~isempty(Tprev), seedTf = seedTf*(Tprev/TN)^dladder.tfExp; end
    tR = tic;
    try
        oR = solveDirect(rv0, rvf, Tnd, cnd, muStar, dladder.N, seedX, seedU, seedTf, dLad);
        perisR = oR.altMinKm + rMoonKm;
        % THE SCREEN. A ladder's premise is that each rung lands NEAR its
        % guess; a rung that does not has left the branch, and seeding the
        % next rung from it walks a different family down. It can look
        % perfect while doing so -- converged, defect 1e-14, clear of the
        % Moon -- so the only thing that catches it is comparing the answer
        % with the guess that produced it.
        % The band applies to a STEP, never to the top rung: that one's guess
        % is a cold round number with no branch behind it, and at 15 N the
        % true answer is ~40x below it. preflight_screen takes seedTf = []
        % for the same reason.
        band = oR.tf/seedTf;
        inBand = isempty(Tprev) || (band <= dladder.tfBand && band >= 1/dladder.tfBand);
        okR = oR.success && oR.maxDefect < tol.defect && oR.termErr < tol.termErr && ...
              oR.maxUnit < tol.unit && perisR >= num.clearKm && inBand;
        dv = cnd*log(1/oR.mf)*lStar/tStar;               % from the mass actually used
        if okR,           note = 'accepted';
        elseif ~inBand,   note = sprintf('REFUSED: %.1fx its own guess -- BASIN JUMP', band);
        elseif ~oR.success, note = ['REFUSED: ' oR.ipoptStatus];
        else,             note = 'REFUSED: defect, terminal, |u| or clearance';
        end
        fprintf('      %8.3f %7.0f  %10.4f  %9.3f  %10.1e  %11.0f  %8.4f  %s (%.0f s)\n', ...
                TN, isp, oR.tf, day(oR.tf), oR.maxDefect, perisR, dv, note, toc(tR));
        if ~isempty(Tprev)
            fprintf('               (guess %.4f ND from the rung above, answer %.2fx it)\n', seedTf, band);
        end
        if ~okR
            fprintf('      the direct ladder STOPS here: this rung cannot seed the next.\n');
            if ~inBand
                fprintf('      A BASIN JUMP is a statement about the STEP, not the problem: the\n');
                fprintf('      remedy is a finer rung ratio, not a looser gate. The shipped\n');
                fprintf('      catalog steps by about 0.75 per rung for exactly this reason.\n');
            end
            break
        end
        R(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', oR.tf, 'defect', oR.maxDefect, ...
                          'perisKm', perisR, 'dv', dv, 'wall', toc(tR));
        seedX = oR.X;  seedU = oR.U;  seedTf = oR.tf;  Tprev = TN;  oPrev = oR;
        S.direct = R;  saveq(outMat, S);
    catch ME
        fprintf('      %8.3f %7.0f  THREW: %s\n', TN, isp, firstline(ME.message));
        break
    end
end
b1 = ~isempty(R);
assert(b1, 'the direct ladder produced no usable rung: nothing below is meaningful');
b2 = R(end).thrustN == dladder.rungs(end, 1);
fprintf('   B1 direct ladder reached %g N of %g N requested                    %s\n', ...
        R(end).thrustN, dladder.rungs(end, 1), pass(b2));
if numel(R) >= 2
    fprintf('   t_f grew %.3f -> %.3f d as thrust fell %g -> %g N; dV %.4f -> %.4f km/s.\n', ...
            day(R(1).tf), day(R(end).tf), R(1).thrustN, R(end).thrustN, R(1).dv, R(end).dv);
    fprintf('   It is the same burn made longer, so dV rises slowly while t_f rises fast.\n');
end

%% ========================================================================
%  4. THE HANDOFF -- direct to indirect, which is where the costates enter.
%
%     The multiplier attached to a dynamics-defect constraint measures how
%     much the objective would improve if that constraint were relaxed,
%     which is the definition of the adjoint variable. So the NLP has been
%     computing costates all along under another name. Harvest them onto a
%     shooting mesh and you have a SEED: a state-and-costate trajectory on
%     K+1 junctions, not a vector. Shoot it and you have an indirect root.
%     anchor_study.m opens every step of this box; here we need the root.
%% ========================================================================
TndB = ndT(R(end).thrustN, op.m0kg);   cndB = ndC(R(end).ispS);
% THE BOTTOM RUNG, RE-SOLVED ONCE IN THE SUNDMAN CHART. The ladder runs in
% plain time, as the shipped catalog did. Under a Sundman change of variable
% the defect system carries the time state as well, so its multipliers carry
% ONE MORE ROW: lambda_t, which the minimum principle fixes at +1 when the
% objective is t_f. That row is a free check on the whole duals-to-costates
% mapping -- station association, sign and scale together -- and it simply
% does not exist without Sundman. One extra warm solve buys it.
try
    dSun = dLad;  dSun.sundman = true;
    oH8 = solveDirect(rv0, rvf, TndB, cndB, muStar, dladder.N, oPrev.X, oPrev.U, oPrev.tf, dSun);
    if ~(oH8.success && oH8.maxDefect < tol.defect), oH8 = []; end
catch
    oH8 = [];
end
oHarv = oPrev;  chart = 'plain time (no lambda_t row)';
if ~isempty(oH8), oHarv = oH8;  chart = 'Sundman (lambda_t exposed)'; end
[seedB, dgB] = harvest_ms_seed(oHarv, num.K);
tH = tic;
[zB, itB] = ms_tfmin(rv0, rvf(1:6), seedB, TndB, cndB, muStar, ...
                     struct('tolR', 1e-11, 'wallSec', num.shootSec, 'maxIter', 100));
h1 = isfinite(itB.normR) && itB.normR < tol.R;
fprintf('\n4. THE HANDOFF at %g N / Isp %g s: %d multipliers -> a %dx%d seed -> a root\n', ...
        R(end).thrustN, R(end).ispS, size(oHarv.lamDef, 2), size(seedB.Y, 1), size(seedB.Y, 2));
fprintf('   harvested from the %s solve\n', chart);
fprintf('   sign vote %.1f%% of stations; lambda_t %s\n', 100*gvd(dgB, 'voteMargin', NaN), ...
        tern(isfinite(gvd(dgB, 'lamT', NaN)), sprintf('%.6f (theory: +1)', gvd(dgB, 'lamT', NaN)), ...
             'NOT EXPOSED without Sundman -- the scale check is unavailable here'));
fprintf('   H1 shooting |R| %9.2e / %-9.0e in %d iterations, %.0f s    %s\n', ...
        itB.normR, tol.R, itB.iters, toc(tH), pass(h1));
fprintf('   t_f: direct %.4f d -> shot %.4f d\n', day(oHarv.tf), day(zB(8)));
fprintf('   lambda(0) = [%s]\n', strjoin(compose('%+.6g', zB(1:7)'), ' '));
assert(h1, 'the handoff shooting residual %.2e is above the gate %.0e: there is no root to walk', ...
       itB.normR, tol.R);
S.handoff = struct('thrustN', R(end).thrustN, 'ispS', R(end).ispS, 'chart', chart, ...
                   'z', zB(:), 'Y', itB.Y, 'tGrid', itB.tGrid, 'normR', itB.normR, ...
                   'voteMargin', gvd(dgB, 'voteMargin', NaN), 'lamT', gvd(dgB, 'lamT', NaN));
saveq(outMat, S);

%% ========================================================================
%  5. THE INDIRECT LADDER -- the rest of the way down.
%
%     Below about 0.5 N the collocation NLP loses the thread: the arc winds
%     more, the mesh has to resolve more revolutions, and the basin gets
%     narrow. Multiple shooting does not care nearly as much, because each
%     of its K segments is short. So the walk continues in the indirect
%     chart: each rung re-cuts the PREVIOUS rung's flown PMP trajectory onto
%     a new time-of-flight, rebuilds the mass row from the all-burn law for
%     the NEW engine, and solves ms_tfmin.
%
%     The per-rung t_f GUESS is swept over several exponents. That sweep is
%     not decoration: at a winding wall the right next solution has a
%     different number of revolutions, and only a different exponent puts
%     the guess in its basin.
%
%     A refused rung is a finding about the problem, not a bug. The recorded
%     chain stalled at 0.067 N. Where THIS run stops is printed and kept.
%% ========================================================================
zCur = zB;  Tcur = R(end).thrustN;  cCur = cndB;  ispCur = R(end).ispS;  itCur = itB;
W = struct('thrustN', {}, 'ispS', {}, 'tf', {}, 'normR', {}, 'flyKm', {}, 'revs', {}, ...
           'tfExp', {}, 'wall', {}, 'z', {}, 'Y', {}, 'tGrid', {});
pool = capped_pool();      % [] without the Parallel Computing Toolbox: then
                           % every solver call below runs UNFENCED and a single
                           % crawling integration can stall the run
if iladder.run
    fprintf('\n5. THE INDIRECT LADDER (multiple shooting, K = %d segments)\n', num.K);
    fprintf('      %8s %7s  %10s  %9s  %10s  %9s  %6s  %s\n', ...
            'T (N)', 'Isp (s)', 't_f (ND)', 't_f (d)', '|R|', 'flown km', 'revs', 'status');
    nMiss = 0;
    for kr = 1:size(iladder.rungs, 1)
        TN = iladder.rungs(kr, 1);  isp = iladder.rungs(kr, 2);
        Tnd = ndT(TN, op.m0kg);  cnd = ndC(isp);
        tR = tic;  best = [];
        [tj, yj] = pumpkyn.cr3bp.tfMinProp(zCur(8), [rv0(:); 1; zCur(1:7)], ...
                                           ndT(Tcur, op.m0kg), cCur, muStar);
        for tfExp = iladder.tfExps
            if toc(tR) > iladder.attemptSec, break, end  % this rung has had its budget
            if TN == Tcur, tfGuess = zCur(8);            % an Isp-only step: keep t_f
            else,          tfGuess = zCur(8)*(Tcur/TN)^tfExp;
            end
            if tfGuess > 0.6*cnd/Tnd                      % the vehicle would run out of mass
                continue
            end
            [Yg, tGs] = flight_to_junctions(tj, yj, num.K, ...
                            struct('tf', tfGuess, 'massLaw', struct('Tnd', Tnd, 'cnd', cnd)));
            [okRun, zt, it] = fenced(pool, @ms_tfmin, 2, num.shootSec + 90, ...
                rv0, rvf(1:6), struct('tf', tfGuess, 'tGrid', tGs, 'Y', Yg), ...
                Tnd, cnd, muStar, struct('wallSec', num.shootSec));
            if ~okRun || ~it.converged, continue, end
            [~, yFly] = pumpkyn.cr3bp.tfMinProp(zt(8), [rv0(:); 1; zt(1:7)], Tnd, cnd, muStar);
            mk = norm(yFly(end, 1:3).' - rvf(1:3))*lStar;
            if mk >= num.gateKm, continue, end
            if isempty(best) || zt(8) < best.z(8)
                best = struct('z', zt, 'it', it, 'mk', mk, 'tfExp', tfExp);
            end
            if TN == Tcur, break, end                     % an Isp step needs no sweep
        end
        if isempty(best)
            nMiss = nMiss + 1;
            fprintf('      %8.3f %7.0f  %10s  %9s  %10s  %9s  %6s  REFUSED (%.0f s) [miss %d]\n', ...
                    TN, isp, '-', '-', '-', '-', '-', toc(tR), nMiss);
            if nMiss >= iladder.maxMiss
                fprintf('      THE WALL: %d consecutive refused rungs at and above %g N.\n', nMiss, TN);
                fprintf('      A wall is a result. The recorded chain met one at 0.067 N, walking\n');
                fprintf('      from the FASTEST 0.5 N cell of the sheet -- deep descent is easier\n');
                fprintf('      on a fast cell, and WHERE you start the deep walk is a choice, not\n');
                fprintf('      a detail. The other ways past a wall are a finer rung ratio, a\n');
                fprintf('      different t_f exponent (already swept here), or a different ROUTE:\n');
                fprintf('      the Isp stage exists for exactly that reason.\n');
                break
            end
            continue
        end
        [~, yb] = pumpkyn.cr3bp.tfMinProp(best.z(8), [rv0(:); 1; best.z(1:7)], Tnd, cnd, muStar);
        ang = unwrap(atan2(yb(:, 2), yb(:, 1) - (1 - muStar)));
        revs = sum(abs(diff(ang)))/(2*pi);
        fprintf('      %8.3f %7.0f  %10.4f  %9.3f  %10.1e  %9.1f  %6.2f  accepted, tfExp %.1f (%.0f s)\n', ...
                TN, isp, best.z(8), day(best.z(8)), best.it.normR, best.mk, revs, best.tfExp, toc(tR));
        % THE ROOT ITSELF IS KEPT, not just its t_f: bare z8 pins t_f only to
        % ~1e-4 ND at many revolutions, so a rung is banked with its full
        % K+1 junction states (the identifiability rule).
        W(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', best.z(8), 'normR', best.it.normR, ...
                          'flyKm', best.mk, 'revs', revs, 'tfExp', best.tfExp, 'wall', toc(tR), ...
                          'z', best.z(:), 'Y', best.it.Y, 'tGrid', best.it.tGrid);
        zCur = best.z;  itCur = best.it;  Tcur = TN;  cCur = cnd;  ispCur = isp;  nMiss = 0;
        S.indirect = W;  saveq(outMat, S);
    end
else
    fprintf('\n5. THE INDIRECT LADDER: skipped (iladder.run = false)\n');
end
reached = struct('thrustN', Tcur, 'ispS', ispCur, 'tfDays', day(zCur(8)), 'z', zCur(:), ...
                 'Y', itCur.Y, 'tGrid', itCur.tGrid);
S.reached = reached;  saveq(outMat, S);
b3 = abs(Tcur - op.thrustN)/op.thrustN < 1e-6 && abs(ispCur - op.ispS) < 1e-6;
fprintf('   B2 the ladder reached %g N / Isp %g s; the target was %g N / %g s    %s\n', ...
        Tcur, ispCur, op.thrustN, op.ispS, pass(b3));

%% ========================================================================
%  6. THE ROOT WE HOLD -- certify it, and compare with the campaign.
%     Branch-blindness means "is it the campaign's root?" is a MEASUREMENT.
%     A certified root at a different t_f is another family, not an error.
%% ========================================================================
TndR = ndT(Tcur, op.m0kg);
seedR = struct('tf', zCur(8), 'tGrid', itCur.tGrid, 'Y', itCur.Y);
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', Tcur, 'ispS', ispCur, 'm0kg', op.m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', TndR, 'cnd', cCur, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
C = certify_root(seedR, rv0, rvf, B, ...
                 struct('pool', pool, 'wallSec', num.certSec, 'sA', phase.sA, 'sD', phase.sD));
fprintf('\n6. CERTIFY the %g N / Isp %g s root: %s\n', Tcur, ispCur, C.reason);
fprintf('   t_f %.4f d, flown miss %.3f km / %.3f m/s, conjugate %s, lift margin %.1fx\n', ...
        C.tfDays, C.flyKm, C.flyVms, tern(gvd(C, 'conj', 0) == 1, 'PASS', 'not passed'), ...
        gvd(C, 'liftMargin', NaN));
if b3
    sameRoot = C.ok && abs(C.tfDays - op.tfDays)/op.tfDays < tol.same;
    fprintf('   this run %.4f d vs the campaign''s %.4f d: %.2f%% apart -- %s\n', ...
            C.tfDays, op.tfDays, 100*abs(C.tfDays - op.tfDays)/op.tfDays, ...
            tern(sameRoot, 'THE SAME ROOT, found without ever being given one', ...
                 'A DIFFERENT ROOT (another family, or another basin)'));
    if cold.run && ~isempty(L) && any([L.ok])
        goodL = L([L.ok]);
        fprintf('   against the lottery: cold found %.4f - %.4f d here; the ladder found %.4f d.\n', ...
                min([goodL.tf])*tStar/86400, max([goodL.tf])*tStar/86400, C.tfDays);
        fprintf('   That difference is the whole argument for continuation.\n');
    end
else
    fprintf('   The ladder stopped short of the target engine, so there is no campaign\n');
    fprintf('   number to compare against at %g N / %g s. WHAT THIS RUN DID SHOW: the\n', Tcur, ispCur);
    fprintf('   mechanism, end to end, with a CERTIFIED root at the depth it reached.\n');
    fprintf('   WHAT IT DID NOT: that this cell can be walked to %g N. The campaign\n', op.thrustN);
    fprintf('   reached %g N from a DIFFERENT cell -- the fastest 0.5 N entry of the\n', op.thrustN);
    fprintf('   sheet -- and the wall depth is a property of the cell, not of the code.\n');
end

%% ========================================================================
%  7. THE BASIN HUNT -- the third mechanism, and the one that found four of
%     the five families. Take the root we now hold and ask the DIRECT solver
%     for a transfer to a different arrival phase, warm-started from its PMP
%     flight. The solver does not know which branch its seed is on. A root
%     that comes back FASTER than the family you started from is a new
%     family worth anchoring.
%% ========================================================================
H = struct('sA', {}, 'tf', {}, 'rel', {}, 'perisKm', {}, 'ok', {});
if hunt.run
    fprintf('\n7. THE BASIN HUNT from the %g N root at sA %.4f (t_f %.3f d)\n', ...
            Tcur, phase.sA, day(zCur(8)));
    fprintf('      %9s  %10s  %10s  %11s  %s\n', 'sA', 't_f (d)', 'vs source', 'peris (km)', 'reading');
    [tjH, yjH] = pumpkyn.cr3bp.tfMinProp(zCur(8), [rv0(:); 1; zCur(1:7)], TndR, cCur, muStar);
    sN = linspace(0, tjH(end), dladder.N + 1);
    X0h = interp1(tjH, yjH(:, 1:7), sN, 'spline').';
    LVh = interp1(tjH, yjH(:, 11:13), sN, 'spline');
    U0h = [(-LVh ./ max(vecnorm(LVh, 2, 2), eps)).'; ones(1, dladder.N + 1)];
    dHunt = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', true, ...
                   'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', hunt.cpuSec);
    for sAh = hunt.sA
        rvh = stateA(sAh);  tHh = tic;
        try
            oH = solveDirect(rv0, rvh(1:6), TndR, cCur, muStar, dladder.N, X0h, U0h, zCur(8), dHunt);
            perisH = oH.altMinKm + rMoonKm;
            okH = oH.success && oH.maxDefect < tol.defect && perisH >= num.clearKm;
            rel = (oH.tf - zCur(8))/zCur(8);
            if ~okH,                      rd = tern(oH.success, 'refused on defect or clearance', ...
                                                    ['no solution: ' oH.ipoptStatus]);
            elseif abs(rel) < tol.basin,  rd = 'same basin, continued';
            elseif rel < 0,               rd = 'ANOTHER BASIN, and FASTER -- anchor it';
            else,                         rd = 'another basin, slower';
            end
            fprintf('      %9.4f  %10.3f  %+9.1f%%  %11.0f  %s (%.0f s)\n', ...
                    sAh, day(oH.tf), 100*rel, perisH, rd, toc(tHh));
            H(end+1) = struct('sA', sAh, 'tf', oH.tf, 'rel', rel, 'perisKm', perisH, 'ok', okH);
        catch ME
            fprintf('      %9.4f  THREW: %s\n', sAh, firstline(ME.message));
        end
        S.hunt = H;  saveq(outMat, S);
    end
    nOther = nnz(arrayfun(@(h) h.ok && abs(h.rel) >= tol.basin, H));
    fprintf('   %d of %d probes landed in another basin at a %.0f%% threshold.\n', ...
            nOther, numel(H), 100*tol.basin);
    fprintf('   This is the move that found fast2, direct18 and direct11 (FINDINGS 61, 63, 68):\n');
    fprintf('   solve where the known family looks slow, and read what comes back.\n');
else
    fprintf('\n7. THE BASIN HUNT: skipped (hunt.run = false)\n');
end

%% ------------------------------------------------------------------------
%  SELF-CHECK and HAND-OFF
%
%  B1b, B2 and C1 are NOT machinery gates. A ladder that stops short -- at a
%  basin jump, or at a wall -- and a root that does not certify at the far
%  end are findings about the problem, and the script reports them without
%  throwing. B1 and H1 are machinery: if the direct ladder produced no rung
%  at all, or the harvest cannot be shot to a root, nothing above can be
%  believed and the script fails.
%% ------------------------------------------------------------------------
gateStatus = { ...
    'B1 direct ladder ran',      b1; ...
    'B1b direct ladder finished', b2; ...
    'H1 handoff root',           h1; ...
    'B2 target engine reached',  b3; ...
    'C1 certified',              C.ok};
failed = gateStatus(~[gateStatus{:, 2}], 1);
fprintf('\n   SELF-CHECK: %d of %d gates passed', nnz([gateStatus{:, 2}]), size(gateStatus, 1));
if isempty(failed), fprintf('; a root was built from nothing and certified at the target engine.\n');
else, fprintf('; NOT passed: %s\n', strjoin(failed', ', ')); end
fprintf('   run saved: %s\n', outMat);
fprintf(['   NEXT: anchor_study.m turns a root into a SEED and a new root;\n' ...
         '         transfer_study.m asks whether a root is a minimum;\n' ...
         '         run_phase_torus.m builds a library from anchors like this one.\n']);
assert(b1 && h1, 'root_origins_study:machinery', ...
       'the ladder / harvest / shooting machinery did not work: %s', strjoin(failed', ', '));

%% ------------------------------------------------------------------------
function varargout = fenced(pool, fh, nout, capSec, varargin)
% FENCED  One call under a HARD wall-clock cap when a pool exists, direct
% otherwise. capped_pool returns [] without the Parallel Computing Toolbox,
% and run_capped hands that straight to parfeval, which refuses it.
%
% INPUTS:
%   pool     - parallel pool or [] [object]
%   fh       - the function to call [function_handle]
%   nout     - outputs to request [scalar]
%   capSec   - hard wall-clock cap [scalar, s]
%   varargin - arguments passed to fh
%
% OUTPUTS:
%   varargout - {ok, fh's nout outputs}
varargout = cell(1, nout + 1);
if isempty(pool)
    [varargout{2:nout+1}] = feval(fh, varargin{:});
    varargout{1} = true;
    return
end
[varargout{1}, varargout{2:nout+1}] = run_capped(pool, fh, nout, capSec, varargin{:});
end

function o = solveDirect(rv0, rvf, Tnd, cnd, mu, N, X0, U0, tf0, dopts)
% SOLVEDIRECT  One direct collocation solve, model stripped.
%
% INPUTS:
%   rv0, rvf - departure / arrival states [6x1 or 1x6, ND]
%   Tnd, cnd - ND thrust acceleration at unit mass, ND exhaust speed [scalar]
%   mu       - CR3BP mass ratio [scalar]
%   N        - collocation intervals [scalar]
%   X0, U0   - warm start [7x(N+1)], [4x(N+1)]; EMPTY = a genuine COLD start
%   tf0      - time-of-flight guess [scalar, ND]
%   dopts    - casadi_mintime_dro options [struct]
%
% OUTPUTS:
%   o        - solver output without the CasADi model [struct]
o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tnd, cnd, mu, N, X0, U0, tf0, dopts);
if isfield(o, 'model'), o = rmfield(o, 'model'); end
end

function saveq(f, S)
% SAVEQ  Save the running record, never fatal.
%
% INPUTS:
%   f - output .mat path [char]
%   S - the record [struct]
%
% OUTPUTS:
%   none
try
    if ~isfolder(fileparts(f)), mkdir(fileparts(f)); end
    save(f, '-struct', 'S');
catch
end
end

function s = firstline(m)
% FIRSTLINE  First line of a message, truncated.
%
% INPUTS:
%   m - message [char]
%
% OUTPUTS:
%   s - at most 80 characters of its first line [char]
c = strsplit(m, newline);
s = c{1};
if numel(s) > 80, s = s(1:80); end
end

function s = pass(c)
% PASS  Verdict text.
%
% INPUTS:
%   c - condition [logical]
%
% OUTPUTS:
%   s - 'PASS' or 'FAIL' [char]
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function v = gvd(s, f, d)
% GVD  Struct field with a default.
%
% INPUTS:
%   s - struct
%   f - field name [char]
%   d - default value
%
% OUTPUTS:
%   v - s.(f) if present and non-empty, else d
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = tern(c, a, b)
% TERN  Ternary selection.
%
% INPUTS:
%   c - condition [logical]
%   a - value if true
%   b - value if false
%
% OUTPUTS:
%   s - a or b
if c, s = a; else, s = b; end
end
