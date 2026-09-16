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
    'tfExps', [1.0 1.3 1.7 2.2 3.0], ...   % per-rung t_f guess sweep; the fastest
    'attemptSec', 600, ...                 %   converged candidate is kept (a search)
    'maxPropFrac', 0.6, ...                % seed policy: skip guesses whose all-burn
    'maxMiss', 2);                         %   mass row would burn more than this
                                           % attemptSec: per-rung budget over the sweep
                                           % maxMiss: consecutive empty rungs = stop

%% C. the basin hunt at whatever operating point we reached
%     A hunt probe is a lottery ticket, not a required result, so it gets a
%     short budget: a probe that has not converged in a few minutes is
%     telling you there is no nearby basin, and the next phase is cheaper to
%     try than this one is to finish.
hunt = struct( ...
    'run',    true, ...
    'cpuSec', 300, ...     % per-probe IPOPT budget
    'sA',     [0.0754 + 6/24, 0.0754 + 12/24, 0.0754 + 17/24]);   % arrival phases to probe

%% budgets. clearKm, gateKm and op.m0kg are FORWARDED to certify_root in
%  section 6: that function reads them from its own options (defaults 1900,
%  100, 150) and does NOT take them from the problem struct, so editing them
%  here without forwarding would silently leave certification on defaults.
num = struct('maxIter', 3000, 'maxCpuSec', 900, 'coldCpuSec', 1800, ...
             'clearKm', 1900, 'K', 24, 'shootSec', 240, 'certSec', 900, ...
             'gateKm', 100, 'gateVms', 10);

%% tolerances. ACCEPTANCE thresholds and SOLVE targets are different things
%  and are kept apart: a solver is ASKED for more than the gate requires, so
%  a solve that plateaus above its own target can still be inside the gate.
tol = struct( ...
    'closure',   1e-7,  ...  % orbit periodicity
    'seamDeriv', 1e-6,  ...  % interpolant derivative mismatch across s = 0
    'defect',    1e-8,  ...  % collocation defect on an accepted rung
    'interp',    1e-8,  ...  % Hermite-Simpson interpolation residual (NOT in maxDefect)
    'tfSpread',  1e-8,  ...  % lifted-time continuity across the nodes
    'termErr',   1e-7,  ...  % terminal boundary residual
    'unit',      1e-8,  ...  % |u| = 1 on the solved directions
    'thrMin',    1e-6,  ...  % 1 - throttle: the throttle is FREE, so saturation is CHECKED
    'R',         1e-8,  ...  % multiple-shooting residual, inf-norm (acceptance)
    'Rsolve',    1e-11, ...  % what the shooting solver is ASKED for (tighter than the gate)
    'tfCluster', 0.02,  ...  % ND-relative spread within which two t_f are one cluster
    'same',      1e-2);      % |t_f - the reference's|/t_f to call the flight times equal

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
%     Both orbits are generated from their own parameters (get_family_orbit
%     starts from a catalogued family seed and refines it to closure), so
%     the run is reproducible from section 0 alone. No TRANSFER is loaded
%     anywhere in this script; that is the point of it.
%% ========================================================================
[tD, rvD, infoD] = get_family_orbit(dep.family, struct('tau', dep.tau, 'muStar', muStar));
[tT, rvT, infoT] = get_family_orbit(arr.family, struct('Np', arr.Np, 'pm', arr.pm, 'muStar', muStar));
closD = norm(rvD(end, :) - rvD(1, :));   closA = norm(rvT(end, :) - rvT(1, :));
[stateD, seamD] = phase_state(tD, rvD);
[stateA, seamA] = phase_state(tT, rvT);
rv0 = stateD(phase.sD);   rv0 = rv0(1:6);
rvf = stateA(phase.sA);
orbOK = closD < tol.closure && closA < tol.closure && max(seamD.deriv, seamA.deriv) < tol.seamDeriv;
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
%     The recorded result on this exact cell (2026-08-03, FINDINGS). Those
%     numbers were logged as ALTITUDE above the lunar surface; the table this
%     section prints is Moon-CENTRE radius, so they are converted here (+1737.4
%     km) and both are labelled. A radius cannot be negative; an altitude can,
%     and a negative one means the arc went through the Moon.
%        N = 400  -> 4.3807 ND, +9.1%,  alt 4819 km  (radius 6556 km)
%        N = 800  -> 4.4507 ND, +10.8%, alt -343 km  (radius 1394 km, THROUGH the Moon)
%        N = 1600 -> 4.7833 ND, +19.1%, alt 4680 km  (radius 6417 km), and a
%                    7.8 m WORST-INTERVAL position error -- an accuracy
%                    measure, not an endpoint error and not a minimality
%                    certificate
%     The N = 1600 row is the decisive one: an accurate, lunar-safe,
%     minimum-time extremal that is 19% slower than the reference. The
%     accuracy machinery works cold. BASIN SELECTION DOES NOT.
%% ========================================================================
dCold = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', false, ...
               'returnModel', true, 'maxCpuSec', num.coldCpuSec);
if cold.floor, dCold.minAltKm = floorKm; end
L = struct('N', {}, 'tf', {}, 'defect', {}, 'nodeRadiusKm', {}, 'rel', {}, 'ok', {}, ...
           'status', {}, 'why', {}, 'wall', {});
outcomeA = 'NOT RUN';
if cold.run
    fprintf('\n2. THE COLD LOTTERY at %.0f mN / Isp %g s -- one solve per mesh, same problem\n', ...
            op.thrustN*1000, op.ispS);
    fprintf('   (the "peris" column is the smallest Moon-centre radius over the NODES:\n');
    fprintf('    a screen, not a certified periselene -- the arc between nodes is unchecked)\n');
    fprintf('      %6s  %10s  %9s  %8s  %10s  %11s  %s\n', ...
            'N', 't_f (ND)', 't_f (d)', 'vs ref', 'defect', 'peris (km)', 'status');
    for N = cold.meshes
        tM = tic;
        try
            oM = solveDirect(rv0, rvf, TndOp, cndOp, muStar, N, [], [], cold.tf0, dCold);
            [okM, whyM] = directOK(oM, tol);
            relM = (oM.tf*tStar/86400 - op.tfDays)/op.tfDays;
            L(end+1) = struct('N', N, 'tf', oM.tf, 'defect', oM.maxDefect, ...
                              'nodeRadiusKm', oM.altMinKm + rMoonKm, 'rel', relM, 'ok', okM, ...
                              'status', oM.ipoptStatus, 'why', whyM, 'wall', toc(tM));
            fprintf('      %6d  %10.4f  %9.3f  %+7.1f%%  %10.1e  %11.0f  %s (%.0f s)\n', ...
                    N, oM.tf, day(oM.tf), 100*relM, oM.maxDefect, oM.altMinKm + rMoonKm, ...
                    tern(okM, 'converged', [oM.ipoptStatus ' / ' whyM]), toc(tM));
        catch ME
            L(end+1) = struct('N', N, 'tf', NaN, 'defect', NaN, 'nodeRadiusKm', NaN, ...
                              'rel', NaN, 'ok', false, 'status', 'THREW', 'why', firstline(ME.message), ...
                              'wall', toc(tM));
            fprintf('      %6d  THREW: %s\n', N, firstline(ME.message));
        end
        S.lottery = L;  saveq(outMat, S);
    end
    good = L([L.ok]);
    fprintf('   %d of %d meshes returned a converged, feasible discrete solution; the others are\n', ...
            numel(good), numel(L));
    fprintf('   ITERATES, printed for the record and excluded from what follows.\n');
    if numel(good) >= 2
        % CLUSTERING, stated exactly: two flight times are one cluster when they
        % differ by less than tol.tfCluster of the SHORTEST converged t_f (an
        % absolute ND tolerance, so adding a long outlier cannot regroup the
        % others -- uniquetol's default scales by the largest element). A
        % cluster of flight times is not a basin: two roots can share a t_f, and
        % one continuous family can span several. It is the cheap first read.
        spread   = (max([good.tf]) - min([good.tf]))/min([good.tf]);
        nCluster = numel(uniquetol([good.tf], tol.tfCluster*min([good.tf]), 'DataScale', 1));
        below    = nnz([good.nodeRadiusKm] < num.clearKm);
        agree    = nnz(abs([good.rel]) < tol.same);
        fprintf('   t_f spread %.1f%%; %d flight-time cluster(s) at %.0f%% of the shortest;\n', ...
                100*spread, nCluster, 100*tol.tfCluster);
        fprintf('   %d within %.0f%% of the reference t_f; %d with a node radius under the %g km floor.\n', ...
                agree, 100*tol.same, below, num.clearKm);
        fprintf('   Each converged row is a converged DISCRETE solution of the same NLP family.\n');
        fprintf('   Solver status and defect size are the same for all of them; the flight times\n');
        fprintf('   are not. Nothing in the exit status labels which one you would want.\n');
        if nCluster > 1 && agree == 0
            outcomeA = 'SUPPORTED: mesh density changed the answer, and none matched the reference';
        elseif nCluster > 1
            outcomeA = 'SUPPORTED: mesh density changed the answer';
        else
            outcomeA = 'NOT OBSERVED on these meshes (the recorded 2026-08-03 table did observe it)';
            fprintf('   (These meshes agreed. That is a result about THESE meshes, not a property\n');
            fprintf('    of cold solving -- see the recorded table above.)\n');
        end
    else
        outcomeA = 'INCONCLUSIVE: fewer than two meshes converged';
        fprintf('   fewer than two meshes converged: the lottery cannot be read from this run\n');
    end
    fprintf('   A: %s\n', outcomeA);
else
    fprintf('\n2. THE COLD LOTTERY: skipped (cold.run = false)\n');
end

%% ========================================================================
%  3. THE DIRECT LADDER -- from where the problem is easy toward where we
%     want it. At 15 N the transfer is short (measured below: half a day,
%     under one revolution) and the NLP converges cold in seconds. Each
%     rung below it is a fresh collocation solve warm-started from the rung
%     above, with a time-of-flight guess scaled by T^-exp. The wall times
%     and the guess-to-answer ratios are printed per rung so "starts next
%     to its answer" is a number, not a slogan.
%
%     THE THROTTLE IS FREE AND ITS SATURATION IS CHECKED. The direct solver
%     leaves the throttle as an unknown on purpose (its header: pinning it
%     would be simpler, but checking that it saturates is a genuine test of
%     the formulation). Minimum time is all-burn, so a rung whose throttle
%     dips is not the problem the shooting solver solves later. directOK
%     gates on thrMin for every direct solve in this script.
%% ========================================================================
dLad = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', false, ...
              'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', num.maxCpuSec);
fprintf(['\n3. THE DIRECT LADDER (N = %d, t_f guess scaled by T^-%.1f, clearance floor ON,\n' ...
         '   a rung refused if its t_f leaves [1/%g, %g] x its own guess)\n'], ...
        dladder.N, dladder.tfExp, dladder.tfBand, dladder.tfBand);
fprintf('      %8s %7s  %10s  %9s  %10s  %11s  %8s  %s\n', ...
        'T (N)', 'Isp (s)', 't_f (ND)', 't_f (d)', 'defect', 'peris (km)', 'dV km/s', 'status');
R = struct('thrustN', {}, 'ispS', {}, 'tf', {}, 'defect', {}, 'nodeRadiusKm', {}, 'dv', {}, ...
           'mf', {}, 'thrMin', {}, 'guess', {}, 'band', {}, 'ok', {}, 'why', {}, 'wall', {});
seedX = [];  seedU = [];  seedTf = cold.tf0;  Tprev = [];  oPrev = [];  kAccepted = 0;
for kr = 1:size(dladder.rungs, 1)
    TN = dladder.rungs(kr, 1);  isp = dladder.rungs(kr, 2);
    Tnd = ndT(TN, op.m0kg);  cnd = ndC(isp);
    if ~isempty(Tprev), seedTf = seedTf*(Tprev/TN)^dladder.tfExp; end
    tR = tic;
    try
        oR = solveDirect(rv0, rvf, Tnd, cnd, muStar, dladder.N, seedX, seedU, seedTf, dLad);
        radR = oR.altMinKm + rMoonKm;                     % node-sampled, Moon centre
        [feasR, whyR] = directOK(oR, tol);
        % THE SCREEN. A ladder's premise is that each rung lands NEAR its
        % guess. A rung that lands far from it has, at the least, left the
        % branch the guess was on -- and it can look perfect while doing so:
        % converged, defect 1e-14, clear of the Moon. Comparing the answer
        % with the guess that produced it is the only cheap thing that
        % catches it. It is a plausibility heuristic, not a branch detector:
        % a sharp turn on the same branch can trip it, and a jump to a
        % nearby-time branch can pass it.
        % The band applies to a STEP, never to the top rung: that one's guess
        % is a cold round number with no branch behind it, and at 15 N the
        % true answer is ~40x below it. preflight_screen takes seedTf = []
        % for the same reason.
        band = oR.tf/seedTf;
        inBand = isempty(Tprev) || (band <= dladder.tfBand && band >= 1/dladder.tfBand);
        clearR = radR >= num.clearKm;
        okR = feasR && clearR && inBand;
        dv = cnd*log(1/oR.mf)*lStar/tStar;               % ideal dV from the mass actually used
        % REASONS IN ORDER OF PRECEDENCE: a solve that did not converge is
        % reported as that, never as a branch change.
        if okR,             note = 'accepted';
        elseif ~feasR,      note = ['REFUSED: ' whyR];
        elseif ~clearR,     note = sprintf('REFUSED: node radius %.0f km under the %g km floor', radR, num.clearKm);
        else,               note = sprintf('REFUSED: %.1fx its own guess -- outside the band, possible branch change', band);
        end
        fprintf('      %8.3f %7.0f  %10.4f  %9.3f  %10.1e  %11.0f  %8.4f  %s (%.0f s)\n', ...
                TN, isp, oR.tf, day(oR.tf), oR.maxDefect, radR, dv, note, toc(tR));
        if ~isempty(Tprev)
            fprintf('               (guess %.4f ND from the rung above, answer %.2fx it; 1 - throttle %.1e)\n', ...
                    seedTf, band, 1 - oR.thrMin);
        end
        R(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', oR.tf, 'defect', oR.maxDefect, ...
                          'nodeRadiusKm', radR, 'dv', dv, 'mf', oR.mf, 'thrMin', oR.thrMin, ...
                          'guess', seedTf, 'band', band, 'ok', okR, 'why', note, 'wall', toc(tR));
        S.direct = R;  saveq(outMat, S);
        if ~okR
            fprintf('      the direct ladder STOPS here: this rung cannot seed the next.\n');
            if feasR && clearR
                fprintf('      A band failure is a statement about the STEP, not the problem: the\n');
                fprintf('      remedy is a finer rung ratio, not a looser gate. The shipped\n');
                fprintf('      catalog steps by about 0.75 per rung for exactly this reason.\n');
            end
            break
        end
        kAccepted = kr;
        seedX = oR.X;  seedU = oR.U;  seedTf = oR.tf;  Tprev = TN;  oPrev = oR;
    catch ME
        R(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', NaN, 'defect', NaN, 'nodeRadiusKm', NaN, ...
                          'dv', NaN, 'mf', NaN, 'thrMin', NaN, 'guess', seedTf, 'band', NaN, ...
                          'ok', false, 'why', ['THREW: ' firstline(ME.message)], 'wall', toc(tR));
        S.direct = R;  saveq(outMat, S);
        fprintf('      %8.3f %7.0f  THREW: %s\n', TN, isp, firstline(ME.message));
        break
    end
end
acc = R([R.ok]);
b1 = ~isempty(acc);
assert(b1, 'the direct ladder produced no usable rung: nothing below is meaningful');
b1b = kAccepted == size(dladder.rungs, 1);                % the LAST ROW was accepted, not just its thrust
fprintf('   B1  direct ladder produced %d accepted rung(s)                          %s\n', numel(acc), pass(b1));
fprintf('   B1b direct ladder finished: reached %g N / %g s of %g N / %g s requested   %s\n', ...
        acc(end).thrustN, acc(end).ispS, dladder.rungs(end, 1), dladder.rungs(end, 2), pass(b1b));
if numel(acc) >= 2
    fprintf('   t_f %.3f -> %.3f d as thrust fell %g -> %g N; ideal dV %.4f -> %.4f km/s;\n', ...
            day(acc(1).tf), day(acc(end).tf), acc(1).thrustN, acc(end).thrustN, acc(1).dv, acc(end).dv);
    fprintf('   mass fraction used %.3f -> %.3f. Both directions are MEASURED here, not assumed:\n', ...
            1 - acc(1).mf, 1 - acc(end).mf);
    fprintf('   at fixed exhaust speed the propellant used goes as T*t_f, and if t_f ~ T^-0.6\n');
    fprintf('   then T*t_f ~ T^0.4 FALLS as thrust falls -- lower thrust, longer trip, less dV.\n');
end

%% ========================================================================
%  4. THE HANDOFF -- direct to indirect, which is where the costates enter.
%
%     The multiplier attached to a dynamics-defect constraint measures how
%     much the objective would improve if that constraint were relaxed,
%     which is the definition of the adjoint variable. So the NLP has been
%     computing costates all along under another name. Harvest them onto a
%     shooting mesh and you have a SEED: a state-and-costate trajectory on
%     K+1 junction times, not a vector. Shoot it and you have an indirect
%     root. anchor_study.m opens every step of this box; here we need the
%     root.
%
%     ONE CAVEAT THE READER SHOULD CARRY. The direct rung was solved with a
%     lunar-clearance PATH CONSTRAINT; the shooting solver solves the
%     UNCONSTRAINED PMP. If the floor were active on the harvested arc, its
%     multipliers would not be costates of the unconstrained problem. The
%     margin is printed: a node radius well above the floor is the evidence
%     that the constraint was slack.
%% ========================================================================
TndB = ndT(acc(end).thrustN, op.m0kg);   cndB = ndC(acc(end).ispS);
pool = capped_pool();      % [] without the Parallel Computing Toolbox: then
                           % every solver call below runs UNFENCED and a single
                           % crawling integration can stall the run
% THE BOTTOM RUNG, RE-SOLVED ONCE IN THE SUNDMAN CHART. The ladder runs in
% plain time, as the shipped catalog did. Under a Sundman change of variable
% the defect system carries the time state as well, so its multipliers carry
% ONE MORE ROW: lambda_t, which the minimum principle fixes at +1 when the
% objective is t_f. That row is a free check on the whole duals-to-costates
% mapping -- station association, sign and scale together -- and it simply
% does not exist without Sundman. One extra warm solve buys it.
sunNote = '';
try
    dSun = dLad;  dSun.sundman = true;
    oH8 = solveDirect(rv0, rvf, TndB, cndB, muStar, dladder.N, oPrev.X, oPrev.U, oPrev.tf, dSun);
    [feasS, whyS] = directOK(oH8, tol);
    % the re-solve must land on the SAME solution as the plain-time rung:
    % same t_f (to the band a re-mesh can move it) and the same node radius
    sameTf = abs(oH8.tf - oPrev.tf)/oPrev.tf < tol.same;
    if ~feasS,       sunNote = ['Sundman re-solve refused: ' whyS];  oH8 = [];
    elseif ~sameTf,  sunNote = sprintf('Sundman re-solve moved t_f by %.1f%%: not the same solution', ...
                                       100*abs(oH8.tf - oPrev.tf)/oPrev.tf);  oH8 = [];
    end
catch ME
    sunNote = ['Sundman re-solve THREW: ' firstline(ME.message)];  oH8 = [];
end
oHarv = oPrev;  chart = 'plain time (no lambda_t row)';
if ~isempty(oH8), oHarv = oH8;  chart = 'Sundman (lambda_t exposed)'; end
% the harvest needs the multipliers the solver extracts in a try/catch of
% its own: an empty lamDef is a named failure here, not an index error
assert(isfield(oHarv, 'lamDef') && size(oHarv.lamDef, 1) >= 7 && ~isempty(oHarv.Um), ...
       'root_origins_study:duals', 'the direct solve returned no defect multipliers: nothing to harvest');
[seedB, dgB] = harvest_ms_seed(oHarv, num.K);
tH = tic;
[okB, zB, itB] = fenced(pool, @ms_tfmin, 2, num.shootSec + 90, rv0, rvf(1:6), seedB, TndB, cndB, muStar, ...
                        struct('tolR', tol.Rsolve, 'wallSec', num.shootSec, 'maxIter', 100));
if ~okB, itB = struct('normR', NaN, 'iters', 0, 'converged', false, 'Y', [], 'tGrid', []); zB = nan(8, 1); end
voteB = gvd(dgB, 'voteMargin', NaN);  lamTB = gvd(dgB, 'lamT', NaN);
h1 = isfinite(itB.normR) && itB.normR < tol.R;
h2 = isfinite(voteB) && voteB >= 0.9;                                  % sign vote: 90% of stations agree
h3 = ~isfinite(lamTB) || abs(lamTB - 1) < 1e-2;                        % unchecked is not failed; it is unchecked
fprintf('\n4. THE HANDOFF at %g N / Isp %g s: %d multipliers -> a %dx%d seed -> a root\n', ...
        acc(end).thrustN, acc(end).ispS, size(oHarv.lamDef, 2), size(seedB.Y, 1), size(seedB.Y, 2));
fprintf('   harvested from the %s solve%s\n', chart, tern(isempty(sunNote), '', ['  (' sunNote ')']));
if ~isempty(oH8)
    fprintf('   plain-time rung t_f %.4f d -> Sundman re-solve %.4f d (%.2e relative); node radius %.0f -> %.0f km\n', ...
            day(oPrev.tf), day(oH8.tf), abs(oH8.tf - oPrev.tf)/oPrev.tf, ...
            oPrev.altMinKm + rMoonKm, oH8.altMinKm + rMoonKm);
end
fprintf('   clearance floor slack on the harvested arc: node radius %.0f km vs the %g km floor (%.0f km margin)\n', ...
        oHarv.altMinKm + rMoonKm, num.clearKm, oHarv.altMinKm + rMoonKm - num.clearKm);
fprintf('   H2 sign vote      %5.1f%% of stations agree / 90%%                     %s\n', 100*voteB, pass(h2));
if isfinite(lamTB)
    fprintf('   H3 lambda_t       %9.6f / 1 +- 1e-2                          %s   (%s)\n', lamTB, pass(h3), ...
            tern(gvd(dgB, 'lamTOK', true), 'the mapping''s own check agrees', 'the mapping''s own check DISAGREES'));
    fprintf('       (a constant +1 row checks sign and scale; it cannot see a half-step station\n');
    fprintf('        shift -- midpoint association is a RULE duals_to_costates applies, not a\n');
    fprintf('        thing this number proves)\n');
else
    fprintf('   H3 lambda_t       NOT EXPOSED without Sundman: the scale check is unavailable here\n');
end
fprintf('   H1 shooting |R|   %9.2e / %-9.0e in %d iterations, %.0f s    %s\n', ...
        itB.normR, tol.R, itB.iters, toc(tH), pass(h1));
fprintf('   t_f: direct %.4f d -> shot %.4f d\n', day(oHarv.tf), day(zB(8)));
fprintf('   lambda(0) = [%s]\n', strjoin(compose('%+.6g', zB(1:7)'), ' '));
S.handoff = struct('thrustN', acc(end).thrustN, 'ispS', acc(end).ispS, 'chart', chart, 'sunNote', sunNote, ...
                   'z', zB(:), 'Y', itB.Y, 'tGrid', itB.tGrid, 'normR', itB.normR, ...
                   'voteMargin', voteB, 'lamT', lamTB, 'lamTOK', gvd(dgB, 'lamTOK', NaN), ...
                   'nodeRadiusKm', oHarv.altMinKm + rMoonKm);
saveq(outMat, S);
assert(h1, 'the handoff shooting residual %.2e is above the gate %.0e: there is no root to walk', ...
       itB.normR, tol.R);

%% ========================================================================
%  5. THE INDIRECT LADDER -- the rest of the way down.
%
%     Below about 0.5 N the collocation NLP loses the thread on this problem
%     (an empirical observation from the shipped catalogs, not a theorem):
%     the arc winds more, the mesh has to resolve more revolutions, and the
%     basin gets narrow. Multiple shooting suffers less, because each of its
%     K segments is short. So the walk continues in the indirect chart:
%     each rung rebuilds the PREVIOUS rung's trajectory FROM ITS BANKED
%     JUNCTIONS, segment by segment, re-cuts it onto a new time-of-flight,
%     rebuilds the mass row from the all-burn law for the NEW engine, and
%     solves ms_tfmin. Segment-wise reconstruction is the reason junctions
%     are banked at all: re-flying the whole arc from lambda(0) alone would
%     bring back the amplification that multiple shooting exists to avoid.
%
%     The per-rung t_f GUESS is swept over several exponents and the
%     FASTEST converged candidate is kept. Say what that is: a small SEARCH
%     over guesses at each rung, not branch-preserving continuation. It is
%     the campaign's own rule (probe_deep_rungs), it is what got the
%     recorded chain past a winding wall, and it can in principle switch
%     branches -- which is why the turn count is printed beside each rung.
%
%     A refused rung is a finding about THIS SEARCH POLICY on this problem
%     within this budget, not a bug and not, by itself, a proof that no
%     root exists there. Every attempt is recorded with its reason.
%% ========================================================================
zCur = zB;  Tcur = acc(end).thrustN;  cCur = cndB;  ispCur = acc(end).ispS;  itCur = itB;
W = struct('thrustN', {}, 'ispS', {}, 'tf', {}, 'normR', {}, 'flyKm', {}, 'flyVms', {}, ...
           'turns', {}, 'tfExp', {}, 'ok', {}, 'why', {}, 'tried', {}, 'wall', {}, ...
           'z', {}, 'Y', {}, 'tGrid', {});
if iladder.run
    fprintf('\n5. THE INDIRECT LADDER (multiple shooting, K = %d segments)\n', num.K);
    fprintf('   ("turns" = total swept Moon-centred azimuth in the rotating frame / 2 pi:\n');
    fprintf('    reversals add, so it is an angular-excursion diagnostic, not a winding number)\n');
    fprintf('      %8s %7s  %10s  %9s  %10s  %9s  %6s  %s\n', ...
            'T (N)', 'Isp (s)', 't_f (ND)', 't_f (d)', '|R|', 'flown km', 'turns', 'status');
    nMiss = 0;
    for kr = 1:size(iladder.rungs, 1)
        TN = iladder.rungs(kr, 1);  isp = iladder.rungs(kr, 2);
        Tnd = ndT(TN, op.m0kg);  cnd = ndC(isp);
        tR = tic;  best = [];
        tried = struct('timeout', 0, 'notConverged', 0, 'flownMiss', 0, 'propSkip', 0, 'bad', 0);
        % the previous root, rebuilt from its banked junctions (fenced: a
        % segment that parks near the Moon can crawl)
        [okFly, tj, yj] = fenced(pool, @flyFromJunctions, 2, 120, itCur, rv0, ndT(Tcur, op.m0kg), cCur, muStar);
        if ~okFly
            fprintf('      %8.3f %7.0f  the previous root could not be re-flown within 120 s: stopping\n', TN, isp);
            break
        end
        ispOnly = TN == Tcur;
        expList = iladder.tfExps;  if ispOnly, expList = 1.0; end   % the exponent has no effect on an Isp step
        for tfExp = expList
            if toc(tR) > iladder.attemptSec, break, end  % this rung has had its budget
            if ispOnly, tfGuess = zCur(8);
            else,       tfGuess = zCur(8)*(Tcur/TN)^tfExp;
            end
            % SEED PROPELLANT POLICY: the all-burn mass row is m = 1 - Tnd*t/cnd,
            % so a guess longer than maxPropFrac*cnd/Tnd would seed a vehicle
            % that has burned more than maxPropFrac of itself. That is a
            % choice about which seeds to try, not a physical exhaustion bound
            % (exhaustion is at cnd/Tnd), and a skip is recorded as a skip.
            if tfGuess > iladder.maxPropFrac*cnd/Tnd, tried.propSkip = tried.propSkip + 1; continue, end
            [Yg, tGs] = flight_to_junctions(tj, yj, num.K, ...
                            struct('tf', tfGuess, 'massLaw', struct('Tnd', Tnd, 'cnd', cnd)));
            capLeft = max(30, iladder.attemptSec - toc(tR));
            [okRun, zt, it] = fenced(pool, @ms_tfmin, 2, min(num.shootSec, capLeft) + 90, ...
                rv0, rvf(1:6), struct('tf', tfGuess, 'tGrid', tGs, 'Y', Yg), ...
                Tnd, cnd, muStar, struct('tolR', tol.Rsolve, 'wallSec', min(num.shootSec, capLeft)));
            if ~okRun,                       tried.timeout = tried.timeout + 1;  continue, end
            if ~all(isfinite(zt)) || ~(isfinite(it.normR) && it.normR < tol.R)
                tried.notConverged = tried.notConverged + 1;  continue
            end
            % the flown witness: from z8 alone, must REACH t_f with a positive
            % mass and land in position AND velocity (fail-closed on NaN)
            [okW, tw, yFly] = fenced(pool, @(a, b, c, d, e) pumpkyn.cr3bp.tfMinProp(a, b, c, d, e), 2, 120, ...
                                     zt(8), [rv0(:); 1; zt(1:7)], Tnd, cnd, muStar);
            reached = okW && ~isempty(tw) && all(isfinite(yFly(end, :))) && ...
                      abs(tw(end) - zt(8)) < 1e-9*max(1, zt(8)) && yFly(end, 7) > 0;
            if ~reached, tried.bad = tried.bad + 1; continue, end
            mk  = norm(yFly(end, 1:3).' - rvf(1:3))*lStar;
            mv  = norm(yFly(end, 4:6).' - rvf(4:6))*lStar/tStar*1000;
            if ~(isfinite(mk) && mk < num.gateKm && isfinite(mv) && mv < num.gateVms)
                tried.flownMiss = tried.flownMiss + 1;  continue
            end
            if isempty(best) || zt(8) < best.z(8)
                best = struct('z', zt, 'it', it, 'mk', mk, 'mv', mv, 'tfExp', tfExp, 'yFly', yFly);
            end
        end
        why = sprintf('timeout %d, unconverged %d, flown-miss %d, seed-policy skip %d, bad flight %d', ...
                      tried.timeout, tried.notConverged, tried.flownMiss, tried.propSkip, tried.bad);
        if isempty(best)
            nMiss = nMiss + 1;
            W(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', NaN, 'normR', NaN, 'flyKm', NaN, 'flyVms', NaN, ...
                              'turns', NaN, 'tfExp', NaN, 'ok', false, 'why', why, 'tried', tried, ...
                              'wall', toc(tR), 'z', [], 'Y', [], 'tGrid', []);
            S.indirect = W;  saveq(outMat, S);
            fprintf('      %8.3f %7.0f  %10s  %9s  %10s  %9s  %6s  no converged candidate in %.0f s [%s]\n', ...
                    TN, isp, '-', '-', '-', '-', '-', toc(tR), why);
            if nMiss >= iladder.maxMiss
                fprintf('      THE WALK STOPS: %d consecutive rungs produced no converged candidate.\n', nMiss);
                fprintf('      That is a statement about this search policy -- these rung ratios,\n');
                fprintf('      these exponents, this budget, this seed-propellant cap -- on this\n');
                fprintf('      cell. It is not by itself a proof that no root exists below. The\n');
                fprintf('      recorded chain stopped at 0.067 N walking from the sheet''s FASTEST\n');
                fprintf('      0.5 N cell; where the deep walk starts is a choice, not a detail.\n');
                fprintf('      The ways on are a finer rung ratio, arclength steps in thrust, a\n');
                fprintf('      wider exponent sweep or budget, or a different ROUTE (the Isp stage).\n');
                break
            end
            continue
        end
        ang = unwrap(atan2(best.yFly(:, 2), best.yFly(:, 1) - (1 - muStar)));
        turns = sum(abs(diff(ang)))/(2*pi);
        fprintf('      %8.3f %7.0f  %10.4f  %9.3f  %10.1e  %9.1f  %6.2f  accepted, tfExp %.1f (%.0f s)\n', ...
                TN, isp, best.z(8), day(best.z(8)), best.it.normR, best.mk, turns, best.tfExp, toc(tR));
        % THE ROOT ITSELF IS KEPT, not just its t_f: bare z8 pins t_f only to
        % ~1e-4 ND at many revolutions, so a rung is banked with its K
        % junction STARTS (ms_bvp's info.Y; the K+1-th column is the flown
        % endpoint and is reconstructed, never stored -- the identifiability
        % rule). "accepted" here means: converged, flown witness inside the
        % position and velocity gates, mass positive. It is a usable
        % continuation seed; physical validation is section 6's job.
        W(end+1) = struct('thrustN', TN, 'ispS', isp, 'tf', best.z(8), 'normR', best.it.normR, ...
                          'flyKm', best.mk, 'flyVms', best.mv, 'turns', turns, 'tfExp', best.tfExp, ...
                          'ok', true, 'why', why, 'tried', tried, 'wall', toc(tR), ...
                          'z', best.z(:), 'Y', best.it.Y, 'tGrid', best.it.tGrid);
        zCur = best.z;  itCur = best.it;  Tcur = TN;  cCur = cnd;  ispCur = isp;  nMiss = 0;
        S.indirect = W;  saveq(outMat, S);
    end
    nAccI = nnz([W.ok]);
    fprintf('   %d of %d indirect rungs attempted were accepted.\n', nAccI, numel(W));
else
    fprintf('\n5. THE INDIRECT LADDER: skipped (iladder.run = false)\n');
end
b2 = abs(Tcur - op.thrustN)/op.thrustN < 1e-6 && abs(ispCur - op.ispS) < 1e-6;
if b2,               outcomeB = 'SUPPORTED: the ladder reached the target engine from a cold 15 N start';
elseif iladder.run,  outcomeB = sprintf('PARTIAL: the ladder reached %g N / %g s of %g N / %g s', Tcur, ispCur, op.thrustN, op.ispS);
else,                outcomeB = sprintf('PARTIAL: direct ladder only, %g N / %g s', Tcur, ispCur);
end
fprintf('   B2 the ladder reached %g N / Isp %g s; the target was %g N / %g s    %s\n', ...
        Tcur, ispCur, op.thrustN, op.ispS, pass(b2));
reached = struct('thrustN', Tcur, 'ispS', ispCur, 'tfDays', day(zCur(8)), 'z', zCur(:), ...
                 'Y', itCur.Y, 'tGrid', itCur.tGrid, 'polished', false);
S.reached = reached;  saveq(outMat, S);

%% ========================================================================
%  6. THE ROOT WE HOLD -- certify it, and compare with the reference.
%     The certifier polishes the root, flies it, runs the pointwise PMP
%     checks, the foreign-solver witness, the conjugate test and the
%     hypothesis gates. transfer_study.m opens that box gate by gate; here
%     we take its verdict, forward every setting it reads from ITS options,
%     and -- if it passes -- adopt ITS polished root for everything after.
%% ========================================================================
TndR = ndT(Tcur, op.m0kg);
seedR = struct('tf', zCur(8), 'tGrid', itCur.tGrid, 'Y', itCur.Y);    % K starts: what ms_bvp reads
B = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, ...
                             'thrustN', Tcur, 'ispS', ispCur, 'm0kg', op.m0kg, ...
                             'tauDRO', dep.tau, 'NpTulip', arr.Np, 'pmTulip', arr.pm), ...
           'Tnd', TndR, 'cnd', cCur, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
certOpts = struct('pool', pool, 'wallSec', num.certSec, 'sA', phase.sA, 'sD', phase.sD, ...
                  'm0kg', op.m0kg, 'gateKm', num.gateKm, 'gateVms', num.gateVms, 'moonKmMin', num.clearKm);
C = certify_root(seedR, rv0, rvf, B, certOpts);
S.certificate = C;  S.certOpts = certOpts;  saveq(outMat, S);
fprintf('\n6. CERTIFY the %g N / Isp %g s root (m0 %g kg, gates %g km / %g m/s, floor %g km): %s\n', ...
        Tcur, ispCur, op.m0kg, num.gateKm, num.gateVms, num.clearKm, C.reason);
if C.ok
    fprintf('   t_f %.4f d, flown miss %.3f km / %.3f m/s, conjugate %s, lift margin %.1fx\n', ...
            C.tfDays, C.flyKm, C.flyVms, tern(gvd(C, 'conj', 0) == 1, 'PASS', 'not passed'), ...
            gvd(C, 'liftMargin', NaN));
    % ADOPT THE POLISHED ROOT: the certificate belongs to C.z / C.Y, not to
    % the pre-polish candidate, so the hunt and the record use those.
    zCur = C.z(:);
    if isfield(C, 'Y') && ~isempty(C.Y), itCur.Y = C.Y; end
    reached.z = zCur;  reached.Y = itCur.Y;  reached.polished = true;  reached.tfDays = C.tfDays;
    S.reached = reached;  saveq(outMat, S);
    if b2
        relRef = abs(C.tfDays - op.tfDays)/op.tfDays;
        fprintf('   this run %.4f d vs the reference''s %.4f d at this cell: %.2f%% apart -- %s\n', ...
                C.tfDays, op.tfDays, 100*relRef, ...
                tern(relRef < tol.same, 'flight times AGREE within 1%', 'flight times DIFFER'));
        fprintf('   (a flight-time match is not a root identity; that needs the states and costates\n');
        fprintf('    compared, or a continuation connecting the two -- anchor_study.m does the first)\n');
        if cold.run && ~isempty(L) && any([L.ok])
            goodL = L([L.ok]);
            fprintf('   against the lottery: cold found %.4f - %.4f d here; the ladder found %.4f d.\n', ...
                    min([goodL.tf])*tStar/86400, max([goodL.tf])*tStar/86400, C.tfDays);
        end
    else
        fprintf('   The ladder stopped short of the target engine, so the reference t_f does not\n');
        fprintf('   apply at %g N / %g s. WHAT THIS RUN DID SHOW: the mechanism, end to end, with\n', Tcur, ispCur);
        fprintf('   a certified root at the depth it reached. WHAT IT DID NOT: that this cell can\n');
        fprintf('   be walked to %g N under this search policy.\n', op.thrustN);
    end
else
    fprintf('   NOT CERTIFIED. The candidate at %g N / %g s is a converged shooting root that\n', Tcur, ispCur);
    fprintf('   did not pass the gate stack (%s). No family or reference comparison is made\n', C.reason);
    fprintf('   from an uncertified candidate; the hunt below starts from it anyway, because a\n');
    fprintf('   hunt needs only a warm start, and says so.\n');
end

%% ========================================================================
%  7. THE BASIN HUNT -- the third mechanism, and the one that found four of
%     the five families. Take the root we hold and ask the DIRECT solver for
%     a transfer to a DIFFERENT arrival phase, warm-started from its PMP
%     flight. The solver does not know which branch its seed is on.
%
%     WHAT THIS DOES AND DOES NOT MEASURE. Each probe is a different
%     boundary-value problem (the arrival state moved), so its t_f is not
%     comparable to the source's as a basin test: one family's t_f varies
%     with phase by far more than a few percent, and two families can share
%     a t_f. What is printed is "a converged transfer at phase sA, X% faster
%     or slower than the source root at ITS phase" -- a candidate. Turning
%     a candidate into a basin verdict needs a SAME-PHASE baseline: continue
%     the held family to sA (arclength_arrival) and compare the two roots'
%     states and costates there. That is what the campaign did (FINDINGS
%     61, 63, 68). Each candidate is banked as a harvested seed so it can be.
%% ========================================================================
H = struct('sA', {}, 'tf', {}, 'rel', {}, 'nodeRadiusKm', {}, 'ok', {}, 'status', {}, 'why', {}, ...
           'seed', {}, 'wall', {});
outcomeC = 'NOT RUN';
if hunt.run
    fprintf('\n7. THE BASIN HUNT from the %g N root at sA %.4f (t_f %.3f d%s)\n', ...
            Tcur, phase.sA, day(zCur(8)), tern(C.ok, ', certified', ', NOT certified: a warm start only'));
    fprintf('      %9s  %10s  %10s  %11s  %s\n', 'sA', 't_f (d)', 'vs source', 'peris (km)', 'reading');
    dHunt = struct('maxIter', num.maxIter, 'scheme', 'hermite-simpson', 'sundman', true, ...
                   'returnModel', true, 'minAltKm', floorKm, 'maxCpuSec', hunt.cpuSec);
    X0h = [];  U0h = [];  huntSeedOK = false;
    try
        [okFly, tjH, yjH] = fenced(pool, @flyFromJunctions, 2, 120, itCur, rv0, TndR, cCur, muStar);
        assert(okFly, 'the source root could not be re-flown within 120 s');
        [tjH, iu] = unique(tjH(:), 'stable');  yjH = yjH(iu, :);    % the propagator can repeat a time sample
        sN  = linspace(0, tjH(end), dladder.N + 1);
        X0h = interp1(tjH, yjH(:, 1:7), sN, 'spline').';
        LVh = interp1(tjH, yjH(:, 11:13), sN, 'spline');
        nLV = vecnorm(LVh, 2, 2);
        assert(all(nLV > 1e-12), 'the primer vanishes on the source flight: no direction to seed');
        U0h = [(-LVh ./ nLV).'; ones(1, dladder.N + 1)];
        huntSeedOK = true;
    catch ME
        fprintf('   the hunt''s warm start could not be built: %s\n', firstline(ME.message));
    end
    for sAh = hunt.sA
        if ~huntSeedOK, break, end
        rvh = stateA(sAh);  tHh = tic;
        try
            oH = solveDirect(rv0, rvh(1:6), TndR, cCur, muStar, dladder.N, X0h, U0h, zCur(8), dHunt);
            radH = oH.altMinKm + rMoonKm;
            [feasH, whyH] = directOK(oH, tol);
            okH = feasH && radH >= num.clearKm;
            rel = (oH.tf - zCur(8))/zCur(8);
            seedH = [];
            if okH && isfield(oH, 'lamDef') && size(oH.lamDef, 1) >= 7, seedH = harvest_ms_seed(oH, num.K); end
            if ~okH,          rd = tern(oH.success, ['refused: ' whyH], ['no converged candidate: ' oH.ipoptStatus]);
            elseif rel < 0,   rd = sprintf('candidate, %.1f%% FASTER than the source at its own phase -- worth a same-phase baseline', -100*rel);
            else,             rd = sprintf('candidate, %.1f%% slower than the source at its own phase', 100*rel);
            end
            fprintf('      %9.4f  %10.3f  %+9.1f%%  %11.0f  %s (%.0f s)\n', ...
                    sAh, day(oH.tf), 100*rel, radH, rd, toc(tHh));
            H(end+1) = struct('sA', sAh, 'tf', oH.tf, 'rel', rel, 'nodeRadiusKm', radH, 'ok', okH, ...
                              'status', oH.ipoptStatus, 'why', rd, 'seed', seedH, 'wall', toc(tHh));
        catch ME
            fprintf('      %9.4f  THREW: %s\n', sAh, firstline(ME.message));
            H(end+1) = struct('sA', sAh, 'tf', NaN, 'rel', NaN, 'nodeRadiusKm', NaN, 'ok', false, ...
                              'status', 'THREW', 'why', firstline(ME.message), 'seed', [], 'wall', toc(tHh));
        end
        S.hunt = H;  saveq(outMat, S);
    end
    nConv = nnz([H.ok]);  nFast = nnz(arrayfun(@(h) h.ok && h.rel < 0, H));
    fprintf('   %d probes attempted, %d converged candidates, %d of them faster than the source.\n', ...
            numel(H), nConv, nFast);
    fprintf('   Each converged candidate is banked as a harvested seed (S.hunt(k).seed): shoot it,\n');
    fprintf('   certify it, and continue the source family to the same phase to learn whether it\n');
    fprintf('   is the same family -- the move that found fast2, direct18 and direct11.\n');
    if ~huntSeedOK,     outcomeC = 'INCONCLUSIVE: no warm start could be built';
    elseif nConv == 0,  outcomeC = 'NOT OBSERVED: no probe converged within its budget';
    else,               outcomeC = sprintf('CANDIDATES: %d converged at other phases (%d faster); basin identity NOT established here', nConv, nFast);
    end
    fprintf('   C: %s\n', outcomeC);
else
    fprintf('\n7. THE BASIN HUNT: skipped (hunt.run = false)\n');
end

%% ------------------------------------------------------------------------
%  SELF-CHECK and HAND-OFF
%
%  Two kinds of line. GATES are machinery: B1 (the direct ladder produced a
%  rung) and H1 (the harvest shot to a root) must pass or nothing above can
%  be believed, and the script fails. B1b, H2, H3, B2 and C1 are reported,
%  never fatal: a ladder that stops short, an unchecked lambda_t and a root
%  that does not certify at the far end are findings. OUTCOMES are what the
%  three experiments showed, in their own words, because a search that
%  found nothing is not a machinery failure.
%% ------------------------------------------------------------------------
gateStatus = { ...
    'B1  direct ladder produced a rung', b1; ...
    'B1b direct ladder finished',        b1b; ...
    'H1  handoff shot to a root',        h1; ...
    'H2  harvest sign vote',             h2; ...
    'H3  harvest lambda_t',              h3; ...
    'B2  target engine reached',         b2; ...
    'C1  certified',                     C.ok};
outcomes = {'A cold lottery', outcomeA; 'B ladder', outcomeB; 'C basin hunt', outcomeC};
S.gateStatus = gateStatus;  S.outcomes = outcomes;
[savedOK, saveMsg] = saveq(outMat, S);
failed = gateStatus(~[gateStatus{:, 2}], 1);
fprintf('\n   SELF-CHECK: %d of %d gates passed', nnz([gateStatus{:, 2}]), size(gateStatus, 1));
if isempty(failed), fprintf('.\n'); else, fprintf('; NOT passed: %s\n', strjoin(strtrim(failed'), ', ')); end
for k = 1:size(outcomes, 1), fprintf('   %-15s %s\n', outcomes{k, 1}, outcomes{k, 2}); end
fprintf('   record: %s\n', tern(savedOK, outMat, ['NOT SAVED -- ' saveMsg]));
fprintf(['   NEXT: anchor_study.m turns a root into a SEED and a new root;\n' ...
         '         transfer_study.m asks whether a root is a minimum;\n' ...
         '         run_phase_torus.m builds a library from anchors like this one.\n']);
assert(b1 && h1, 'root_origins_study:machinery', ...
       'the ladder / harvest / shooting machinery did not work: %s', strjoin(strtrim(failed'), ', '));

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
rv0 = rv0(:).';  rvf = rvf(:).';                  % the solver documents ROW endpoints
o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tnd, cnd, mu, N, X0, U0, tf0, dopts);
if isfield(o, 'model'), o = rmfield(o, 'model'); end
end

function [ok, why] = directOK(o, tol)
% DIRECTOK  The one feasibility predicate every direct solve in this script
% is judged by, so no section accepts less than another.
%
% INPUTS:
%   o   - casadi_mintime_dro output [struct]
%   tol - the script's tolerance block [struct]
%
% OUTPUTS:
%   ok  - every check below passed [logical]
%   why - the FIRST failed check by name, '' if none [char]
checks = { ...
  'not converged',          ~(isfield(o, 'success') && isscalar(o.success) && o.success); ...
  'non-finite fields',      ~(isfinite(o.tf) && isfinite(o.mf) && isfinite(o.maxDefect)); ...
  'collocation defect',     ~(o.maxDefect < tol.defect); ...
  'interpolation residual', ~(gvd(o, 'maxInterp', 0) < tol.interp); ...
  'lifted t_f spread',      ~(gvd(o, 'tfSpread', 0) < tol.tfSpread); ...
  'terminal residual',      ~(o.termErr < tol.termErr); ...
  '|u| = 1',                ~(o.maxUnit < tol.unit); ...
  'throttle not saturated', ~(1 - gvd(o, 'thrMin', 0) < tol.thrMin); ...
  'mass fraction',          ~(o.mf > 0 && o.mf <= 1); ...
  't_f not positive',       ~(o.tf > 0)};
bad = find([checks{:, 2}], 1);
ok = isempty(bad);
if ok, why = ''; else, why = checks{bad, 1}; end
end

function [tj, yj] = flyFromJunctions(it, rv0, Tnd, cnd, mu)
% FLYFROMJUNCTIONS  A shooting root's trajectory rebuilt SEGMENT BY SEGMENT
% from its banked junction starts, so no segment is longer than 1/K of the
% arc and nothing is amplified across the whole flight.
%
% INPUTS:
%   it   - ms_tfmin info: .Y [14 x K] junction starts, .tGrid [1 x K+1] [struct]
%   rv0  - departure state [6x1]; column 1's state part is fixed to it
%   Tnd  - ND thrust acceleration at unit mass [scalar]
%   cnd  - ND exhaust speed [scalar]
%   mu   - CR3BP mass ratio [scalar]
%
% OUTPUTS:
%   tj   - sample times over the whole arc [n x 1]
%   yj   - state and costate samples [n x 14]
K = numel(it.tGrid) - 1;
Y = it.Y;  Y(1:7, 1) = [rv0(:); 1];
tj = [];  yj = [];
for k = 1:K
    dt = it.tGrid(k+1) - it.tGrid(k);
    [ts, ys] = pumpkyn.cr3bp.tfMinProp(dt, Y(:, k), Tnd, cnd, mu);
    if k > 1, ts = ts(2:end); ys = ys(2:end, :); end
    tj = [tj; it.tGrid(k) + ts(:)];
    yj = [yj; ys];
end
end

function [ok, msg] = saveq(f, S)
% SAVEQ  Save the running record. Never fatal, never silent: a failed save
% warns, and the caller can print the last checkpoint's status.
%
% INPUTS:
%   f - output .mat path [char]
%   S - the record [struct]
%
% OUTPUTS:
%   ok  - the save succeeded [logical]
%   msg - the error message when it did not, '' otherwise [char]
ok = true;  msg = '';
try
    if ~isfolder(fileparts(f)), mkdir(fileparts(f)); end
    save(f, '-struct', 'S');
catch ME
    ok = false;  msg = ME.message;
    warning('root_origins_study:save', 'checkpoint NOT saved (%s): %s', f, ME.message);
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
