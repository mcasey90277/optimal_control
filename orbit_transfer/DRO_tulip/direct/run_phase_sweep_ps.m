function run_phase_sweep_ps()
%% Purpose:
%
%  This routine maps the minimum-time DRO -> Tulip transfer over ORBIT PHASING.
%  A satellite can depart from any point on the DRO, and arrive at any point
%  on the Tulip orbit. One number describes each choice: the PHASE, the
%  fraction of the orbit period elapsed past a reference point. Because both
%  phases wrap around at 1, the space of all (departure, arrival) choices is a
%  TORUS. Every point of the torus is a complete minimum-time optimal control
%  problem; this routine solves them one by one and colors the torus with the
%  results.
%
%  The filling strategy has two waves:
%
%    WAVE 0 ("knock on every door"):  every grid point is attempted from a
%    cheap straight-line seed under a small iteration budget. Points where the
%    optimizer converges quickly turn green; stubborn ones fail fast and cheap.
%
%    WAVE 1 ("spread from the doors that opened"):  every green point hands
%    its solved trajectory to its four torus neighbours as a warm start, and
%    the solutions flood outward. The flood stops on its own at "creases" --
%    curves on the torus where the family of optimal transfers reorganizes and
%    no warm start survives. Those creases are findings, not failures.
%
%  A point is called GREEN only if the solution passes a physical test, never
%  on the optimizer's word alone: the solved control history is handed to a
%  high-accuracy integrator and flown once, end to end. If the spacecraft
%  actually arrives (within the map tolerance), the trajectory is real. If
%  not, the point stays red no matter what the solver reported -- a coarse
%  mesh can satisfy its own discrete equations with a "trajectory" no
%  spacecraft could fly.
%
%  Companion routines:
%    run_dro_tulip_ps.m       one transfer, this same style
%    sweep_phasing_direct.m   the industrial version (checkpoints, options)
%    viz/plot_phase_torus.m   the torus pictures this routine ends with
%
%  References:
%  pumpkynPie/demos/lowThrustDRO2Tulip.m   (the style this file mirrors)
%  ../FINDINGS.md                          (why every design choice is here)

%% Constants:
  muStar = 0.012150585609624;          % Mass ratio
   lStar = 389703.264829278;           % Characteristic length (km)
   tStar = 382981.289129055;           % Characteristic time (s)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'lib'), fullfile(here,'certify'), fullfile(here,'viz'));
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));

%% Sweep Settings:
      nD = 6;                          % departure-phase grid points
      nA = 6;                          % arrival-phase grid points
       N = 800;                        % collocation intervals per solve
 floorKm = 500;                        % minimum lunar altitude (km) -- without
                                       % it the problem lets trajectories pass
                                       % through the Moon
wave0Its = 500;                        % iteration budget, wave 0 (fail fast)
wave1Its = 1500;                       % iteration budget, wave 1 (warm starts)
  gateKm = 100;                        % GREEN means: flown control arrives
  gateMs = 10;                         %   within gateKm km and gateMs m/s

%% Initial Orbit (DRO) over one period:
             tauND0 = 1.0;
          [~,rvND0] = pumpkynPie.cr3bp.getDRO(tauND0);
              rvND0 = pumpkyn.cr3bp.cont_np(rvND0,tauND0,muStar,1e-12);
     [tauNDv,rvNDv] = pumpkyn.cr3bp.prop(tauND0,rvND0,muStar);

%% Final Orbit (Tulip) over one period:
                  Np = 7;
              tauNDf = 5*2*pi/6;
           [~,rvNDF] = pumpkyn.cr3bp.getTulip(tauNDf,Np,-1);
               rvNDF = pumpkyn.cr3bp.cont_np(rvNDF,tauNDf,muStar,1e-12);
        [tauND,rvND] = pumpkyn.cr3bp.prop(tauNDf,rvNDF,muStar);

% a state anywhere on either orbit, by phase fraction f in [0,1):
depState = @(f) interp1(tauNDv, rvNDv, mod(f,1)*tauNDv(end), 'spline');
arrState = @(f) interp1(tauND,  rvND,  mod(f,1)*tauND(end),  'spline');

%% Thruster:
 Tmax = 0.07;                                    % Max thrust (N)
   m0 = 150;                                     % Initial mass (kg)
   g0 = 9.80665 * tStar^2 / (1000 * lStar);      % Gravity in ND units
 Tmax = (Tmax / m0) * tStar^2 / (lStar * 1000);  % ND thrust accel
  Isp = 900 / tStar;                             % ND specific impulse
    c = Isp * g0;                                % ND exhaust velocity

%% The Torus Grid:
sD = (0:nD-1)/nD;                      % departure phases (wraps at 1)
sA = (0:nA-1)/nA;                      % arrival phases   (wraps at 1)
TF    = nan(nD,nA);                    % minimum time found (ND)
PASS  = false(nD,nA);                  % green?
GLOBKM= nan(nD,nA);                    % flown-control arrival miss (km)
TRIES = zeros(nD,nA);
seeds = cell(nD,nA);                   % solved trajectories, for wave 1

%% WAVE 0 -- knock on every door:
fprintf('\nWAVE 0: %d points, straight-line seed, %d-iteration budget\n', nD*nA, wave0Its);
for kD = 1:nD
    for kA = 1:nA
        [TF,PASS,GLOBKM,TRIES,seeds] = solveOnePoint(kD, kA, [], wave0Its, ...
            sD, sA, depState, arrState, Tmax, c, muStar, N, floorKm, ...
            gateKm, gateMs, lStar, tStar, TF, PASS, GLOBKM, TRIES, seeds);
    end
end
fprintf('WAVE 0 done: %d of %d green\n', nnz(PASS), nD*nA);

%% WAVE 1 -- spread from the doors that opened:
fprintf('WAVE 1: warm-start flood from the green points\n');
[iw, jw] = find(PASS);
queue = [iw, jw];
while ~isempty(queue)
    kD = queue(1,1);  kA = queue(1,2);  queue(1,:) = [];
    for step = [1 0; -1 0; 0 1; 0 -1].'
        jD = mod(kD-1+step(1), nD)+1;                  % torus wraparound
        jA = mod(kA-1+step(2), nA)+1;
        if ~PASS(jD,jA) && TRIES(jD,jA) < 4
            [TF,PASS,GLOBKM,TRIES,seeds] = solveOnePoint(jD, jA, seeds{kD,kA}, ...
                wave1Its, sD, sA, depState, arrState, Tmax, c, muStar, N, ...
                floorKm, gateKm, gateMs, lStar, tStar, TF, PASS, GLOBKM, TRIES, seeds);
            if PASS(jD,jA), queue(end+1,:) = [jD jA]; end %#ok<AGROW>
        end
    end
end
fprintf('WAVE 1 done: %d of %d green\n', nnz(PASS), nD*nA);

%% Save and Show the Torus:
rd = fullfile(here,'results');  if ~isfolder(rd), mkdir(rd); end
resFile = fullfile(rd, sprintf('phase_sweep_ps_%dx%d.mat', nD, nA));
save(resFile, 'sD','sA','TF','PASS','GLOBKM','TRIES');
plot_phase_torus(resFile, fullfile(rd, sprintf('phase_sweep_ps_%dx%d', nD, nA)));
fprintf('\nDone. %d/%d green. Open the _torus.fig to rotate the torus.\n', nnz(PASS), nD*nA);
end

%% ------------------------------------------------------------------------
function [TF,PASS,GLOBKM,TRIES,seeds] = solveOnePoint(kD, kA, seedIn, maxIts, ...
    sD, sA, depState, arrState, Tmax, c, muStar, N, floorKm, ...
    gateKm, gateMs, lStar, tStar, TF, PASS, GLOBKM, TRIES, seeds)
% SOLVEONEPOINT  One torus point: solve, then VERIFY BY FLYING.
%
%  1. Endpoints come from the two orbits at this point's phases.
%  2. The minimum-time problem is solved by direct collocation
%     (casadi_mintime_dro: Hermite-Simpson + Sundman mesh + altitude floor),
%     from the given warm start, or from its internal straight-line seed.
%  3. The gate: the solved control is FLOWN end-to-end by a high-accuracy
%     integrator (inside certify_dro_mintime). Green requires the flown
%     spacecraft to arrive within gateKm km and gateMs m/s. An optimizer
%     success certificate alone never makes a point green.
%
% INPUTS:  kD,kA grid indices; seedIn ([] = cold) with .X .U .tf; maxIts;
%          the grids, orbit-state handles, physics constants, and the running
%          result arrays.
% OUTPUTS: the updated result arrays (best result per point is kept).
rv0 = depState(sD(kD));  rvf = arrState(sA(kA));
try
    if isempty(seedIn), X0 = []; U0 = []; tf0 = [];
    else,               X0 = seedIn.X; U0 = seedIn.U; tf0 = seedIn.tf; end
    sol = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tmax, c, muStar, N, X0, U0, tf0, ...
            struct('maxIter',maxIts, 'scheme','hermite-simpson', ...
                   'sundman',true, 'minAltKm',floorKm, 'maxCpuSec',300));
    if sol.success && sol.maxDefect < 1e-9      % never fly garbage: verification
        chk = certify_dro_mintime(sol, ...       % integrations are unbounded on
            struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...  % junk iterates
            Tmax, c, struct('tfRef',[], 'verbose',false, 'posTolKm',inf));
        flownKm = chk.globKm;
        flownMs = chk.gates(strcmp({chk.gates.id},'G1bv')).value;
    else
        flownKm = Inf;  flownMs = Inf;
    end
    green = sol.success && flownKm < gateKm && flownMs < gateMs && sol.maxDefect < 1e-9;
    if green && (~PASS(kD,kA) || sol.tf < TF(kD,kA))
        TF(kD,kA) = sol.tf;  PASS(kD,kA) = true;  GLOBKM(kD,kA) = flownKm;
        tu = linspace(0, sol.tf, N+1);             % store for the neighbours
        sd.X = interp1(sol.tNodes, sol.X.', tu, 'spline').';
        sd.U = interp1(sol.tNodes, sol.U.', tu, 'spline').';
        sd.U(1:3,:) = sd.U(1:3,:) ./ max(vecnorm(sd.U(1:3,:),2,1), eps);
        sd.U(4,:)   = min(max(sd.U(4,:),0),1);
        sd.tf = sol.tf;
        seeds{kD,kA} = sd;
    elseif isnan(TF(kD,kA))
        TF(kD,kA) = sol.tf;  GLOBKM(kD,kA) = flownKm;
    end
    fprintf('  (%2d,%2d) t_f = %8.4f   flown miss = %10.2f km   %s\n', ...
        kD, kA, sol.tf, flownKm, ternary(green,'GREEN','red'));
catch ME
    fprintf('  (%2d,%2d) solver error: %s\n', kD, kA, ME.message);
end
TRIES(kD,kA) = TRIES(kD,kA) + 1;
end

%% ------------------------------------------------------------------------
function out = ternary(tf, a, b)
% TERNARY  a if tf else b.  INPUTS: tf logical; a; b.  OUTPUTS: out
if tf, out = a; else, out = b; end
end
