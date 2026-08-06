function run_costate_library(stage, maxCellsIn, batchSecIn)
%% Purpose:
%
%   ONE FRONT DOOR for building a DRO -> tulip costate library: a catalog of
%   converged minimum-time transfer solutions (initial costates + flight
%   time), indexed by orbital phasing and thrust level, that pumpkyn's
%   single-shooting solver accepts unchanged.
%
%   Everything a person needs to change lives in the ADJUSTABLE PARAMETERS
%   block below. Edit those, run this file, and the library regenerates.
%
%   THE PIPELINE, in three steps (see doc/costate_library_process.tex):
%
%     1. SOLVE     A direct collocation solve produces a trajectory AND a
%                  costate trajectory. The solved control is then FLOWN
%                  end-to-end by a high-accuracy integrator: if the
%                  spacecraft does not actually arrive, the case is rejected
%                  no matter what the optimizer reported.
%
%     2. REFINE    Collocation costates are only ~1e-3 accurate, and single
%                  shooting amplifies that ~1000x over a long arc. ms_tfmin
%                  (multiple shooting, analytic Jacobian from pumpkyn's own
%                  segment STMs) drives them to ~1e-12.
%
%     3. ACCEPT    pumpkyn.cr3bp.tfMin must return the refined solution
%                  UNCHANGED. Only then does the entry enter the library.
%
%   Thrust levels are covered by a LADDER: anchor at the highest thrust,
%   where the transfer is nearly impulsive and converges from a cold start in
%   seconds, then walk thrust DOWN with each rung warm-started from the one
%   above. Continuation keeps a phase pair on one solution family; solving
%   each thrust independently does not.
%
%% Inputs:
%
%  maxCellsIn               double                  Optional: process at
%                                                   most this many cells
%                                                   this call, then return
%                                                   (used by the batched
%                                                   driver). Default inf.
%
%  batchSecIn               double                  Optional: stop cleanly
%                                                   between cells after this
%                                                   many seconds. Default
%                                                   inf.
%
%  stage                    string                  Which part to run:
%                                                   'ladder'  build the
%                                                             library
%                                                   'extend'  continue to
%                                                             lower thrust
%                                                   'densify' fill per-rung
%                                                             gaps among
%                                                             solving pairs
%                                                   'package' write the
%                                                             shareable .mat
%                                                   'all'     all three
%                                                   (default 'all')
%
%% Outputs:
%
%   none (files written under results/):
%     thrust_ladder_12x12.mat        raw ladder results (all diagnostics)
%     costate_lib_dro_tulip_v2.mat   the shareable library
%
%   For long unattended runs use the batched driver instead, which survives
%   solver hangs:  ./run_ladder_batched.sh
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('stage','var'), stage = 'all'; end

%% ===================== ADJUSTABLE PARAMETERS ============================
% Everything below this line is meant to be edited.

%% Departure orbit (DRO). tau IS THE PERIOD in nondimensional time:
%  getDRO(tau) returns the DRO whose period equals tau. tau = 1 is 4.43 days.
      tauDRO = 1.0;

%% Arrival orbit (tulip):
     tauTulip = 5*2*pi/6;             % period (ND) -- 5/6 of a lunar month
      NpTulip = 7;                    % number of lobes
      pmTulip = -1;                   % branch of the family

%% Phasing grid. Departure phase x arrival phase, both wrapping (a torus):
           nD = 12;                   % departure phases
           nA = 12;                   % arrival phases
          sD0 = 0;                    % grid origin, departure (fraction)
          sA0 = 0.075378;             % grid origin, arrival -- the tulip's
                                      % max-velocity-angle anchor. KEEP this
                                      % to stay aligned with the 0.07 N
                                      % library, whose cells then correspond
                                      % one-to-one.

%% Propulsion:
   thrustRungs = [15 12 10 7 5 3 2 1.5 1];   % N, HIGH to LOW; the ladder
                                             % walks down this list
      extraRungs = [0.9 0.8 0.7 0.6 0.5];    % used by stage 'extend'
           ispS = 1710;               % specific impulse (s)
          m0kg = 150;                 % initial spacecraft mass (kg)

%% Solver settings:
            N = 400;                  % collocation intervals per solve
      floorKm = 500;                  % minimum lunar altitude (km); without
                                      % it the optimizer will happily fly
                                      % through the Moon
      maxIter = 3000;                 % NLP iteration cap
        msSeg = [12 24];              % multiple-shooting segment ladder:
                                      % more segments for rougher seeds
      msWallS = 120;                  % ms_tfmin budget per attempt (s)

%% Acceptance gates -- a case joins the library only if it clears all three:
       gateKm = 100;                  % flown control must arrive within this
       accTol = 1e-6;                 % ||z - z_tfMin|| after acceptance
      msTolR  = 1e-10;                % multiple-shooting residual

%% Run control (for the in-process runner; the batched driver overrides):
     maxCells = inf;                  % cells to process this call
     batchSec = inf;                  % stop cleanly between cells after (s)
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end

%% ======================= END ADJUSTABLE PARAMETERS ======================

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'direct'), fullfile(here,'direct','lib'), ...
        fullfile(here,'direct','certify'), fullfile(here,'indirect'), ...
        fullfile(getenv('HOME'),'casadi-3.7.0'));
resDir = fullfile(here,'direct','results');
if ~isfolder(resDir), mkdir(resDir); end
ladderFile = fullfile(resDir, sprintf('thrust_ladder_%dx%d.mat', nD, nA));
libFile    = fullfile(resDir, 'costate_lib_dro_tulip_v2.mat');

ladderOpts = struct('rungs',thrustRungs, 'ispS',ispS, 'm0kg',m0kg, ...
    'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'gateKm',gateKm, ...
    'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
    'maxCells',maxCells, 'batchSec',batchSec, ...
    'tauDRO',tauDRO, 'tauTulip',tauTulip, 'NpTulip',NpTulip, 'pmTulip',pmTulip, ...
    'msSeg',msSeg, 'msWallS',msWallS, 'accTol',accTol, 'msTolR',msTolR);

%% Stage 1-3: build the library over the thrust ladder:
if any(strcmpi(stage, {'ladder','all'}))
    fprintf('\n=== LADDER: %d phase pairs x %d thrust rungs ===\n', ...
            nD*nA, numel(thrustRungs));
    ladderOpts.logFile = fullfile(resDir,'ladder_progress.txt');
    thrust_ladder_library(ladderFile, ladderOpts);
end

%% Densify: fill per-rung gaps among pairs that solve somewhere:
if any(strcmpi(stage, {'densify'}))
    fprintf('\n=== DENSIFY: per-rung gaps in the %g-%g N band ===\n', ...
            min(thrustRungs), max(thrustRungs));
    % rungSel = the front door's declared rung set: densify fills only the
    % band this campaign is about; deeper rungs are their own campaign
    densify_ladder(ladderFile, struct('gateKm',gateKm, 'accTol',accTol, ...
        'msWallS',msWallS, 'maxCells',maxCells, 'batchSec',batchSec, ...
        'rungSel',thrustRungs, ...
        'logFile',fullfile(resDir,'densify_progress.txt')));
end

%% Optional: continue the ladder to lower thrust:
if any(strcmpi(stage, {'extend','all'}))
    fprintf('\n=== EXTEND: adding rungs [%s] N ===\n', num2str(extraRungs));
    extend_thrust_ladder(ladderFile, extraRungs, ...
        struct('K',msSeg(end), 'wallSec',msWallS, 'gateKm',gateKm, ...
               'maxCells',maxCells, 'batchSec',batchSec, ...
               'logFile',fullfile(resDir,'extend_progress.txt')));
end

%% Package the shareable library:
if any(strcmpi(stage, {'package','all'}))
    fprintf('\n=== PACKAGE ===\n');
    build_costate_lib_v2(ladderFile, libFile);
    fprintf(['\nSend to a collaborator:\n', ...
             '  %s\n', ...
             '  indirect/costate_lib_pick.m        (lookup by phase + thrust)\n', ...
             '  indirect/costate_lib_v2_example.m  (worked example)\n'], libFile);
end

end
