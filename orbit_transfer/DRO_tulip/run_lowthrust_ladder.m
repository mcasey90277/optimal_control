function run_lowthrust_ladder(stage, maxCellsIn, batchSecIn)
%% Purpose:
%
%   FRONT DOOR for the low-thrust sheet of the DRO -> tulip costate library:
%   climb from the solved 0.07 N library UP to 0.5 N, meeting the
%   high-thrust sheet (15 -> 0.5 N) in the middle. Everything a person needs
%   to change lives in the ADJUSTABLE PARAMETERS block; run this file from
%   MATLAB and the sheet builds, saving after every rung.
%
%   WHY THIS DIRECTION: coverage of the downward ladder collapses where
%   transfers reach one DRO revolution (63 pairs at 1 N -> 20 at 0.5 N) --
%   each revolution gained is a solution-family change that defeats
%   continuation from above. But 0.07 N is ALREADY SOLVED (126 pairs), so
%   the cheap route to the 0.1-0.5 N band starts there and climbs, using
%   direct-solve warm starts, which deform benignly through family changes
%   where costate continuation diverges.
%
%   SUCCESS CRITERIA -- fixed before any run; the pilot evaluates them
%   automatically and prints a verdict:
%
%     S1  re-anchor: 3/3 pilot cells re-converge at 0.07 N / target Isp
%         through all three gates
%     S2  climb: >= 2/3 pilot cells reach 0.5 N with EVERY rung passing all
%         three gates (flown < gateKm, ms residual, tfMin |dz| < accTol)
%     S3  mesh: at least one mesh mode converges through the 1-2.5
%         revolution band; if both, the faster is recorded as the default
%     S4  join: where the high-thrust sheet holds a 0.5 N entry,
%         |tf_up - tf_down| < 1e-5 ND means one continuous family;
%         disagreement is recorded as a family split (finding, not failure)
%     S5  cost: mean <= 6 min per cell, so 126 cells fit in <= ~13 h
%
%   The pilot FAILS if fewer than 2 cells climb out of 0.07 N in BOTH mesh
%   modes; the fallback is then cold flood-solving per rung (the v1 method).
%
%% Inputs:
%
%  stage                    char                    'pilot'  3 cells, BOTH
%                                                            mesh modes,
%                                                            criteria verdict
%                                                   'all'    full campaign
%                                                            (auto mesh)
%                                                   'report' census of the
%                                                            sheet so far
%                                                   (default 'pilot')
%
%  maxCellsIn               double                  Optional batch-driver
%                                                   override: cells per call
%
%  batchSecIn               double                  Optional batch-driver
%                                                   override: clean-exit
%                                                   budget, s
%
%% Outputs:
%
%   none (files under direct/results/):
%     lowthrust_pilot_uniform.mat / _sundman.mat    pilot sheets
%     lowthrust_ladder_12x12.mat                    the full sheet
%
%   Unattended runs:  ./run_ladder_batched.sh lowthrust
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('stage','var'), stage = 'pilot'; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Thrust rungs to climb (N, ASCENDING -- the 0.07 N re-anchor is automatic):
        rungs = [0.1 0.15 0.2 0.3 0.4 0.5];

%% Propulsion (must match the high-thrust sheet to allow the join test):
         ispS = 1710;                 % specific impulse (s)
         m0kg = 150;                  % initial spacecraft mass (kg)

%% Mesh strategy ('auto' scales N, K and Sundman with transfer length --
%  see mesh_policy.m, the one place the rules live):
     meshMode = 'auto';

%% Acceptance gates (identical to both existing sheets):
       gateKm = 100;                  % flown control must arrive within (km)
       accTol = 1e-6;                 % tfMin must accept within ||dz||

%% Budgets:
      msWallS = 180;                  % ms_tfmin budget per attempt (s)
      maxIter = 3000;                 % NLP iteration cap
       cpuSec = 400;                  % IPOPT wall cap per solve (s)
       maxAtt = 2;                    % attempts per cell before retiring

%% Pilot cells (well-separated, all solved in v1):
   pilotCells = [2 5; 1 3; 3 9];

%% Run control (the batched driver overrides these):
     maxCells = inf;
     batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end

%% ======================= END ADJUSTABLE PARAMETERS ======================

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'direct'), fullfile(here,'direct','lib'), ...
        fullfile(here,'direct','certify'), fullfile(here,'indirect'), ...
        fullfile(getenv('HOME'),'casadi-3.7.0'));
resDir = fullfile(here, 'direct', 'results');
base = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, 'gateKm',gateKm, ...
    'accTol',accTol, 'msWallS',msWallS, 'maxIter',maxIter, 'cpuSec',cpuSec, ...
    'maxCells',maxCells, 'batchSec',batchSec, 'maxAtt',maxAtt);

switch lower(stage)

case 'pilot'
    %% Both mesh modes on the pilot cells, then the criteria verdict:
    tPilot = tic;
    for mode = {'uniform','sundman'}
        o = base;
        o.meshMode = mode{1};
        o.cells    = pilotCells;
        o.maxAtt   = 1;
        o.logFile  = fullfile(resDir, sprintf('lowthrust_pilot_%s.txt', mode{1}));
        fprintf('\n=== PILOT, %s mesh ===\n', mode{1});
        lowthrust_ladder(fullfile(resDir, ...
            sprintf('lowthrust_pilot_%s.mat', mode{1})), o);
    end
    evaluate_pilot(resDir, pilotCells, toc(tPilot));

case 'all'
    o = base;
    o.meshMode = meshMode;
    o.logFile  = fullfile(resDir, 'lowthrust_progress.txt');
    lowthrust_ladder(fullfile(resDir, 'lowthrust_ladder_12x12.mat'), o);

case 'report'
    f = fullfile(resDir, 'lowthrust_ladder_12x12.mat');
    assert(isfile(f), 'no sheet yet -- run the pilot, then ''all''');
    Q = load(f);
    fprintf('\nLOW-THRUST SHEET: %d verified entries\n', nnz(Q.OK));
    for kr = 1:numel(Q.rungs)
        fprintf('  %5.2f N: %3d of 144 phase pairs\n', ...
                Q.rungs(kr), nnz(Q.OK(:,:,kr)));
    end
    j = Q.JOINDTF(~isnan(Q.JOINDTF));
    if ~isempty(j)
        fprintf('joins tested: %d;  same family (|dtf|<1e-5): %d\n', ...
                numel(j), nnz(j < 1e-5));
    end

otherwise
    error('unknown stage ''%s'' (pilot | all | report)', stage);
end
end

% ------------------------------------------------------------------------
function evaluate_pilot(resDir, cells, wallSec)
% EVALUATE_PILOT  Score the pilot against success criteria S1-S5.
% INPUTS: resDir results folder; cells [k x 2]; wallSec total pilot wall.
% OUTPUTS: none (printed verdict).
U = load(fullfile(resDir, 'lowthrust_pilot_uniform.mat'));
S = load(fullfile(resDir, 'lowthrust_pilot_sundman.mat'));
nR = numel(U.rungs);
nC = size(cells,1);
reanchOK = 0;  climbU = 0;  climbS = 0;
for k = 1:nC
    iD = cells(k,1);  iA = cells(k,2);
    reanchOK = reanchOK + (U.OK(iD,iA,1) || S.OK(iD,iA,1));
    climbU = climbU + all(U.OK(iD,iA,:));
    climbS = climbS + all(S.OK(iD,iA,:));
end
climbBest = max(climbU, climbS);
wU = sum(U.WALL(~isnan(U.WALL)));  wS = sum(S.WALL(~isnan(S.WALL)));
perCellMin = (wallSec/60) / (2*nC);
fprintf('\n================= PILOT VERDICT =================\n');
crit('S1 re-anchor 3/3', reanchOK == nC, sprintf('%d/%d', reanchOK, nC));
crit('S2 climb >=2/3 to top rung', climbBest >= 2, ...
     sprintf('uniform %d/%d, sundman %d/%d', climbU, nC, climbS, nC));
if climbU > 0 || climbS > 0
    if climbU >= climbS && wU <= wS, rec = 'uniform';
    elseif climbS > climbU || wS < wU, rec = 'sundman';
    else, rec = 'uniform';
    end
    crit('S3 a mesh mode works', true, sprintf( ...
        'uniform %.0fs total, sundman %.0fs -> recommend ''%s''', wU, wS, rec));
else
    crit('S3 a mesh mode works', false, 'neither mode climbed');
end
jAll = [U.JOINDTF(~isnan(U.JOINDTF)); S.JOINDTF(~isnan(S.JOINDTF))];
if isempty(jAll)
    fprintf('  S4 join: no overlap with sheet A among pilot cells (untested)\n');
else
    crit('S4 join same-family', all(jAll < 1e-5), ...
         sprintf('max |dtf| = %.2e ND over %d joins', max(jAll), numel(jAll)));
end
crit('S5 cost <= 6 min/cell', perCellMin <= 6, ...
     sprintf('%.1f min/cell -> ~%.1f h for 126 cells', ...
             perCellMin, perCellMin*126/60));
fprintf('=================================================\n');
fprintf(['If S2 failed in BOTH modes: the upward route is refuted; fall\n', ...
         'back to cold flood-solving per rung (the v1 method).\n']);
end

% ------------------------------------------------------------------------
function crit(name, pass, detail)
% CRIT  Print one criterion verdict.  INPUTS: name; pass; detail.  OUTPUTS: none.
if pass, v = 'PASS'; else, v = 'FAIL'; end
fprintf('  %-28s %s   (%s)\n', name, v, detail);
end
