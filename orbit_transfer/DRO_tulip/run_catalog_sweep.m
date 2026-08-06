function run_catalog_sweep(stage, maxCellsIn, batchSecIn)
%% Purpose:
%
%   FRONT DOOR for the MULTI-ORBIT costate catalog: the thrust-ladder
%   pipeline swept over the admissible box of orbit parameters -- every
%   (DRO period) x (tulip petal count) pair gets its own phasing-torus
%   SHEET, built by the proven three-gate pipeline (direct solve ->
%   ms_tfmin -> pumpkyn tfMin acceptance) walked down the thrust rungs.
%
%   The admissible box comes from survey_period_bounds (Darin's criteria:
%   >= 500 km lunar altitude, within 100 Mm of the Moon):
%       DRO   tau (= period, ND) in [0.05, 3.25]
%       tulip Np in 3..14 -- and pm = +1 is the exact z-mirror of pm = -1
%       (CR3BP symmetry), so only ONE branch is ever computed.
%
%   COARSE FIRST: the default grid below is deliberately coarse (Darin:
%   assemble a coarse dataset, interpolate, refine where customers need
%   it). Densify later by adding values to the lists -- finished sheets
%   are skipped automatically.
%
%   RUN THE PILOT BEFORE 'all': two sheets away from the flown (tau = 1,
%   Np = 7) point, a few cells each. It exists to catch anything silently
%   specific to the original orbit pair before a day of compute is
%   committed. Success = both pilot sheets produce three-gate entries at
%   the anchor thrust and at least one walked rung.
%
%% Inputs:
%
%  stage                    char                    'pilot'  2 sheets x few
%                                                            cells (default)
%                                                   'all'    every sheet in
%                                                            the grid
%                                                   'report' catalog census
%
%  maxCellsIn               double                  Optional batch override:
%                                                   cells per call
%
%  batchSecIn               double                  Optional batch override:
%                                                   clean-exit budget (s)
%
%% Outputs:
%
%   none (files under direct/results/catalog/):
%     ladder_tau<t>_Np<n>.mat     one sheet per orbit pair
%     <same>_progress.txt         per-sheet log
%
%   Unattended:  ./run_ladder_batched.sh catalog
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('stage','var'), stage = 'pilot'; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% The coarse orbit grid (admissible box: tau in [0.05,3.25], Np in 3..14):
   tauDROList = [0.5 1.0 2.0 3.0];   % DRO periods (ND; tau IS the period)
       NpList = [5 7 9 12];          % tulip petal counts
           pm = -1;                  % ONE branch only (pm=+1 is its mirror)

%% Phasing resolution per sheet (coarse; densify later):
           nD = 6;                   % departure phases
           nA = 6;                   % arrival phases

%% Thrust ladder (N, high to low -- the proven set; sub-1 N deferred):
        rungs = [15 12 10 7 5 3 2 1.5 1];
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (identical to the flown tau=1, Np=7 campaign):
            N = 400;                 % collocation intervals
      floorKm = 500;                 % lunar altitude floor (km)
      maxIter = 3000;
       gateKm = 100;                 % flown-arrival gate (km)
       accTol = 1e-6;                % tfMin acceptance

%% Pilot: two sheets chosen to move each orbit axis separately:
   pilotSheets = [2.0 7;             % new DRO period, known tulip
                  1.0 9];            % known DRO, new petal count
 pilotCellsMax = 6;                  % cells per pilot sheet

%% Run control (batched driver overrides):
     maxCells = inf;
     batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end

%% ======================= END ADJUSTABLE PARAMETERS ======================

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'direct'), fullfile(here,'direct','lib'), ...
        fullfile(here,'direct','certify'), fullfile(here,'indirect'), ...
        fullfile(getenv('HOME'),'casadi-3.7.0'));
catDir = fullfile(here, 'direct', 'results', 'catalog');
if ~isfolder(catDir), mkdir(catDir); end

sheetName = @(tau,Np) fullfile(catDir, sprintf('ladder_tau%s_Np%d', ...
    strrep(sprintf('%.2f',tau),'.','p'), Np));
mkOpts = @(tau,Np,mc,bs) struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
    'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'gateKm',gateKm, ...
    'nD',nD, 'nA',nA, 'sD0',0, 'sA0',0, ...
    'tauDRO',tau, 'tauTulip',(Np-2)/(Np-1)*2*pi, 'NpTulip',Np, 'pmTulip',pm, ...
    'accTol',accTol, 'maxCells',mc, 'batchSec',bs, ...
    'logFile',[sheetName(tau,Np) '_progress.txt']);

switch lower(stage)

case 'pilot'
    for ks = 1:size(pilotSheets,1)
        tau = pilotSheets(ks,1);  Np = pilotSheets(ks,2);
        fprintf('\n=== CATALOG PILOT sheet tau=%.2f, Np=%d ===\n', tau, Np);
        thrust_ladder_library([sheetName(tau,Np) '.mat'], ...
            mkOpts(tau, Np, pilotCellsMax, inf));
    end
    fprintf(['\nPilot done. SUCCESS = both sheets hold three-gate entries at\n', ...
             '15 N and at least one lower rung; then launch ''all'' (or\n', ...
             './run_ladder_batched.sh catalog for the unattended run).\n']);

case 'all'
    tAll = tic;
    allDone = true;
    for tau = tauDROList
        for Np = NpList
            if toc(tAll) > batchSec, allDone = false; break, end
            P = thrust_ladder_library([sheetName(tau,Np) '.mat'], ...
                mkOpts(tau, Np, maxCells, max(0, batchSec - toc(tAll))));
            % a sheet is open while any cell is neither finished nor retired
            open_ = false;
            for iD = 1:nD
                for iA = 1:nA
                    if ~P.OK(iD,iA,end) && P.ATT(iD,iA) < 2, open_ = true; end
                end
            end
            if open_, allDone = false; end
        end
        if toc(tAll) > batchSec, allDone = false; break, end
    end
    if allDone
        fid = fopen(fullfile(catDir,'combined_progress.txt'),'a');
        fprintf(fid, 'CATALOG ALL SHEETS COMPLETE\n');  fclose(fid);
    end

case 'report'
    fprintf('\nCATALOG CENSUS (%d x %d sheets, %d rungs, %dx%d phasing):\n', ...
        numel(tauDROList), numel(NpList), numel(rungs), nD, nA);
    for tau = tauDROList
        for Np = NpList
            f = [sheetName(tau,Np) '.mat'];
            if isfile(f)
                Q = load(f);
                fprintf('  tau=%.2f Np=%2d: %4d entries, %2d/%d pairs, %2d full ladders\n', ...
                    tau, Np, nnz(Q.OK), nnz(any(Q.OK,3)), nD*nA, ...
                    nnz(sum(Q.OK,3)==numel(Q.rungs)));
            else
                fprintf('  tau=%.2f Np=%2d: (not started)\n', tau, Np);
            end
        end
    end

otherwise
    error('unknown stage ''%s'' (pilot | all | report)', stage);
end
end
