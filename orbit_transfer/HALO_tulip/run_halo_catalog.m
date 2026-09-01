function run_halo_catalog(stage, maxCellsIn, batchSecIn)
%% Purpose:
%
%   FRONT DOOR for the HALO -> TULIP costate catalog: the proven three-gate
%   pipeline (direct solve -> ms_tfmin refinement -> pumpkyn tfMin
%   acceptance), walked down a thrust ladder, swept over (halo period x
%   tulip petal count) sheets. Edit the ADJUSTABLE PARAMETERS block and run
%   this file from MATLAB; sheets save after every rung and every run is
%   resumable.
%
%   The halo box comes from the measured admissibility survey
%   (HALO_tulip/direct/results/halo_bounds.mat, Darin's criteria: >= 500 km
%   periselene, within 100 Mm of the Moon):
%     - L2 SOUTHERN is continuously admissible for periods 7.1-15.1 days
%       (tau 1.61-3.42); the northern branch is its exact z-mirror.
%     - L1 offers only two clean members; L3 lies outside the box entirely.
%
%   Endpoints are built through costate_common/get_family_orbit, so this
%   front door differs from the DRO catalog's ONLY in the departure-family
%   parameters -- the engines are shared, unmodified.
%
%% Inputs:
%
%  stage                    char                    'pilot'  2 sheets x few
%                                                            cells (default)
%                                                   'all'    every sheet
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
%   none (files under HALO_tulip/direct/results/catalog/):
%     halo_tau<t>_Np<n>.mat + _progress.txt   one per sheet
%
%   Unattended:  ../DRO_tulip/run_ladder_batched.sh is the pattern; this
%   campaign uses ./run_halo_batched.sh
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('stage','var'), stage = 'pilot'; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Departure: L2 SOUTHERN halos, selected by period (measured box 1.61-3.42):
  haloTauList = [1.75 2.2 2.8 3.4];  % periods (ND; 1 ND = 4.43 d)
      haloLpt = 2;                   % Lagrange point (L2 = the workhorse)
       haloPm = -1;                  % southern branch (+1 = its z-mirror)

%% Arrival: tulip petal counts (periods locked by resonance):
       NpList = [5 7 9 12];
      tulipPm = -1;

%% Phasing grid per sheet (coarse; densify later):
           nD = 6;
           nA = 6;

%% Thrust ladder (N, high to low -- the proven set):
        rungs = [15 12 10 7 5 3 2 1.5 1];
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (identical to the DRO catalog):
            N = 400;                 % collocation intervals
      floorKm = 500;                 % lunar altitude floor (km)
      maxIter = 3000;
       gateKm = 100;                 % flown-arrival gate (km)
       accTol = 1e-6;                % tfMin acceptance

%% Pilot: move each orbit axis separately, few cells each:
  pilotSheets = [2.2 7;              % mid-box halo, known tulip
                 2.8 9];             % different halo period AND petal count
pilotCellsMax = 6;

%% Run control (the batched driver overrides):
     maxCells = inf;
     batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end

%% ======================= END ADJUSTABLE PARAMETERS ======================

here = fileparts(mfilename('fullpath'));
droDir = fullfile(fileparts(here), 'DRO_tulip');
addpath(fullfile(droDir,'direct'), fullfile(droDir,'direct','lib'), ...
        fullfile(droDir,'direct','certify'), fullfile(droDir,'indirect'), ...
        fullfile(fileparts(here),'costate_common'), ...
        fullfile(fileparts(fileparts(here)),'oclib'), ...   %oc.local_residual (path-gap fix 2026-08-31)
        fullfile(getenv('HOME'),'casadi-3.7.0'));
catDir = fullfile(here, 'direct', 'results', 'catalog');
if ~isfolder(catDir), mkdir(catDir); end

sheetName = @(tau,Np) fullfile(catDir, sprintf('halo_tau%s_Np%d', ...
    strrep(sprintf('%.2f',tau),'.','p'), Np));
mkOpts = @(tau,Np,mc,bs) struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
    'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'gateKm',gateKm, ...
    'nD',nD, 'nA',nA, 'sD0',0, 'sA0',0, ...
    'depFamily','halo', 'depParams',struct('tau',tau,'Lpt',haloLpt,'pm',haloPm), ...
    'arrFamily','tulip', 'arrParams',struct('Np',Np,'pm',tulipPm), ...
    'NpTulip',Np, 'tauTulip',2*pi*(Np-2)/(Np-1), 'pmTulip',tulipPm, ...
    'tauDRO',tau, ...                 % legacy meta field; records the period
    'accTol',accTol, 'maxCells',mc, 'batchSec',bs, ...
    'logFile',[sheetName(tau,Np) '_progress.txt']);

switch lower(stage)

case 'pilot'
    for ks = 1:size(pilotSheets,1)
        tau = pilotSheets(ks,1);  Np = pilotSheets(ks,2);
        fprintf('\n=== HALO PILOT sheet tau=%.2f (L%d, pm=%d), Np=%d ===\n', ...
                tau, haloLpt, haloPm, Np);
        thrust_ladder_library([sheetName(tau,Np) '.mat'], ...
            mkOpts(tau, Np, pilotCellsMax, inf));
    end
    fprintf(['\nPilot done. SUCCESS = both sheets hold three-gate entries\n', ...
             'at 15 N and at least one lower rung; then stage ''all''.\n']);

case 'all'
    tAll = tic;
    allDone = true;
    for tau = haloTauList
        for Np = NpList
            if toc(tAll) > batchSec, allDone = false; break, end
            P = thrust_ladder_library([sheetName(tau,Np) '.mat'], ...
                mkOpts(tau, Np, maxCells, max(0, batchSec - toc(tAll))));
            open_ = false;
            for iD = 1:nD
                for iA = 1:nA
                    if ~any(P.OK(iD,iA,:)) && P.ATT(iD,iA) < 2, open_ = true; end
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
    fprintf('\nHALO->TULIP CENSUS (%d x %d sheets):\n', ...
            numel(haloTauList), numel(NpList));
    for tau = haloTauList
        for Np = NpList
            f = [sheetName(tau,Np) '.mat'];
            if isfile(f)
                Q = load(f);
                fprintf('  halo tau=%.2f Np=%2d: %4d entries, %2d/%d pairs, %2d full ladders\n', ...
                    tau, Np, nnz(Q.OK), nnz(any(Q.OK,3)), nD*nA, ...
                    nnz(sum(Q.OK,3)==numel(Q.rungs)));
            else
                fprintf('  halo tau=%.2f Np=%2d: (not started)\n', tau, Np);
            end
        end
    end

otherwise
    error('unknown stage ''%s'' (pilot | all | report)', stage);
end
end
