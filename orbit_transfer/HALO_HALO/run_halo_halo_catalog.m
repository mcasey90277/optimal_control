function run_halo_halo_catalog(stage, maxCellsIn, batchSecIn)
%% Purpose:
%
%   FRONT DOOR for the L1 -> L2 HALO-TO-HALO costate catalog: the proven
%   three-gate pipeline (direct solve -> ms_tfmin refinement -> pumpkyn
%   tfMin acceptance), walked down a thrust ladder, swept over (L1 period x
%   L2 period) sheets. First catalog on the SCHEMA-V2 arrival-period axis
%   (sheets keyed by tau_dep x tau_arr; petal count does not apply).
%
%   Departure box: only TWO L1 southern halos are admissible under the
%   standing criteria (>= 500 km periselene, <= 100 Mm): tau = 1.8037
%   (8.00 d, periselene 2555 km) and tau = 2.7433 (12.16 d). Arrival: the
%   L2 southern set proven by the HALO_tulip campaign. Feasibility probes
%   (probe_l1_l2_halo): near pair 7/7 rungs, far pair 6/7.
%
%   Endpoints are built through costate_common/get_family_orbit; the
%   engines are shared, unmodified. The L2 -> L1 direction is a separate
%   stage ('allB') because no symmetry maps one direction onto the other.
%
%% Inputs:
%
%  stage                    char                    'pilot'  2 sheets x few
%                                                            cells (default)
%                                                   'all'    L1 -> L2 sheets
%                                                   'allB'   L2 -> L1 sheets
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
%   none (files under HALO_HALO/direct/results/catalog/):
%     hh_d<t1>_a<t2>.mat + _progress.txt      one per L1->L2 sheet
%     hhB_d<t2>_a<t1>.mat + _progress.txt     one per L2->L1 sheet
%
%   Unattended:  ./run_halo_halo_batched.sh
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('stage','var'), stage = 'pilot'; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% L1 southern halos (the complete admissible set, measured):
   l1TauList = [1.8037 2.7433];
%% L2 southern halos (the HALO_tulip campaign's proven set):
   l2TauList = [1.75 2.2 2.8 3.4];
       haloPm = -1;                 % southern; +1 is the exact z-mirror

%% Phasing grid per sheet (coarse; densify later):
           nD = 6;
           nA = 6;

%% Thrust ladder (N, high to low -- the proven set):
        rungs = [15 12 10 7 5 3 2 1.5 1];
         ispS = 1710;               % specific impulse (s)
         m0kg = 150;                % initial mass (kg)

%% Solver / gates (identical to the tulip catalogs):
            N = 400;                % collocation intervals
      floorKm = 500;                % lunar altitude floor (km)
      maxIter = 3000;
       gateKm = 100;                % flown-arrival gate (km)
       accTol = 1e-6;               % tfMin acceptance

%% Pilot: one near pair, one far pair, few cells each:
  pilotSheets = [1.8037 2.2;        % the 7/7-rung probe pair
                 2.7433 3.4];       % the far/circular pair
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
        fullfile(getenv('HOME'),'casadi-3.7.0'));
catDir = fullfile(here, 'direct', 'results', 'catalog');
if ~isfolder(catDir), mkdir(catDir); end

t2s = @(t) strrep(sprintf('%.4f',t), '.', 'p');
nameAB = @(td,ta) fullfile(catDir, sprintf('hh_d%s_a%s',  t2s(td), t2s(ta)));
nameBA = @(td,ta) fullfile(catDir, sprintf('hhB_d%s_a%s', t2s(td), t2s(ta)));
mkOpts = @(tauD,LptD,tauA,LptA,base,mc,bs) struct( ...
    'rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
    'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'gateKm',gateKm, ...
    'nD',nD, 'nA',nA, 'sD0',0, 'sA0',0, ...
    'depFamily','halo', 'depParams',struct('tau',tauD,'Lpt',LptD,'pm',haloPm), ...
    'arrFamily','halo', 'arrParams',struct('tau',tauA,'Lpt',LptA,'pm',haloPm), ...
    'tauDRO',tauD, ...               % legacy meta field; departure period
    'accTol',accTol, 'maxCells',mc, 'batchSec',bs, ...
    'logFile',[base '_progress.txt']);

switch lower(stage)

case 'pilot'
    for ks = 1:size(pilotSheets,1)
        td = pilotSheets(ks,1);  ta = pilotSheets(ks,2);
        fprintf('\n=== HH PILOT sheet L1 tau=%.4f -> L2 tau=%.2f ===\n', td, ta);
        base = nameAB(td, ta);
        thrust_ladder_library([base '.mat'], ...
            mkOpts(td, 1, ta, 2, base, pilotCellsMax, inf));
    end
    fprintf(['\nPilot done. SUCCESS = both sheets hold three-gate entries\n', ...
             'at 15 N and at least one lower rung; then stage ''all''.\n']);

case {'all', 'allb'}
    ba = strcmpi(stage, 'allb');
    tAll = tic;
    allDone = true;
    for td = l1TauList
        for ta = l2TauList
            if toc(tAll) > batchSec, allDone = false; break, end
            if ba
                base = nameBA(ta, td);
                opts = mkOpts(ta, 2, td, 1, base, maxCells, ...
                              max(0, batchSec - toc(tAll)));
            else
                base = nameAB(td, ta);
                opts = mkOpts(td, 1, ta, 2, base, maxCells, ...
                              max(0, batchSec - toc(tAll)));
            end
            P = thrust_ladder_library([base '.mat'], opts);
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
        if ba, fprintf(fid, 'CATALOG B (L2->L1) ALL SHEETS COMPLETE\n');
        else,  fprintf(fid, 'CATALOG ALL SHEETS COMPLETE\n');
        end
        fclose(fid);
    end

case 'report'
    for pre = {'hh_', 'hhB_'}
        F = dir(fullfile(catDir, [pre{1} 'd*.mat']));
        if isempty(F), continue, end
        fprintf('\n%s sheets:\n', pre{1});
        for k = 1:numel(F)
            Q = load(fullfile(F(k).folder, F(k).name));
            fprintf('  %-26s %4d entries, %2d/%d pairs, %2d full ladders\n', ...
                erase(F(k).name,'.mat'), nnz(Q.OK), nnz(any(Q.OK,3)), ...
                nD*nA, nnz(sum(Q.OK,3)==numel(Q.rungs)));
        end
    end

otherwise
    error('unknown stage ''%s'' (pilot | all | allB | report)', stage);
end
end
