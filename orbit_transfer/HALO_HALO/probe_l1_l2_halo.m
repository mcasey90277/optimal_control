function probe_l1_l2_halo()
%% Purpose:
%
%   FEASIBILITY PROBE for L1 -> L2 halo-to-halo transfers (Darin's ask,
%   2026-08-07): runs ONE phasing cell of the proven three-gate pipeline
%   (direct solve -> ms_tfmin -> pumpkyn tfMin acceptance) with a halo at
%   BOTH ends -- the ladder engine takes departure and arrival families
%   independently, so this is a parameterization, not new machinery.
%
%   Departure: the short-period admissible L1 southern halo (tau = 1.8037
%   ND = 8.00 days, periselene 2555 km -- one of only two L1 members inside
%   Darin's box). Arrival: the L2 southern halo at tau = 2.2 ND, a
%   mid-catalog member of the L2 sheet set.
%
%% Inputs:
%
%   none (edit the ADJUSTABLE PARAMETERS block)
%
%% Outputs:
%
%   none (sheet file direct/results/l1l2/probe_L1_L2.mat + progress log)
%
%% Revision History:
%  M. Casey                                                   (c) 08/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% ===================== ADJUSTABLE PARAMETERS ============================
   tauL1 = 1.8037;                 % departure L1 southern halo period (ND)
   tauL2 = 2.2;                    % arrival L2 southern halo period (ND)
   rungs = [15 10 5 3 2 1.5 1];    % thrust ladder (N)
    ispS = 1710;                   % specific impulse (s)
    m0kg = 150;                    % initial mass (kg)
       N = 400;                    % collocation intervals
 floorKm = 500;                    % lunar altitude floor (km)
  gateKm = 100;                    % flown-arrival gate (km)
  accTol = 1e-6;                   % tfMin acceptance
%% ======================= END ADJUSTABLE PARAMETERS ======================

here = fileparts(mfilename('fullpath'));
droDir = fullfile(fileparts(here), 'DRO_tulip');
addpath(fullfile(droDir,'direct'), fullfile(droDir,'direct','lib'), ...
        fullfile(droDir,'direct','certify'), fullfile(droDir,'indirect'), ...
        fullfile(fileparts(here),'costate_common'), ...
        fullfile(getenv('HOME'),'casadi-3.7.0'));
outDir = fullfile(here, 'direct', 'results', 'l1l2');
if ~isfolder(outDir), mkdir(outDir); end

opts = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
    'N',N, 'floorKm',floorKm, 'maxIter',3000, 'gateKm',gateKm, ...
    'nD',1, 'nA',1, 'sD0',0, 'sA0',0, ...
    'depFamily','halo', 'depParams',struct('tau',tauL1,'Lpt',1,'pm',-1), ...
    'arrFamily','halo', 'arrParams',struct('tau',tauL2,'Lpt',2,'pm',-1), ...
    'tauDRO',tauL1, ...              % legacy meta field; departure period
    'accTol',accTol, 'maxCells',inf, 'batchSec',inf, ...
    'logFile',fullfile(outDir,'probe_L1_L2_progress.txt'));

P = thrust_ladder_library(fullfile(outDir,'probe_L1_L2.mat'), opts);
fprintf('\nL1->L2 probe: %d of %d rungs passed all three gates\n', ...
        nnz(P.OK), numel(rungs));
end
