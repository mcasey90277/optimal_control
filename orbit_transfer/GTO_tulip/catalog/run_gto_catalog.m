function run_gto_catalog(sheetSel, maxCellsIn, batchSecIn)
%% Purpose:
%
%   FRONT DOOR for the GTO -> TULIP costate catalog: the proven three-gate
%   pipeline (direct solve -> ms_tfmin refinement -> pumpkyn tfMin
%   acceptance), walked down a thrust ladder, swept over (GTO departure
%   orientation x tulip petal count) sheets. THE SWATH DRIVER:
%   run_gto_catalog() solves all 16 sheets; run_gto_catalog(sheetSel)
%   solves the subset indexed by sheetSel. Resume/regeneration of unsolved
%   cells is automatic (thrust_ladder_library revisits no-OK cells on
%   every call; a batched shell driver is how an unattended run survives a
%   hung solver call -- see run_gto_batched.sh).
%
%   Endpoints are built through costate_common/get_family_orbit, so this
%   front door differs from the DRO/HALO catalogs' ONLY in the departure-
%   family parameters ('gto', keyed by orientDeg) -- the engine is shared,
%   unmodified (clone of HALO_tulip/run_halo_catalog.m's 'all' stage).
%
%   SHEET KEY RULE: opts.tauDRO = orientDeg for every sheet, NOT a Kepler
%   period -- the picker (and every consumer downstream) selects sheets by
%   nearest (tauDRO, Np); the GTO family's own Kepler period rides in
%   depParams.orientDeg / get_family_orbit's info, never in tauDRO.
%
%   DEPARTURE-PHASE GRID CAVEAT: thrust_ladder_library only accepts a
%   UNIFORM departure-phase grid (opts.sD0 + opts.nD -> mod(sD0+(0:nD-1)/
%   nD,1)); it has no explicit-fraction-vector input. A review requested a
%   grid uniform in TRUE ANOMALY (denser near perigee on this e~0.72
%   ellipse) mapped through Kepler's equation to non-uniform time
%   fractions -- that would require an engine change, which this task does
%   NOT make (house rule: house engines are not modified in place to serve
%   one campaign). sD0=0, nD=12 UNIFORM IN TIME FRACTION is used instead;
%   this still lands cleanly on perigee (index 1) and apogee (index 7 of
%   12, since 6/12=0.5) so both regimes are sampled, just not with extra
%   density at perigee. See task-3-report.md.
%
%% Inputs:
%
%  sheetSel                 [1 x k]                 Index vector into the
%                                                   16-sheet list (default
%                                                   1:16, i.e. all sheets)
%
%  maxCellsIn               double                  Optional batch
%                                                   override: fresh cells
%                                                   per sheet this call
%
%  batchSecIn                double                 Optional batch
%                                                   override: clean-exit
%                                                   wall budget (s), shared
%                                                   across the selected
%                                                   sheets this call
%
%% Outputs:
%
%   none (files under results/catalog/, relative to this file's folder):
%     gto_o<ORIENT3>_Np<NP>.mat + _progress.txt   one pair per sheet
%     combined_progress.txt                       driver-level log +
%                                                  completion marker
%
%   Unattended: ./run_gto_batched.sh is the batched driver (OS-kill
%   backstop for a hung solver call, per the matlab-campaign pattern).
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('sheetSel','var') || isempty(sheetSel), sheetSel = 1:16; end
maxCells = inf;  batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Departure: GTO pseudo-family, selected by orientation (sheet axis 1):
      orients = [0 90 180 270];               % deg

%% Arrival: tulip petal counts (sheet axis 2; periods locked by resonance):
       Nps     = [5 7 9 12];
       tulipPm = -1;

%% Phasing grid per sheet (departure x arrival; see the module-header
%% caveat above -- uniform-only, sD0=0/sA0=0 origins):
           nD = 12;
           nA = 6;
           sD0 = 0;
           sA0 = 0;

%% Thrust ladder (N, high to low -- the proven set):
        rungs = [15 12 10 7 5 3 2 1.5 1];
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (identical to the DRO/HALO catalogs):
            N = 400;                 % collocation intervals
      floorKm = 500;                 % lunar altitude floor (km)
      maxIter = 3000;
       gateKm = 100;                 % flown-arrival gate (km)
       accTol = 1e-6;                % tfMin acceptance

%% ======================= END ADJUSTABLE PARAMETERS ======================

% Review fix (c): assert the departure-phase grid the engine will build
% internally has no duplicate phase (no sD=1 aliasing sD=0):
sDCheck = mod(sD0 + (0:nD-1)/nD, 1);
assert(numel(unique(sDCheck)) == nD, 'run_gto_catalog:sDgrid', ...
    'departure-phase grid has a duplicate phase');

setup_paths();
here = fileparts(mfilename('fullpath'));
catDir = fullfile(here, 'results', 'catalog');
if ~isfolder(catDir), mkdir(catDir); end

sheets = {};
for od = orients
    for np = Nps
        sheets{end+1} = struct('orientDeg', od, 'Np', np, ...
            'tag', sprintf('gto_o%03d_Np%d', od, np));
    end
end
nSheets = numel(sheets);
assert(nSheets == 16, 'run_gto_catalog:sheetCount', ...
    'expected 16 sheets, got %d', nSheets);

sheetSel = sheetSel(sheetSel >= 1 & sheetSel <= nSheets);
assert(~isempty(sheetSel), 'run_gto_catalog:sheetSel', ...
    'sheetSel selects no valid sheet (1..%d)', nSheets);

tAll = tic;
allDone = true;
for ks = sheetSel
    if toc(tAll) > batchSec, allDone = false; break, end
    s = sheets{ks};
    outMat  = fullfile(catDir, [s.tag '.mat']);
    logFile = fullfile(catDir, [s.tag '_progress.txt']);
    fprintf('\n=== GTO sheet %d/%d: %s (orient=%d deg, Np=%d) ===\n', ...
            ks, nSheets, s.tag, s.orientDeg, s.Np);
    opts = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
        'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'gateKm',gateKm, ...
        'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
        'depFamily','gto', 'depParams',struct('orientDeg',s.orientDeg), ...
        'arrFamily','tulip', 'arrParams',struct('Np',s.Np,'pm',tulipPm), ...
        'NpTulip',s.Np, 'tauTulip',2*pi*(s.Np-2)/(s.Np-1), 'pmTulip',tulipPm, ...
        'tauDRO',s.orientDeg, ...          % SHEET KEY RULE: orientDeg
        'accTol',accTol, 'maxCells',maxCells, ...
        'batchSec',max(0, batchSec - toc(tAll)), 'logFile',logFile);
    P = thrust_ladder_library(outMat, opts);
    open_ = false;
    for iD = 1:nD
        for iA = 1:nA
            if ~any(P.OK(iD,iA,:)) && P.ATT(iD,iA) < 2, open_ = true; end
        end
    end
    if open_, allDone = false; end
end

logComb = fullfile(catDir, 'combined_progress.txt');
if allDone && numel(sheetSel) == nSheets
    fid = fopen(logComb, 'a');
    fprintf(fid, 'CATALOG ALL SHEETS COMPLETE\n');  fclose(fid);
elseif allDone
    fid = fopen(logComb, 'a');
    fprintf(fid, 'SELECTED SHEETS COMPLETE (%s)\n', mat2str(sheetSel));
    fclose(fid);
end
end
