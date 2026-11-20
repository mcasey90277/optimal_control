function sheetMat = fixture_dro_elfo_sheet(outMat, cut)
%% Purpose:
%
%   A small REAL DRO -> ELFO sheet for the Task 7 tests, written to outMat
%   (a tempname path; nothing is written into the repo). The entries are
%   CUT from the finished matched-experiment sheet
%   results/compare/dro_elfo_tau1_raan000.mat (6 x 6 phases, 11 rungs),
%   read-only: the chosen departure rows, arrival columns and leading rungs
%   are copied with every engine array, sheet_key and meta, and meta.sD0 /
%   meta.sA0 are set to the cut's first phases. The cut must form a uniform
%   engine grid (phases sD0 + (0:n-1)/n), so the engine and the repair tool
%   accept the cut sheet as their own: e.g. rows [2 5] of 6 give
%   sD = [1/6 2/3] = 1/6 + (0:1)/2.
%
%   No source sheet (a fresh clone, or results/compare removed): the
%   fixture is SOLVED, not skipped -- through the front door
%   run_dro_elfo_torus on the cut's grid and rungs, unfenced, cold. It
%   first prints ONE line saying so, with the cell and rung count and the
%   expected time (about 3 min per cell at 2 rungs on an idle machine, e.g.
%   about 12 min for the chain test's 2 x 2 cut, 18 min for a 6 x 1 cut),
%   so a long test run is never silent.
%
%% Inputs:
%
%  outMat                   char                    Where to write the sheet
%                                                   (its name is the sheet
%                                                   tag + .mat for the
%                                                   chain test)
%
%  cut                      struct                  .iD [1 x nD] departure
%                                                   rows, .iA [1 x nA]
%                                                   arrival columns, .nR
%                                                   leading rungs (>= 2)
%                                                   kept; (.src optional
%                                                   source sheet path)
%
%% Outputs:
%
%  sheetMat                 char                    outMat
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  fix round 1: the solve fallback announces itself (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
src  = fullfile(fileparts(here), 'results', 'compare', 'dro_elfo_tau1_raan000.mat');
if isfield(cut, 'src'), src = cut.src; end
if ~isfolder(fileparts(outMat)), mkdir(fileparts(outMat)); end
sheetMat = outMat;

if isfile(src)
    Q = load(src);
    nD0 = size(Q.OK, 1);  nA0 = size(Q.OK, 2);
    iD = cut.iD(:).';  iA = cut.iA(:).';  kr = 1:cut.nR;
    assert(all(iD >= 1 & iD <= nD0) && all(iA >= 1 & iA <= nA0) && cut.nR <= numel(Q.rungs), ...
           'fixture_dro_elfo_sheet:cut', 'cut outside the source sheet');
    sD = Q.sD(iD);  sA = Q.sA(iA);
    assert(max(abs(sD - mod(sD(1) + (0:numel(iD) - 1)/numel(iD), 1))) < 1e-12 && ...
           max(abs(sA - mod(sA(1) + (0:numel(iA) - 1)/numel(iA), 1))) < 1e-12, ...
           'fixture_dro_elfo_sheet:grid', 'the cut phases are not a uniform engine grid');
    P = Q;
    for f = {'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL', 'OK'}
        P.(f{1}) = Q.(f{1})(iD, iA, kr);
    end
    P.Z8    = Q.Z8(:, iD, iA, kr);
    P.ATT   = Q.ATT(iD, iA);
    P.rungs = Q.rungs(kr);
    P.sD    = sD;  P.sA = sA;
    P.meta.sD0 = sD(1);  P.meta.sA0 = sA(1);
    save(outMat, '-struct', 'P');
    return
end

% fallback: solve the cut grid through the front door (slow)
nD = numel(cut.iD);  nA = numel(cut.iA);
fprintf(['fixture_dro_elfo_sheet: source sheet %s not found -- SOLVING the fixture instead: ', ...
         '%d x %d cells x %d rungs, about %d min (3 min per cell)\n'], src, nD, nA, cut.nR, 3*nD*nA);
plan = run_dro_elfo_torus([], struct('planOnly', true));
rungs = plan(1).optsA.rungs(1:cut.nR);
wd = tempname;  mkdir(wd);
run_dro_elfo_torus(1, struct('sheets', struct('target', 'elfo', 'droTau', 1, 'raanDeg', 0), ...
    'outDir', wd, 'nD', nD, 'nA', nA, 'sD0', (cut.iD(1) - 1)/6, 'sA0', (cut.iA(1) - 1)/6, ...
    'rungs', rungs, 'msHardCapSec', 0, 'keepDirect', false));
copyfile(fullfile(wd, 'dro_elfo_tau1_raan000.mat'), outMat);
rmdir(wd, 's');
end
