function results = audit_gto_entries(catDirIn)
%% Purpose:
%
%   TASK 5 SPOT-AUDIT: end-to-end replay of 3 solved GTO -> tulip costate
%   entries, picked PROGRAMMATICALLY from the sheet .mats under
%   results/catalog/ (engine-native OK/TF/Z8/ATT/rungs/sD/sA/meta arrays --
%   the packaged catalog format does not exist for this campaign yet, so
%   this audit reads the sheets thrust_ladder_library itself writes).
%
%   For each picked entry:
%     1. rebuild BOTH endpoints from THE SHEET'S OWN meta.depParams /
%        meta.arrParams via costate_common/get_family_orbit -- this proves
%        sheet-local reconstruction (a future packager must not homogenize
%        orientations away);
%     2. fly the stored z8 with pumpkyn.cr3bp.tfMinProp and ASSERT the
%        flown arrival miss < 100 km (the campaign's own flown-arrival
%        gate, gateKm);
%     3. hand z8 to pumpkyn.cr3bp.tfMin and ASSERT it returns the SAME
%        point, |zAccepted - z8| < 1e-6 (the campaign's own acceptance
%        tolerance, accTol);
%     4. report (not gate) the entry's minimum Earth altitude over the
%        flown arc -- Earth held fixed at (-muStar,0,0) ND, radius 6378 km.
%
%   PICK RULE (GPT review, Stage B, no hard-coded coordinates): candidates
%   are every (sheet, iD, iA, kr) triple with OK==true across the 16
%   sheets, walked in a DETERMINISTIC order (sheet name, then iD, iA, kr)
%   to a 3-entry set spanning >=2 DEPARTURE ORIENTATIONS
%   (meta.depParams.orientDeg) and >=2 THRUST RUNGS, with at least one pick
%   at a LOW rung (5 N or 7 N). See the pick_three() subfunction.
%
%% Inputs:
%
%  catDirIn                 char                    (optional) sheet
%                                                   directory override
%                                                   [results/catalog next
%                                                   to this file]
%
%% Outputs:
%
%  results                  [3 x 1] struct          One row per picked
%                                                   entry: .sheet .iD .iA
%                                                   .rungN .tf_nd .flownKm
%                                                   .acceptDz
%                                                   .minEarthAltKm
%                                                   .passFly .passAccept
%
%% Revision History:
%  M. Casey                                                   (c) 08/29/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(catDirIn)
    catDirIn = fullfile(here, 'results', 'catalog');
end
setup_paths();

earthRadiusKm = 6378;                % km, spherical Earth
gateKm        = 100;                 % campaign flown-arrival gate (km)
accTol        = 1e-6;                % campaign tfMin acceptance tolerance

%% 1. The 16 official sheets, named exactly as run_gto_catalog.m builds
%% them (never a glob -- a stray backup like *_pre_densify.mat must not
%% enter the candidate pool):
orients = [0 90 180 270];
Nps     = [5 7 9 12];
nSheets = numel(orients)*numel(Nps);
sheetName = cell(nSheets,1);
ks = 0;
for od = orients
    for np = Nps
        ks = ks + 1;
        sheetName{ks} = sprintf('gto_o%03d_Np%d.mat', od, np);
    end
end

%% 2. Load each sheet's bookkeeping (OK/rungs/meta) once: ------------------
sheetOK    = cell(nSheets,1);
sheetRungs = cell(nSheets,1);
sheetMeta  = cell(nSheets,1);
sheetPath  = cell(nSheets,1);
totalCells = 0;
for ks = 1:nSheets
    sp = fullfile(catDirIn, sheetName{ks});
    assert(isfile(sp), 'audit_gto_entries:missingSheet', ...
        'expected sheet %s not found in %s', sheetName{ks}, catDirIn);
    Q = load(sp, 'OK', 'rungs', 'meta');
    sheetOK{ks}    = Q.OK;
    sheetRungs{ks} = Q.rungs;
    sheetMeta{ks}  = Q.meta;
    sheetPath{ks}  = sp;
    totalCells     = totalCells + numel(Q.OK);
end

%% 3. Flatten every solved (sheet,iD,iA,kr) triple into candidate arrays,
%% preallocated to the exact upper bound so nothing grows in the loop:
candSheet  = zeros(totalCells,1);
candOrient = zeros(totalCells,1);
candRungN  = zeros(totalCells,1);
candID     = zeros(totalCells,1);
candIA     = zeros(totalCells,1);
candKr     = zeros(totalCells,1);
cnt = 0;
for ks = 1:nSheets
    OKs = sheetOK{ks};
    rungsHere = sheetRungs{ks};
    orientDeg = sheetMeta{ks}.depParams.orientDeg;
    [nD, nA, nR] = size(OKs);
    for iD = 1:nD
        for iA = 1:nA
            for kr = 1:nR
                if ~OKs(iD,iA,kr), continue, end
                cnt = cnt + 1;
                candSheet(cnt)  = ks;
                candOrient(cnt) = orientDeg;
                candRungN(cnt)  = rungsHere(kr);
                candID(cnt)     = iD;
                candIA(cnt)     = iA;
                candKr(cnt)     = kr;
            end
        end
    end
end
candSheet  = candSheet(1:cnt);
candOrient = candOrient(1:cnt);
candRungN  = candRungN(1:cnt);
candID     = candID(1:cnt);
candIA     = candIA(1:cnt);
candKr     = candKr(1:cnt);
assert(cnt > 0, 'audit_gto_entries:noSolvedCells', ...
    'no solved (OK) cells found across the 16 sheets');

%% 4. Pick 3 entries spanning >=2 orientations and >=2 rungs, one low-rung:
pick = pick_three(candOrient, candRungN);

%% 5. Replay each pick: rebuild endpoints, fly z8, re-accept, report Earth
%% clearance: --------------------------------------------------------------
results = repmat(struct('sheet','', 'iD',0, 'iA',0, 'rungN',0, ...
    'tf_nd',0, 'flownKm',0, 'acceptDz',0, 'minEarthAltKm',0, ...
    'passFly',false, 'passAccept',false), 3, 1);

for kp = 1:3
    p  = pick(kp);
    ks = candSheet(p);
    iD = candID(p);  iA = candIA(p);  kr = candKr(p);
    meta = sheetMeta{ks};

    Qz = load(sheetPath{ks}, 'Z8', 'sD', 'sA');
    z8    = Qz.Z8(:,iD,iA,kr);
    rungN = sheetRungs{ks}(kr);

    % Sheet-local endpoint reconstruction (brief step 1): the family and
    % parameters come from THIS sheet's own meta, never a shared default.
    [tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
    [tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
    rv0 = interp1(tD, rvD, mod(Qz.sD(iD),1)*tD(end), 'spline');
    rvf = interp1(tA, rvA, mod(Qz.sA(iA),1)*tA(end), 'spline');

    % Same ND thrust/exhaust conversion thrust_ladder_library uses, from
    % THIS sheet's own thruster metadata:
    g0  = 9.80665*meta.tStar^2/(1000*meta.lStar);
    cnd = (meta.ispS/meta.tStar)*g0;
    Tnd = (rungN/meta.m0kg)*meta.tStar^2/(meta.lStar*1000);

    % --- Gate 1: fly the stored z8, measure the arrival miss ------------
    [~, y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6)'; 1; z8(1:7)], ...
                Tnd, cnd, meta.muStar);
    flownKm = vecnorm(y(end,1:3) - rvf(1:3), 2, 2);
    flownKm = flownKm * meta.lStar;
    assert(flownKm < gateKm, 'audit_gto_entries:flownMiss', ...
        '%s (iD=%d,iA=%d,kr=%d): flown miss %.2f km >= gate %.0f km', ...
        sheetName{ks}, iD, iA, kr, flownKm, gateKm);

    % --- Gate 2: hand z8 to tfMin, require an unchanged accepted point --
    % (evalc suppresses fsolve's own 'Display','iter' chatter, exactly as
    % thrust_ladder_library does at its own acceptance step)
    evalc(['zAccepted = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z8, ' ...
           'Tnd, cnd, meta.muStar);']);
    acceptDz = vecnorm(zAccepted - z8);
    assert(acceptDz < accTol, 'audit_gto_entries:acceptDrift', ...
        '%s (iD=%d,iA=%d,kr=%d): |zAccepted-z8| %.2e >= tol %.1e', ...
        sheetName{ks}, iD, iA, kr, acceptDz, accTol);

    % --- Report-only: minimum Earth altitude over the flown arc ---------
    earthPos = [-meta.muStar, 0, 0];
    altKm = vecnorm(y(:,1:3) - earthPos, 2, 2)*meta.lStar - earthRadiusKm;
    minEarthAltKm = min(altKm);

    results(kp) = struct('sheet', sheetName{ks}, 'iD', iD, 'iA', iA, ...
        'rungN', rungN, 'tf_nd', z8(8), 'flownKm', flownKm, ...
        'acceptDz', acceptDz, 'minEarthAltKm', minEarthAltKm, ...
        'passFly', flownKm < gateKm, 'passAccept', acceptDz < accTol);
end

%% 6. Print the 3-row table: ------------------------------------------------
fprintf('\n%-20s %-9s %6s %10s %10s %10s %14s\n', ...
    'sheet', 'cell', 'rung', 'tf_nd', 'flownKm', 'acceptDz', 'minEarthAltKm');
for kp = 1:3
    r = results(kp);
    fprintf('%-20s (%2d,%2d)  %4gN %10.5f %10.3f %10.2e %14.1f\n', ...
        r.sheet, r.iD, r.iA, r.rungN, r.tf_nd, r.flownKm, r.acceptDz, ...
        r.minEarthAltKm);
end
fprintf('\nall %d rows: flown < %g km = %d, |zAccepted-z8| < %.0e = %d\n', ...
    numel(results), gateKm, all([results.passFly]), accTol, ...
    all([results.passAccept]));
end

% ------------------------------------------------------------------------
function pick = pick_three(candOrient, candRungN)
% PICK_THREE  Deterministic 3-of-N pick spanning >=2 orientations and >=2
% rung classes, with at least one LOW rung (5 or 7 N). Candidates are
% already in the caller's canonical (sheet, iD, iA, kr) order -- "first"
% below always means first in that fixed order, so the pick never depends
% on anything but the sheet contents themselves.
%
% INPUTS: candOrient [cnt x 1]; candRungN [cnt x 1] (parallel arrays).
% OUTPUTS: pick [3 x 1] indices into the candidate arrays.
cnt = numel(candOrient);
isLow = (candRungN == 5) | (candRungN == 7);

idx1 = find(isLow, 1, 'first');
assert(~isempty(idx1), 'audit_gto_entries:noLowRung', ...
    'no solved entry found at a low rung (5 N or 7 N)');

idx2 = find(candOrient ~= candOrient(idx1), 1, 'first');
assert(~isempty(idx2), 'audit_gto_entries:noSecondOrient', ...
    'no solved entry found at an orientation different from the first pick');

taken = false(cnt,1);  taken([idx1 idx2]) = true;
rungsSoFar = [candRungN(idx1); candRungN(idx2)];
if numel(unique(rungsSoFar)) >= 2
    % rung diversity already satisfied by the first two picks; prefer a
    % THIRD distinct orientation as a bonus, else any remaining candidate:
    orientsSoFar = [candOrient(idx1); candOrient(idx2)];
    idx3 = find(~taken & ~ismember(candOrient, orientsSoFar), 1, 'first');
    if isempty(idx3)
        idx3 = find(~taken, 1, 'first');
    end
else
    idx3 = find(~taken & candRungN ~= rungsSoFar(1), 1, 'first');
    assert(~isempty(idx3), 'audit_gto_entries:noThirdRung', ...
        'no solved entry found with a rung class different from the first two picks');
end

pick = [idx1; idx2; idx3];
assert(numel(unique(candOrient(pick))) >= 2, ...
    'audit_gto_entries:pickOrient', 'picked set spans < 2 orientations');
assert(numel(unique(candRungN(pick))) >= 2, ...
    'audit_gto_entries:pickRung', 'picked set spans < 2 rung classes');
assert(any(candRungN(pick) == 5 | candRungN(pick) == 7), ...
    'audit_gto_entries:pickLow', 'picked set has no low-rung (5/7 N) entry');
end
