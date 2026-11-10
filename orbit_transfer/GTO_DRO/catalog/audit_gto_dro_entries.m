function results = audit_gto_dro_entries(catDirIn, sheetMats, floorEarthKm)
%% Purpose:
%
%   SPOT-AUDIT of the GTO -> DRO costate catalog: end-to-end replay of up to
%   3 solved entries, picked PROGRAMMATICALLY from the sheet .mats (the
%   engine-native OK/TF/Z8/rungs/sD/sA/meta arrays), the DRO-arrival twin of
%   GTO_tulip/catalog/audit_gto_entries.m.
%
%   For each picked entry:
%     1. rebuild BOTH endpoints from THE SHEET'S OWN meta.depParams /
%        meta.arrParams via costate_common/get_family_orbit (proves
%        sheet-local reconstruction);
%     2. fly the stored z8 with pumpkyn.cr3bp.tfMinProp and ASSERT the
%        flown arrival miss < 100 km (the campaign's gateKm);
%     3. hand z8 to pumpkyn.cr3bp.tfMin and ASSERT it returns the SAME
%        point, |zAccepted - z8| < 1e-6 (the campaign's accTol);
%     4. GATE the entry's minimum Earth altitude over the flown arc against
%        floorEarthKm (assert audit_gto_dro_entries:earthFloor). This
%        closes, for the audited entries, the certifier's missing
%        propagated Earth clearance (earth_clearance_sweep closes it for
%        every entry). The gated value is the MINIMUM of the replay's own
%        altitude (Earth fixed at (-muStar,0,0) ND, radius 6378 km) and
%        the sweep's sidecar value (<sheet>_earth.mat EARTHKM, used only
%        where the sidecar's Z8SWEPT equals the sheet's current z8); a
%        non-finite altitude fails the gate.
%
%   Like the tulip audit the three gates are asserts, so a returned row has
%   passFly = passAccept = passEarth = true; the flags exist so the chain
%   can carry the same fields whatever the audit does.
%
%   PICK RULE. When the sheets can satisfy the tulip rule (3 entries
%   spanning >= 2 orientations and >= 2 rung classes, at least one at 5 or
%   7 N) that rule is used unchanged. Otherwise the audit DEGRADES to up to
%   3 distinct OK entries preferring distinct rungs and distinct sheets
%   (deterministic walk in (sheet, iD, iA, kr) order), and prints that the
%   reduced rule was used (results(1).rule = 'reduced').
%
%% Inputs:
%
%  catDirIn                 char                    (optional) sheet folder
%                                                   [../results/catalog];
%                                                   used only to build the
%                                                   default sheet list
%
%  sheetMats                cell                    (optional) explicit
%                                                   sheet .mat paths; [] =
%                                                   the 16 official names
%                                                   under catDirIn that
%                                                   exist
%
%  floorEarthKm             double                  (optional) Earth
%                                                   altitude floor [200]
%
%% Outputs:
%
%  results                  [n x 1] struct          n = 1..3 rows: .sheet
%                                                   .iD .iA .rungN .tf_nd
%                                                   .flownKm .acceptDz
%                                                   .minEarthAltKm
%                                                   .earthSweepKm (NaN when
%                                                   no valid sidecar value)
%                                                   .passFly .passAccept
%                                                   .passEarth .rule
%                                                   .kr (rung index) and
%                                                   .z8 [8 x 1] (the
%                                                   replayed entry: the
%                                                   chain stores them as
%                                                   the audit's fingerprint
%                                                   of the sheets)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  fix round 1: rows carry .kr and .z8 (the chain's audit fingerprint)
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
if nargin < 1 || isempty(catDirIn), catDirIn = fullfile(fileparts(here), 'results', 'catalog'); end
if nargin < 3 || isempty(floorEarthKm), floorEarthKm = 200; end
if nargin < 2 || isempty(sheetMats)
    sheetMats = {};
    for od = [0 90 180 270]
        for ta = [0.5 1 2 3]
            p = fullfile(catDirIn, sprintf('gto_dro_o%03d_tau%g.mat', od, ta));
            if isfile(p), sheetMats{end+1} = p; end
        end
    end
end

earthRadiusKm = 6378;                % km, spherical Earth
gateKm        = 100;                 % campaign flown-arrival gate (km)
accTol        = 1e-6;                % campaign tfMin acceptance tolerance

%% 1. Candidates: every OK (sheet, iD, iA, kr), deterministic order:
nSheets = numel(sheetMats);
sheetName = cell(nSheets, 1);
nMax = 0;
for ks = 1:nSheets
    Q = load(sheetMats{ks}, 'OK');
    nMax = nMax + nnz(Q.OK);
end
cand = zeros(nMax, 4);  candOrient = zeros(nMax, 1);  candRungN = zeros(nMax, 1);
cnt = 0;
for ks = 1:nSheets
    Q = load(sheetMats{ks}, 'OK', 'rungs', 'meta');
    [~, sheetName{ks}] = fileparts(sheetMats{ks});
    [nD, nA, nR] = size(Q.OK);
    for iD = 1:nD
        for iA = 1:nA
            for kr = 1:nR
                if ~Q.OK(iD, iA, kr), continue, end
                cnt = cnt + 1;
                cand(cnt, :)     = [ks iD iA kr];
                candOrient(cnt)  = Q.meta.depParams.orientDeg;
                candRungN(cnt)   = Q.rungs(kr);
            end
        end
    end
end
assert(cnt > 0, 'audit_gto_dro_entries:noSolvedCells', ...
    'no solved (OK) cells found across the %d sheet(s)', nSheets);
cand = cand(1:cnt, :);  candOrient = candOrient(1:cnt);  candRungN = candRungN(1:cnt);

%% 2. Pick: the tulip rule when satisfiable, else the reduced rule:
pick = pick_tulip_rule(candOrient, candRungN);
rule = 'tulip';
if isempty(pick)
    pick = pick_reduced(cand(:, 1), candRungN);
    rule = 'reduced';
    fprintf(['audit_gto_dro_entries: the tulip pick rule (5/7 N rung, 2 orientations) cannot be ', ...
             'satisfied by these sheets; reduced rule used: %d distinct OK entr%s\n'], ...
            numel(pick), pick_plural(numel(pick)));
end
nPick = numel(pick);

%% 3. Replay each pick:
results = repmat(struct('sheet', '', 'iD', 0, 'iA', 0, 'rungN', 0, 'tf_nd', 0, ...
    'flownKm', 0, 'acceptDz', 0, 'minEarthAltKm', 0, 'earthSweepKm', NaN, ...
    'passFly', false, 'passAccept', false, 'passEarth', false, 'rule', rule, ...
    'kr', 0, 'z8', zeros(8, 1)), nPick, 1);

for kp = 1:nPick
    ks = cand(pick(kp), 1);  iD = cand(pick(kp), 2);  iA = cand(pick(kp), 3);  kr = cand(pick(kp), 4);
    Q = load(sheetMats{ks}, 'Z8', 'sD', 'sA', 'rungs', 'meta');
    meta = Q.meta;
    z8 = Q.Z8(:, iD, iA, kr);
    rungN = Q.rungs(kr);

    % Sheet-local endpoint reconstruction: family and parameters come from
    % THIS sheet's own meta, never a shared default.
    [tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
    [tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
    depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
    rv0 = depAt(Q.sD(iD)).';
    rvf = arrAt(Q.sA(iA)).';

    ndp = nd_propulsion(rungN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
    cnd = ndp.cnd;   Tnd = ndp.Tnd;

    % --- Gate 1: fly the stored z8, measure the arrival miss ------------
    [~, y] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(1:6)'; 1; z8(1:7)], Tnd, cnd, meta.muStar);
    flownKm = vecnorm(y(end, 1:3) - rvf(1:3), 2, 2) * meta.lStar;
    assert(flownKm < gateKm, 'audit_gto_dro_entries:flownMiss', ...
        '%s (iD=%d,iA=%d,kr=%d): flown miss %.2f km >= gate %.0f km', ...
        sheetName{ks}, iD, iA, kr, flownKm, gateKm);

    % --- Gate 2: tfMin must return the same point -----------------------
    evalc(['zAccepted = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z8, ' ...
           'Tnd, cnd, meta.muStar);']);
    acceptDz = vecnorm(zAccepted - z8);
    assert(acceptDz < accTol, 'audit_gto_dro_entries:acceptDrift', ...
        '%s (iD=%d,iA=%d,kr=%d): |zAccepted-z8| %.2e >= tol %.1e', ...
        sheetName{ks}, iD, iA, kr, acceptDz, accTol);

    % --- Gate 3: minimum Earth altitude over the flown arc --------------
    earthPos = [-meta.muStar, 0, 0];
    altKm = vecnorm(y(:, 1:3) - earthPos, 2, 2) * meta.lStar - earthRadiusKm;
    replayKm = min(altKm);
    sweepKm = NaN;
    sideMat = [sheetMats{ks}(1:end-4) '_earth.mat'];
    if isfile(sideMat)
        Se = load(sideMat);
        if isfield(Se, 'Z8SWEPT') && isfield(Se, 'done') && Se.done(iD, iA, kr) ...
                && isequaln(Se.Z8SWEPT(:, iD, iA, kr), z8)
            sweepKm = Se.EARTHKM(iD, iA, kr);
        end
    end
    altSet = [replayKm sweepKm];
    if any(~isfinite(altSet)), minEarthAltKm = NaN; else, minEarthAltKm = min(altSet); end
    assert(isfinite(minEarthAltKm) && minEarthAltKm >= floorEarthKm, 'audit_gto_dro_entries:earthFloor', ...
        '%s (iD=%d,iA=%d,kr=%d): Earth altitude %.1f km (replay %.1f, sweep %.1f) below floor %.0f km', ...
        sheetName{ks}, iD, iA, kr, minEarthAltKm, replayKm, sweepKm, floorEarthKm);

    results(kp) = struct('sheet', sheetName{ks}, 'iD', iD, 'iA', iA, 'rungN', rungN, ...
        'tf_nd', z8(8), 'flownKm', flownKm, 'acceptDz', acceptDz, ...
        'minEarthAltKm', minEarthAltKm, 'earthSweepKm', sweepKm, ...
        'passFly', flownKm < gateKm, 'passAccept', acceptDz < accTol, ...
        'passEarth', minEarthAltKm >= floorEarthKm, 'rule', rule, 'kr', kr, 'z8', z8);
end

%% 4. Table:
fprintf('\n%-22s %-9s %6s %10s %10s %10s %14s\n', ...
    'sheet', 'cell', 'rung', 'tf_nd', 'flownKm', 'acceptDz', 'minEarthAltKm');
for kp = 1:nPick
    r = results(kp);
    fprintf('%-22s (%2d,%2d)  %4gN %10.5f %10.3f %10.2e %14.1f\n', ...
        r.sheet, r.iD, r.iA, r.rungN, r.tf_nd, r.flownKm, r.acceptDz, r.minEarthAltKm);
end
fprintf('\npick rule %s; %d rows: flown < %g km = %d, |zAccepted-z8| < %.0e = %d, Earth >= %g km = %d\n', ...
    rule, nPick, gateKm, all([results.passFly]), accTol, all([results.passAccept]), ...
    floorEarthKm, all([results.passEarth]));
end

% ------------------------------------------------------------------------
function pick = pick_tulip_rule(candOrient, candRungN)
% PICK_TULIP_RULE  The tulip audit's 3-of-N rule (>= 2 orientations, >= 2
% rung classes, one 5 or 7 N), or [] when the candidates cannot satisfy it.
% INPUTS: candOrient [cnt x 1]; candRungN [cnt x 1].  OUTPUTS: pick [3 x 1] or [].
cnt = numel(candOrient);
pick = [];
isLow = (candRungN == 5) | (candRungN == 7);
idx1 = find(isLow, 1, 'first');
if isempty(idx1), return, end
idx2 = find(candOrient ~= candOrient(idx1), 1, 'first');
if isempty(idx2), return, end
taken = false(cnt, 1);  taken([idx1 idx2]) = true;
rungsSoFar = [candRungN(idx1); candRungN(idx2)];
if numel(unique(rungsSoFar)) >= 2
    orientsSoFar = [candOrient(idx1); candOrient(idx2)];
    idx3 = find(~taken & ~ismember(candOrient, orientsSoFar), 1, 'first');
    if isempty(idx3), idx3 = find(~taken, 1, 'first'); end
else
    idx3 = find(~taken & candRungN ~= rungsSoFar(1), 1, 'first');
end
if isempty(idx3), return, end
p = [idx1; idx2; idx3];
if numel(unique(candOrient(p))) >= 2 && numel(unique(candRungN(p))) >= 2 ...
        && any(candRungN(p) == 5 | candRungN(p) == 7)
    pick = p;
end
end

function pick = pick_reduced(candSheet, candRungN)
% PICK_REDUCED  Up to 3 distinct candidates, greedy in canonical order,
% each step taking the first candidate that adds the most new rungs and new
% sheets.  INPUTS: candSheet [cnt x 1]; candRungN [cnt x 1].  OUTPUTS: pick [n x 1].
cnt = numel(candSheet);
nPickMax = min(3, cnt);
pick = zeros(nPickMax, 1);
for kk = 1:nPickMax
    score = -inf(cnt, 1);
    for kc = 1:cnt
        if any(pick(1:kk-1) == kc), continue, end
        score(kc) = double(~any(candRungN(pick(1:kk-1)) == candRungN(kc))) ...
                  + double(~any(candSheet(pick(1:kk-1)) == candSheet(kc)));
    end
    [~, kbest] = max(score);
    pick(kk) = kbest;
end
end

function s = pick_plural(n)
% PICK_PLURAL  'y' or 'ies'.  INPUTS: n.  OUTPUTS: s char.
if n == 1, s = 'y'; else, s = 'ies'; end
end
