function S = earth_clearance_sweep(sheetMat, opts)
% EARTH_CLEARANCE_SWEEP  Flown Earth altitude of every OK entry in a sheet.
%
%   The direct solve enforced the Earth floor at its nodes and midpoints,
%   and the engine's G1b gate flew each accepted z8 for the arrival miss,
%   but nothing recorded the flown arc's minimum Earth altitude. This sweep
%   closes that gap on the sheet (campaign side, no engine change): every
%   OK entry is flown with fly_transfer, earth_clearance measures its
%   minimum altitude, the result goes to a resumable sidecar, and an
%   optional writeback demotes entries below the floor to OK = false after
%   backing the sheet up once. Audit and packaging read the sidecar.
%
%   Sidecar rule: done/EARTHKM are keyed on Z8 IDENTITY only, never on OK.
%   The sidecar stores Z8SWEPT, the z8 each altitude was measured for; on
%   load an entry stays done only where the sheet's current Z8 column is
%   isequaln to the stored one, so an entry re-solved by the engine's redo
%   is re-flown. Entries this sweep demoted keep their record (OK = false
%   but Z8 unchanged). The census masks (nDone, nTodo, nBelow, minAltKm)
%   use the sheet's current OK. A non-finite altitude (failed or garbage
%   flight, or fly_transfer throwing) on an OK entry counts as BELOW the
%   floor. The altitude is the minimum over the flown output SAMPLES
%   (earth_clearance does not interpolate), so a dip between samples is
%   not seen: adequate where the margin is large (demo cell 13,888 km
%   against 200 km), a known limit near the floor.
%
%   nBelow / demoted / nNonFinite come from the Z8-keyed sidecar record, not
%   from the current OK: every entry whose record is done (Z8 matches the
%   sheet) and whose altitude is below the floor or non-finite, whether the
%   sheet has it OK (demoted by this writeback) or OK = false (demoted by
%   an earlier writeback), so repeated writebacks keep the record. Entries
%   that are OK = false and never swept stay out. nDone, nTodo and
%   minAltKm are defined on the currently-OK entries. Writeback reloads the
%   sheet, refuses (earth_clearance_sweep:sheetChanged) if its Z8 differs
%   from what the census flew, and applies the demotions to the fresh load.
%
% INPUTS:
%   sheetMat - sheet .mat written by thrust_ladder_library        [char]
%   opts     - (optional) .floorKm [meta.floorEarthKm, or 200 if the sheet
%              has none], .sideMat [<sheet>_earth.mat], .batchSec [inf],
%              .maxEntries [inf], .writeback [false], .logFile ['']
%
% OUTPUTS:
%   S        - .done, .nDone, .nTodo, .nBelow, .minAltKm, .sideMat,
%              .demoted [k x 3] (iD, iA, kr) when written back, .floorKm,
%              .nFlown (flights done in THIS call), .nNonFinite
%
% M. Casey                                                   (c) 11/10/2026
% Copyright Coorbital Inc.

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'));
    sheetMat = fixture_demo_sheet();  opts = struct();
end
if nargin < 2, opts = struct(); end
Q = load(sheetMat);
meta = Q.meta;
floorDef = 200;
if isfield(meta, 'floorEarthKm') && ~isempty(meta.floorEarthKm), floorDef = meta.floorEarthKm; end
floorKm    = fdef(opts, 'floorKm', floorDef);
sideMat    = fdef(opts, 'sideMat', [sheetMat(1:end-4) '_earth.mat']);
batchSec   = fdef(opts, 'batchSec', inf);
maxEntries = fdef(opts, 'maxEntries', inf);
writeback  = fdef(opts, 'writeback', false);
logFile    = fdef(opts, 'logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

[nD, nA, nR] = size(Q.OK);
EARTHKM = nan(nD, nA, nR);  done = false(nD, nA, nR);
Z8SWEPT = nan(8, nD, nA, nR);
if isfile(sideMat)
    E = load(sideMat);
    if isfield(E, 'Z8SWEPT') && isequal(size(E.EARTHKM), [nD nA nR]) ...
            && isequal(size(E.Z8SWEPT), size(Q.Z8)) && E.floorKm == floorKm
        EARTHKM = E.EARTHKM;  done = E.done;  Z8SWEPT = E.Z8SWEPT;
        for kk = 1:numel(done)
            if done(kk) && ~isequaln(Q.Z8(:, kk), Z8SWEPT(:, kk))
                done(kk) = false;  EARTHKM(kk) = NaN;  Z8SWEPT(:, kk) = NaN;
            end
        end
    end
end

[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
dep = phase_state(tD, rvD);  arr = phase_state(tA, rvA);
pEarth = struct('muStar', meta.muStar, 'lStar', meta.lStar);

tAll = tic;  nThis = 0;  stopNow = false;  nFlown = 0;
for iD = 1:nD
    for iA = 1:nA
        for kr = 1:nR
            if stopNow, break, end
            if ~Q.OK(iD, iA, kr) || done(iD, iA, kr), continue, end
            if toc(tAll) > batchSec || nThis >= maxEntries, stopNow = true;  break, end
            TN  = Q.rungs(kr);
            ndp = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
            phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', meta.muStar, ...
                          'lStar', meta.lStar, 'tStar', meta.tStar, 'm0kg', meta.m0kg);
            rv0 = dep(Q.sD(iD));  rvf = arr(Q.sA(iA));
            try
                F = fly_transfer(Q.Z8(:, iD, iA, kr), rv0(:), rvf(:), phys);
                hMin = earth_clearance(F.t(:), F.Y, floorKm, pEarth);
                if any(~isfinite(F.Y(:, 1:3)), 'all'), hMin = NaN; end
            catch err
                hMin = NaN;
                lg('  (%2d,%2d) T=%5.1f  flight threw: %s', iD, iA, TN, err.message);
            end
            EARTHKM(iD, iA, kr) = hMin;  done(iD, iA, kr) = true;
            Z8SWEPT(:, iD, iA, kr) = Q.Z8(:, iD, iA, kr);
            nThis = nThis + 1;  nFlown = nFlown + 1;
            if ~isfinite(hMin), tag = 'NON-FINITE (counts as below)';
            elseif hMin < floorKm, tag = '< FLOOR';
            else, tag = '';
            end
            lg('  (%2d,%2d) T=%5.1f  min Earth alt %.1f km %s', iD, iA, TN, hMin, tag);
            save(sideMat, 'EARTHKM', 'done', 'floorKm', 'Z8SWEPT');
        end
    end
end
save(sideMat, 'EARTHKM', 'done', 'floorKm', 'Z8SWEPT');   % always leave a sidecar

okMask  = Q.OK;
todo    = okMask & ~done;
below   = done & ~(EARTHKM >= floorKm);                   % record-based; NaN counts as below
altOK   = EARTHKM(okMask & done & isfinite(EARTHKM));
S = struct('done', ~any(todo(:)), 'nDone', nnz(okMask & done), 'nTodo', nnz(todo), ...
           'nBelow', nnz(below), 'minAltKm', min(altOK, [], 'all'), ...
           'sideMat', sideMat, 'demoted', zeros(0, 3), 'floorKm', floorKm, ...
           'nFlown', nFlown, 'nNonFinite', nnz(done & ~isfinite(EARTHKM)));
if isempty(S.minAltKm), S.minAltKm = NaN; end
lg('earth sweep: %d done, %d to do, %d below %.0f km, min %.1f km', S.nDone, S.nTodo, S.nBelow, floorKm, S.minAltKm);

if writeback
    assert(S.done, 'earth_clearance_sweep:incomplete', 'writeback needs a complete census (%d to do)', S.nTodo);
    bak = [sheetMat(1:end-4) '.bak_earth'];
    R = load(sheetMat);                                   % fresh load, not the stale copy
    if ~isequaln(R.Z8, Q.Z8)
        error('earth_clearance_sweep:sheetChanged', 'sheet %s changed during the census; rerun', sheetMat);
    end
    if ~isfile(bak), copyfile(sheetMat, bak); end
    [dI, dJ, dK] = ind2sub(size(below), find(below));
    S.demoted = [dI, dJ, dK];
    Q = R;
    Q.OK(below) = false;
    Q.meta.earthSweep = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'floorKm', floorKm, ...
                               'nDemoted', nnz(below), 'demoted', S.demoted);
    save(sheetMat, '-struct', 'Q');
    lg('earth sweep writeback: %d entr%s demoted, backup %s', nnz(below), pick(nnz(below) == 1, 'y', 'ies'), bak);
end
if nargin == 0, disp(S), end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f; s.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
