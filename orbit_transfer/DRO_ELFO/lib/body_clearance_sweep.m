function S = body_clearance_sweep(sheetMat, opts)
%% Purpose:
%
%   THE LOCATED CLEARANCE of every accepted entry of a DRO -> ELFO sheet,
%   above the Moon and above the Earth, with the time of each minimum. A
%   copy of GTO_DRO/lib/earth_clearance_sweep.m with two changes: both
%   bodies are measured, and the number is the LOCATED minimum (shared
%   costate_common/located_min_altitude: the range-rate root on the true
%   min-time dynamics, so a perilune pass that falls between two flight
%   samples is caught), not the sample minimum.
%
%   Every OK entry is flown (fly_transfer, pumpkyn tfMinProp), measured
%   twice (Moon centre [1-mu 0 0], radius 1737.4 km; Earth centre [-mu 0 0],
%   radius 6378.137 km), and the result goes to a resumable sidecar
%   <tag>_clearance.mat saved after every entry. With opts.demote (default)
%   and a complete census, entries below a floor are demoted to OK = false
%   after the sheet is backed up ONCE to <sheet>.bak_clearance (never
%   overwritten), and the sheet is written atomically (temp file + move).
%
%   SIDECAR RULE (copied): records are keyed on Z8 IDENTITY only, never on
%   OK. Z8SWEPT holds the z8 each altitude was measured for; on load an
%   entry stays done only where the sheet's current Z8 column is isequaln
%   to the stored one, so a re-solved entry is re-flown and a stale record
%   is never reused. The altitudes are stored RAW (independent of the
%   floors), so the floors are applied at census time and changing a floor
%   needs no new flights. A minimum that cannot be located (a flight that
%   throws or is non-finite, or located_min_altitude's nFailed > 0) is
%   stored as NaN and counts as BELOW the floor.
%
%   COUNTS: nBelow / nNonFinite come from the record (done and Z8 matching
%   the sheet), whether the entry is still OK (demoted now) or OK = false
%   (demoted by an earlier call), so a repeated call keeps the record.
%   nDone, nTodo, minMoonKm and minEarthKm are defined on the currently OK
%   entries. Demotion reloads the sheet and refuses
%   (body_clearance_sweep:sheetChanged) if its Z8 differs from what the
%   census flew. Every other variable of the sheet file (sheet_key
%   included) is carried through the rewrite.
%
%% Inputs:
%
%  sheetMat                 char                    Sheet .mat written by
%                                                   thrust_ladder_library
%                                                   (OK, Z8, rungs, sD, sA,
%                                                   meta)
%
%  opts                     struct (optional)       .demote [true]
%                                                   .floorMoonKm
%                                                   [meta.floorKm, else 500]
%                                                   .floorEarthKm [200]
%                                                   ([] or -inf = no floor)
%                                                   .sideMat
%                                                   [<sheet>_clearance.mat]
%                                                   .batchSec [inf]
%                                                   .maxEntries [inf]
%                                                   .logFile ['' = stdout]
%                                                   .propagate [true]
%                                                   locate on the true
%                                                   dynamics (false =
%                                                   Hermite on the samples)
%                                                   (.flyFcn TEST SEAM
%                                                   [@fly_transfer], handle
%                                                   (z8, rv0, rvf, phys) ->
%                                                   struct .t [M x 1], .Y
%                                                   [M x >= 6])
%
%% Outputs:
%
%  S                        struct                  .done .nDone .nTodo
%                                                   .nBelow .nBelowMoon
%                                                   .nBelowEarth .nNonFinite
%                                                   .minMoonKm .minEarthKm
%                                                   .demoted [k x 3] (iD,
%                                                   iA, kr) demoted in THIS
%                                                   call .nFlown (flights in
%                                                   this call) .floorMoonKm
%                                                   .floorEarthKm .sideMat
%                                                   .backup
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    % self-demo on a copy of the matched-experiment sheet, 5 entries, no demotion
    here = fileparts(fileparts(mfilename('fullpath')));
    src  = fullfile(here, 'results', 'compare', 'dro_elfo_tau1_raan000.mat');
    if ~isfile(src)
        fprintf('body_clearance_sweep: self-demo needs %s\n', src);  S = [];  return
    end
    wd = tempname;  mkdir(wd);
    sheetMat = fullfile(wd, 'dro_elfo_tau1_raan000.mat');  copyfile(src, sheetMat);
    opts = struct('maxEntries', 5, 'demote', false);
end
if nargin < 2, opts = struct(); end
Q = load(sheetMat);
meta = Q.meta;
floorMoonDef = 500;
if isfield(meta, 'floorKm') && ~isempty(meta.floorKm), floorMoonDef = meta.floorKm; end
 floorMoonKm = nofloor(fdef(opts, 'floorMoonKm', floorMoonDef));
floorEarthKm = nofloor(fdef(opts, 'floorEarthKm', 200));
     sideMat = fdef(opts, 'sideMat', [sheetMat(1:end-4) '_clearance.mat']);
    batchSec = fdef(opts, 'batchSec', inf);
  maxEntries = fdef(opts, 'maxEntries', inf);
      demote = fdef(opts, 'demote', true);
     logFile = fdef(opts, 'logFile', '');
   propagate = fdef(opts, 'propagate', true);
      flyFcn = fdef(opts, 'flyFcn', @fly_transfer);
 moonRadius  = 1737.4;                       % km, the shared library's Moon
 earthRadius = 6378.137;                     % km, the shared library's Earth
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

%% Sidecar (resume; Z8 identity per entry):
[nD, nA, nR] = size(Q.OK);
MOONKM = nan(nD, nA, nR);  EARTHKM = MOONKM;  TMOON = MOONKM;  TEARTH = MOONKM;
done = false(nD, nA, nR);  Z8SWEPT = nan(8, nD, nA, nR);
if isfile(sideMat)
    E = load(sideMat);
    if all(isfield(E, {'MOONKM', 'EARTHKM', 'TMOON', 'TEARTH', 'done', 'Z8SWEPT'})) ...
            && isequal(size(E.MOONKM), [nD nA nR]) && isequal(size(E.Z8SWEPT), size(Q.Z8))
        MOONKM = E.MOONKM;  EARTHKM = E.EARTHKM;  TMOON = E.TMOON;  TEARTH = E.TEARTH;
        done = E.done;  Z8SWEPT = E.Z8SWEPT;
        nStale = 0;
        for kk = 1:numel(done)
            if done(kk) && ~isequaln(Q.Z8(:, kk), Z8SWEPT(:, kk))
                done(kk) = false;  Z8SWEPT(:, kk) = NaN;  nStale = nStale + 1;
                MOONKM(kk) = NaN;  EARTHKM(kk) = NaN;  TMOON(kk) = NaN;  TEARTH(kk) = NaN;
            end
        end
        if nStale > 0, lg('clearance sweep: %d stale record(s) (z8 changed) dropped, re-flown', nStale); end
    else
        lg('clearance sweep: sidecar %s does not fit this sheet -- started over', sideMat);
    end
end

[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
dep = phase_state(tD, rvD);  arr = phase_state(tA, rvA);
moonC  = [1 - meta.muStar, 0, 0];
earthC = [-meta.muStar, 0, 0];

tAll = tic;  nFlown = 0;  stopNow = false;
for iD = 1:nD
    for iA = 1:nA
        for kr = 1:nR
            if stopNow, break, end
            if ~Q.OK(iD, iA, kr) || done(iD, iA, kr), continue, end
            if toc(tAll) > batchSec || nFlown >= maxEntries, stopNow = true;  break, end
            TN   = Q.rungs(kr);
            ndp  = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
            phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', meta.muStar, ...
                          'lStar', meta.lStar, 'tStar', meta.tStar, 'm0kg', meta.m0kg);
            rv0 = dep(Q.sD(iD));  rvf = arr(Q.sA(iA));
            hM = NaN;  hE = NaN;  tM = NaN;  tE = NaN;  why = '';
            try
                F = flyFcn(Q.Z8(:, iD, iA, kr), rv0(:), rvf(:), phys);
                if ~all(isfinite(F.Y(:, 1:6)), 'all')
                    why = 'non-finite flight';
                else
                    lo = struct('lStar', meta.lStar);
                    if propagate, lo.phys = phys; end
                    Cm = located_min_altitude(F.t(:), F.Y, moonC, moonRadius, -inf, lo);
                    Ce = located_min_altitude(F.t(:), F.Y, earthC, earthRadius, -inf, lo);
                    if Cm.nFailed == 0, hM = Cm.hMinKm;  tM = Cm.tAtMin; else, why = Cm.reason; end
                    if Ce.nFailed == 0, hE = Ce.hMinKm;  tE = Ce.tAtMin; else, why = [why ' ' Ce.reason]; end
                end
            catch err
                why = ['flight threw: ' err.message];
            end
            MOONKM(iD, iA, kr) = hM;  EARTHKM(iD, iA, kr) = hE;
            TMOON(iD, iA, kr)  = tM;  TEARTH(iD, iA, kr)  = tE;
            done(iD, iA, kr) = true;  Z8SWEPT(:, iD, iA, kr) = Q.Z8(:, iD, iA, kr);
            nFlown = nFlown + 1;
            tag = '';
            if ~(isfinite(hM) && isfinite(hE)), tag = ['NOT LOCATED (counts as below): ' why];
            elseif hM < floorMoonKm || hE < floorEarthKm, tag = '< FLOOR';
            end
            lg('  (%2d,%2d) T=%5.2f  Moon %.1f km at %.5f, Earth %.1f km at %.5f %s', ...
               iD, iA, TN, hM, tM, hE, tE, tag);
            save_side(sideMat, MOONKM, EARTHKM, TMOON, TEARTH, done, Z8SWEPT);
        end
    end
end
save_side(sideMat, MOONKM, EARTHKM, TMOON, TEARTH, done, Z8SWEPT);   % always leave a sidecar

%% Census (floors applied to the raw record):
okMask     = Q.OK;
todo       = okMask & ~done;
belowMoon  = done & ~(MOONKM >= floorMoonKm);                 % NaN counts as below
belowEarth = done & ~(EARTHKM >= floorEarthKm);
below      = belowMoon | belowEarth;
cur        = okMask & done;
S = struct('done', ~any(todo(:)), 'nDone', nnz(cur), 'nTodo', nnz(todo), ...
           'nBelow', nnz(below), 'nBelowMoon', nnz(belowMoon), 'nBelowEarth', nnz(belowEarth), ...
           'nNonFinite', nnz(done & ~(isfinite(MOONKM) & isfinite(EARTHKM))), ...
           'minMoonKm', minor_nan(MOONKM(cur & isfinite(MOONKM))), ...
           'minEarthKm', minor_nan(EARTHKM(cur & isfinite(EARTHKM))), ...
           'demoted', zeros(0, 3), 'nFlown', nFlown, 'floorMoonKm', floorMoonKm, ...
           'floorEarthKm', floorEarthKm, 'sideMat', sideMat, 'backup', '');
lg('clearance sweep: %d done, %d to do, %d below (Moon %d / %.0f km, Earth %d / %.0f km), min Moon %.1f km, min Earth %.1f km', ...
   S.nDone, S.nTodo, S.nBelow, S.nBelowMoon, floorMoonKm, S.nBelowEarth, floorEarthKm, ...
   S.minMoonKm, S.minEarthKm);

%% Demotion (complete census only; backup once; atomic rewrite):
if demote && ~S.done
    lg('clearance sweep: census incomplete (%d to do) -- no demotion this call', S.nTodo);
elseif demote && any(below(:) & okMask(:))
    R = load(sheetMat);                                       % fresh load, every variable
    if ~isequaln(R.Z8, Q.Z8)
        error('body_clearance_sweep:sheetChanged', 'sheet %s changed during the census; rerun', sheetMat);
    end
    bak = [sheetMat(1:end-4) '.bak_clearance'];
    if ~isfile(bak), copyfile(sheetMat, bak); end
    S.backup = bak;
    dem = below & R.OK;
    [dI, dJ, dK] = ind2sub(size(dem), find(dem));
    S.demoted = [dI, dJ, dK];
    R.OK(dem) = false;
    R.meta.clearanceSweep = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
        'floorMoonKm', floorMoonKm, 'floorEarthKm', floorEarthKm, ...
        'nDemoted', nnz(below), 'demoted', S.demoted, 'sideMat', sideMat);
    tmp = sprintf('%s_tmp%06d.mat', sheetMat(1:end-4), randi(999999));
    save(tmp, '-struct', 'R');
    movefile(tmp, sheetMat, 'f');
    lg('clearance sweep: %d entr%s demoted (OK = false), backup %s', size(S.demoted, 1), ...
       pick(size(S.demoted, 1) == 1, 'y', 'ies'), bak);
end
if nargin == 0, disp(S), end
end

% ------------------------------------------------------------------------
function save_side(sideMat, MOONKM, EARTHKM, TMOON, TEARTH, done, Z8SWEPT)
% SAVE_SIDE  Write the sidecar atomically (temp file + move).
% INPUTS: sideMat [1 x n] char; MOONKM, EARTHKM, TMOON, TEARTH [nD x nA x nR];
%         done logical [nD x nA x nR]; Z8SWEPT [8 x nD x nA x nR].
% OUTPUTS: none (file written).
tmp = sprintf('%s_tmp%06d.mat', sideMat(1:end-4), randi(999999));
save(tmp, 'MOONKM', 'EARTHKM', 'TMOON', 'TEARTH', 'done', 'Z8SWEPT');
movefile(tmp, sideMat, 'f');
end

function f = nofloor(f)
% NOFLOOR  [] means no floor (-inf).  INPUTS: f [] or [1 x 1] km.
% OUTPUTS: f [1 x 1] km.
if isempty(f), f = -inf; end
end

function m = minor_nan(v)
% MINOR_NAN  min of v, NaN when v is empty.  INPUTS: v [k x 1] (k >= 0).
% OUTPUTS: m [1 x 1].
if isempty(v), m = NaN; else, m = min(v(:)); end
end

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s [1 x 1] struct; f [1 x n]
% char; v0 [any].  OUTPUTS: v [any].
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c [1 x 1] logical; a, b [any].
% OUTPUTS: v [any] (a when c, else b).
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f [1 x n]
% char log path ('' = stdout); s [1 x m] char.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
