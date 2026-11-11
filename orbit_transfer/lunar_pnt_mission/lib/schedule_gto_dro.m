function S = schedule_gto_dro(tf, sD, sA, sheetOrient, tauG, tauD, sgn, opts)
%% Purpose:
%
%   Schedule leg 1 of the lunar PNT mission: six transfer stages leaving
%   ONE inertially fixed GTO and arriving on the six DRO slots, from the
%   catalog min-time table tf(iD, iA, iS) (departure phase sD(iD), arrival
%   phase sA(iA), sheet iS at GTO orientation sheetOrient(iS)) and whole-
%   GTO-period waits. A pure function of tables: no file access, no solver.
%
%   GLOBAL TRIPLE (orient0, gtoPhase0, delta) is searched over opts grids:
%   orient0 [deg] is the rotating-frame orientation of the inertial GTO at
%   leg-1 time 0, gtoPhase0 [fraction] the spacecraft's GTO phase at time
%   0, delta [fraction] the DRO phase offset of the slot grid (slot j is at
%   DRO phase mod(delta + (j-1)/6 + t/tauD, 1)).
%
%   CANDIDATE RULE for slot j, cell (iD, iA, iS) and wait n in 0:nMaxRevs:
%       tDep      = tauG*(mod(sD(iD) - gtoPhase0, 1) + n)
%       orientDeg = orient0 + sgn*(180/pi)*tDep
%       tArr      = tDep + tf(iD, iA, iS)
%       resid     = wrap(sA(iA) - mod(delta + (j-1)/6 + tArr/tauD, 1))
%   admissible iff abs(wrap180(orientDeg - sheetOrient(iS))) <= orientBound
%   and abs(resid) <= tolPhase (wrap into (-0.5, 0.5]; a NaN cell is never
%   admissible). The departure-phase difference is snapped to 0 when it
%   lies within 1e-9 below 1 (floating noise would cost a whole period).
%
%   PER TRIPLE the six slots are independent: each takes its earliest
%   admissible candidate not in opts.exclude (sort key [tArr to 1e-12, tf,
%   cell order]); two stages may share a departure instant. The triple is
%   infeasible if any slot has none. lastArr = max(tArr) and T0 is the
%   smallest t >= lastArr with mod(delta + t/tauD, 1/6) == 0 (lastArr
%   itself when it satisfies this to 1e-12). In units of 1/6 DRO period:
%   a lastArr within 1e-12 below OR above a grid instant is treated as ON
%   it; one clearly above (e.g. 1e-9) jumps to the next instant.
%
%   OBJECTIVE over triples: minimise T0; ties within 1e-12 go to less total
%   propellant (propKgPerNd * sum of the six chosen tf), then (ties within
%   1e-9 deg) to the smaller maxOrientOffsetDeg = max over the six chosen
%   departures of abs(wrap180(orientDeg - sheetOrient(iS))) (the centred
%   triple), then to the lexicographically smaller [orient0 gtoPhase0
%   delta] (grids are sorted and the search is in that order).
%
%   LOCKED TRIPLE (opts.lockTriple): only that triple is evaluated and the
%   three grids are ignored. Re-scheduling after an exclude then leaves
%   every OTHER satellite's candidate (iD, iA, iS, n, tDep, tArr)
%   unchanged. T0 and the leg-2 numbering k MAY change when the replacement
%   arrives later than the old lastArr, so the caller must re-read S.T0
%   and every k. If the triple is infeasible under exclude, throws
%   schedule_gto_dro:noAssignment.
%
%   RELABEL: output index k is the leg-2 spacecraft number. Slot j has DRO
%   phase mod(delta + (j-1)/6 + T0/tauD, 1) at T0, a multiple of 1/6; its
%   satellite is labelled k = 6*phase + 1 (mod 6), so the satellite at
%   phase 0 at T0 is leg-2 spacecraft 1.
%
%   S.cands{k} lists EVERY admissible candidate of that satellite's slot
%   under the chosen triple (including excluded ones, flagged .excluded),
%   sorted by tArr; the chosen one is the first with excluded == false.
%
%   Throws schedule_gto_dro:noAssignment when no triple fills all six
%   slots (message: the triple with most slots filled and its six per-slot
%   candidate counts) and schedule_gto_dro:exclude for a malformed exclude
%   list (not numeric, not m x 4, non-integer, or outside the table / wait
%   range 0:nMaxRevs).
%
%   Self-demo (nargin == 0): schedules the real leg-1 table.
%
%% Inputs:
%
%  tf          nD x nA x nS Catalog min-time table [ND], NaN = unusable
%  sD          1 x nD       Catalog departure (GTO) phases [fraction]
%  sA          1 x nA       Catalog arrival (DRO) phases [fraction]
%  sheetOrient 1 x nS       Sheet GTO orientations [deg]
%  tauG        scalar       GTO Kepler period [ND]
%  tauD        scalar       DRO period [ND]
%  sgn         scalar       +1 or -1: orientation drift sign (leg1_inputs)
%  opts        struct, optional
%    .orient0Grid     1 x m  [-180:5:175] deg
%    .gtoPhase0Grid   1 x m  [(0:11)/12] fraction
%    .deltaGrid       1 x m  [(0:11)/72] fraction
%    .nMaxRevs        scalar [12] largest wait, whole GTO periods
%    .tolPhase        scalar [0.02] arrival-phase tolerance [fraction]
%    .orientBound     scalar [25] sheet orientation bound [deg]
%    .propKgPerNd     scalar [1] propellant per ND of burn [kg/ND]
%    .exclude         m x 4  Rows [iD iA iS n] barred for every slot
%                            [zeros(0, 4)]
%    .lockTriple      struct .orient0 .gtoPhase0 .delta: evaluate only
%                            this triple (grids ignored) [absent]
%
%% Outputs:
%
%  S.orient0 S.gtoPhase0 S.delta  scalars  The chosen triple
%  S.maxOrientOffsetDeg scalar  max over the six chosen departures of
%                        abs(wrap180(orientDeg - sheetOrient(iS))) [deg]
%  S.T0         scalar   Leg-2 time origin on the leg-1 clock [ND]
%  S.lastArr    scalar   Last arrival [ND]
%  S.iD S.iA S.iS S.n S.tDep S.tArr S.tf S.orientDeg S.slotPhaseAtDep
%  S.resid      1 x 6    Chosen candidate per satellite k (leg-2 numbering)
%                        (slotPhaseAtDep: the slot's DRO phase at tDep)
%  S.slot       1 x 6    Slot index j of satellite k
%  S.propKg     1 x 6    propKgPerNd*tf per satellite [kg]
%  S.totalPropKg scalar  Sum of S.propKg [kg]
%  S.cands      {1 x 6}  Struct arrays (fields iD iA iS n tDep tArr tf
%                        orientDeg slotPhaseAtDep resid excluded), by tArr
%  S.nCand      1 x 6    numel of each S.cands{k}
%  S.nFeasible  scalar   Number of feasible triples
%  S.opts       struct   Options as used
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  centred tie-break, lockTriple, exclude validation (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    L = leg1_inputs(struct());
    o = struct('propKgPerNd', L.m0kg*L.phys.Tnd/L.phys.cnd);
    t0 = tic;
    S = schedule_gto_dro(L.tf, L.sD, L.sA, L.sheetOrient, L.tauG, L.tauD, L.sgn, o);
    wall = toc(t0);
    fprintf('triple: orient0 %g deg  gtoPhase0 %g  delta %g   T0 %.3f d  wall %.2f s\n', ...
            S.orient0, S.gtoPhase0, S.delta, S.T0*L.tStar/86400, wall);
    fprintf('%-3s %-5s %-3s %-3s %-3s %-3s %9s %9s %9s %8s %7s\n', 'k', 'slot', ...
            'iD', 'iA', 'iS', 'n', 'tDep[d]', 'tArr[d]', 'orient', 'resid', 'nCand');
    for kk = 1:6
        fprintf('%-3d %-5d %-3d %-3d %-3d %-3d %9.3f %9.3f %9.2f %8.4f %7d\n', kk, ...
                S.slot(kk), S.iD(kk), S.iA(kk), S.iS(kk), S.n(kk), ...
                S.tDep(kk)*L.tStar/86400, S.tArr(kk)*L.tStar/86400, ...
                S.orientDeg(kk), S.resid(kk), S.nCand(kk));
    end
    return
end

if nargin < 8, opts = struct(); end
dflt = struct('orient0Grid', -180:5:175, 'gtoPhase0Grid', (0:11)/12, ...
              'deltaGrid', (0:11)/72, 'nMaxRevs', 12, 'tolPhase', 0.02, ...
              'orientBound', 25, 'propKgPerNd', 1, 'exclude', zeros(0, 4));
fn = fieldnames(dflt);
for kf = 1:numel(fn)
    if ~isfield(opts, fn{kf}), opts.(fn{kf}) = dflt.(fn{kf}); end
end
ex = opts.exclude;
if isempty(ex), ex = zeros(0, 4); end
if ~isnumeric(ex) || ~ismatrix(ex) || size(ex, 2) ~= 4
    error('schedule_gto_dro:exclude', ...
          'opts.exclude must be a numeric m x 4 matrix [iD iA iS n] (got %s)', mat2str(size(ex)));
end

if isfield(opts, 'lockTriple') && ~isempty(opts.lockTriple)
    lk = opts.lockTriple;
    opts.orient0Grid = lk.orient0;  opts.gtoPhase0Grid = lk.gtoPhase0;
    opts.deltaGrid = lk.delta;
end
orientGrid = sort(unique(opts.orient0Grid(:))).';
gtoGrid    = sort(unique(opts.gtoPhase0Grid(:))).';
deltaGrid  = sort(unique(opts.deltaGrid(:))).';
sD = sD(:);  sA = sA(:).';  sheetOrient = sheetOrient(:);
[nD, nA, nS] = size(tf);
assert(numel(sD) == nD && numel(sA) == nA && numel(sheetOrient) == nS, ...
       'schedule_gto_dro:size', 'table size');
nW   = opts.nMaxRevs + 1;
nn   = reshape(0:opts.nMaxRevs, 1, 1, 1, nW);
nCell = nD*nA*nS;
tfv  = tf(:);
tol  = opts.tolPhase;
slotOff = (0:5)/6;
wrap180 = @(x) mod(x + 180, 360) - 180;

% barred candidates, [nD nA nS nW]
dimsEx = [nD nA nS opts.nMaxRevs + 1];
lo = [1 1 1 1];
excl = false(nD, nA, nS, nW);
for kr = 1:size(ex, 1)
    row = ex(kr, :) + [0 0 0 1];                       % wait n -> index n + 1
    if any(~isfinite(row)) || any(row ~= round(row)) || any(row < lo) || any(row > dimsEx)
        error('schedule_gto_dro:exclude', ...
              'opts.exclude row %d [%s] is not an integer [iD iA iS n] inside the table (%d x %d x %d, n in 0:%d)', ...
              kr, num2str(ex(kr, :)), nD, nA, nS, opts.nMaxRevs);
    end
    excl(row(1), row(2), row(3), row(4)) = true;
end

% per gtoPhase0: departure times and arrival times
nG = numel(gtoGrid);  nDl = numel(deltaGrid);  nO = numel(orientGrid);
tDepG = cell(1, nG);  tArrG = cell(1, nG);
for kg = 1:nG
    dph = mod(sD - gtoGrid(kg), 1);
    dph(dph > 1 - 1e-9) = 0;
    tDepG{kg} = tauG*(dph + nn);                       % nD x 1 x 1 x nW
    tArrG{kg} = tDepG{kg} + tf;                        % nD x nA x nS x nW
end

% per (gtoPhase0, delta, slot): phase residual and admissibility
resG = cell(nG, nDl, 6);  okG = cell(nG, nDl, 6);
for kg = 1:nG
    for kd = 1:nDl
        for js = 1:6
            ph = mod(deltaGrid(kd) + slotOff(js) + tArrG{kg}/tauD, 1);
            r  = mod(sA - ph, 1);
            r(r > 0.5) = r(r > 0.5) - 1;
            resG{kg, kd, js} = r;
            okG{kg, kd, js}  = abs(r) <= tol;
        end
    end
end

% search over triples, lexicographic order (strict improvement only)
bestT0 = inf;  bestProp = inf;  bestOff = inf;  bestTrip = [0 0 0];  nFeas = 0;
bestFilled = -1;  bestCounts = zeros(1, 6);  bestFilledTrip = [0 0 0];
for ko = 1:nO
    for kg = 1:nG
        orientDeg = orientGrid(ko) + sgn*(180/pi)*tDepG{kg};          % nD x 1 x 1 x nW
        sheetOK = abs(wrap180(orientDeg - reshape(sheetOrient, 1, 1, nS))) ...
                  <= opts.orientBound;                                % nD x 1 x nS x nW
        avail = sheetOK & ~excl;                                      % nD x nA x nS x nW
        for kd = 1:nDl
            counts = zeros(1, 6);  tPick = nan(1, 6);  fPick = nan(1, 6);  oPick = nan(1, 6);
            for js = 1:6
                idx = find(avail & okG{kg, kd, js});
                counts(js) = numel(idx);
                if ~isempty(idx)
                    cellIdx = mod(idx - 1, nCell) + 1;
                    key = [round(tArrG{kg}(idx)*1e12), tfv(cellIdx), idx];
                    [~, ord] = min_row(key);
                    tPick(js) = tArrG{kg}(idx(ord));
                    fPick(js) = tfv(cellIdx(ord));
                    [idP, ~, isP, iwP] = ind2sub([nD nA nS nW], idx(ord));
                    oPick(js) = abs(wrap180(orientGrid(ko) + sgn*(180/pi)*tDepG{kg}(idP, 1, 1, iwP) ...
                                            - sheetOrient(isP)));
                end
            end
            nFilled = nnz(counts > 0);
            if nFilled > bestFilled
                bestFilled = nFilled;  bestCounts = counts;
                bestFilledTrip = [orientGrid(ko) gtoGrid(kg) deltaGrid(kd)];
            end
            if nFilled < 6, continue, end
            nFeas = nFeas + 1;
            t0 = t0_of(max(tPick), deltaGrid(kd), tauD);
            prop = opts.propKgPerNd*sum(fPick);
            offs = max(oPick);
            if t0 < bestT0 - 1e-12 || (abs(t0 - bestT0) <= 1e-12 && ...
               (prop < bestProp - 1e-12 || (abs(prop - bestProp) <= 1e-12 && offs < bestOff - 1e-9)))
                bestT0 = t0;  bestProp = prop;  bestOff = offs;
                bestTrip = [orientGrid(ko) gtoGrid(kg) deltaGrid(kd)];
            end
        end
    end
end

if nFeas == 0
    error('schedule_gto_dro:noAssignment', ...
          ['no triple fills all six slots; best triple (orient0 %g, gtoPhase0 %g, ' ...
           'delta %g) fills %d of 6; per-slot candidate counts [%s]'], ...
          bestFilledTrip(1), bestFilledTrip(2), bestFilledTrip(3), bestFilled, ...
          strtrim(sprintf('%d ', bestCounts)));
end

% rebuild the chosen triple in full
ko = find(orientGrid == bestTrip(1), 1);
kg = find(gtoGrid == bestTrip(2), 1);
kd = find(deltaGrid == bestTrip(3), 1);
delta = bestTrip(3);
orientDeg = orientGrid(ko) + sgn*(180/pi)*tDepG{kg};
sheetOK = abs(wrap180(orientDeg - reshape(sheetOrient, 1, 1, nS))) <= opts.orientBound;
candsSlot = cell(1, 6);  pickSlot = zeros(1, 6);
for js = 1:6
    idx = find(sheetOK & okG{kg, kd, js});
    cellIdx = mod(idx - 1, nCell) + 1;
    key = [round(tArrG{kg}(idx)*1e12), tfv(cellIdx), idx];
    [~, ord] = sortrows(key);
    idx = idx(ord);
    [id, ia, is, iw] = ind2sub([nD nA nS nW], idx);
    tdv = reshape(tDepG{kg}(sub2ind([nD 1 1 nW], id, ones(size(id)), ones(size(id)), iw)), [], 1);
    cs = struct('iD', num2cell(id), 'iA', num2cell(ia), 'iS', num2cell(is), ...
                'n', num2cell(iw - 1), 'tDep', num2cell(tdv), ...
                'tArr', num2cell(tArrG{kg}(idx)), 'tf', num2cell(tfv(mod(idx - 1, nCell) + 1)), ...
                'orientDeg', num2cell(orientGrid(ko) + sgn*(180/pi)*tdv), ...
                'slotPhaseAtDep', num2cell(mod(delta + slotOff(js) + tdv/tauD, 1)), ...
                'resid', num2cell(resG{kg, kd, js}(idx)), ...
                'excluded', num2cell(excl(idx)));
    candsSlot{js} = reshape(cs, 1, []);
    pickSlot(js) = find(~[cs.excluded], 1);
end
lastArr = -inf;
for js = 1:6
    lastArr = max(lastArr, candsSlot{js}(pickSlot(js)).tArr);
end
T0 = t0_of(lastArr, delta, tauD);

S.orient0 = bestTrip(1);  S.gtoPhase0 = bestTrip(2);  S.delta = delta;
S.T0 = T0;  S.lastArr = lastArr;
fl = {'iD', 'iA', 'iS', 'n', 'tDep', 'tArr', 'tf', 'orientDeg', 'slotPhaseAtDep', 'resid'};
for kf = 1:numel(fl), S.(fl{kf}) = zeros(1, 6); end
S.slot = zeros(1, 6);
S.cands = cell(1, 6);
for js = 1:6
    ph = mod(delta + slotOff(js) + T0/tauD, 1);
    kk = mod(round(6*ph), 6) + 1;
    c  = candsSlot{js}(pickSlot(js));
    for kf = 1:numel(fl), S.(fl{kf})(kk) = c.(fl{kf}); end
    S.slot(kk)  = js;
    S.cands{kk} = candsSlot{js};
end
S.maxOrientOffsetDeg = max(abs(wrap180(S.orientDeg - reshape(sheetOrient(S.iS), 1, []))));
S.propKg = opts.propKgPerNd*S.tf;
S.totalPropKg = sum(S.propKg);
S.nCand = cellfun(@numel, S.cands);
S.nFeasible = nFeas;
S.opts = opts;
end

function t0 = t0_of(lastArr, delta, tauD)
% Smallest t >= lastArr with mod(delta + t/tauD, 1/6) == 0, taking lastArr
% itself when it satisfies this to 1e-12.
kq = ceil(6*(delta + lastArr/tauD) - 1e-12);
t0 = tauD*(kq/6 - delta);
end

function [m, ord] = min_row(key)
% Lexicographically smallest row of key (ties keep the first).
ord = 1;
for kr = 2:size(key, 1)
    a = key(kr, :);  b = key(ord, :);
    d = find(a ~= b, 1);
    if ~isempty(d) && a(d) < b(d), ord = kr; end
end
m = key(ord, :);
end
