function G = mission_gdop(H, M, window, grid, opts)
%% Purpose:
%
%   Lunar south-pole GDOP of the six-spacecraft mission through proj7's
%   engine gdop_grid_from_states, and the paper's requirements R1-R3 by
%   the recipe of run_requirements_table, on the proj7 scorers
%   max_min_access_bounds / max_min_access_block_bounds. Two call forms:
%
%     G = mission_gdop(H, M, window, grid, opts)   mission window
%     V = mission_gdop('validate', M, opts)        validation gate
%
%   EPOCHS. tq = t0 + (0:Nt-1)*dt, dt = grid.dtSec/tStar, Nt =
%   round((t1 - t0)/dt) (opts.Nt overrides: dt = (t1 - t0)/Nt; the
%   validation uses Nt = 5760 over one tulip period exactly as the
%   reference did).
%
%   ON-STATION SET. Spacecraft k is a navigation source at tq iff tq >=
%   H.tArr(k) (always, for NaN tArr). It is decided from tArr, not by
%   interpolating the int8 phase; H.phase is checked to agree with tArr at
%   every H.t sample (mission_gdop:phase otherwise).
%
%   POSITIONS. A spacecraft on station sits on its slot, so its position
%   at any epoch is evaluated EXACTLY from the slot definition the history
%   used: M.arrState(mod(M.slotPhase0(H.assign(k)) + t/tauT, 1)). Before
%   use, H.r is checked to equal that slot state at every on-station H.t
%   sample (1e-9 ND, else mission_gdop:slot), so H and M cannot disagree.
%   (Interpolating H.r instead would need samples before tArr, i.e. the
%   transfer, inside the spline stencil at a run's start.)
%
%   STITCH. tq is split into maximal runs over which the on-station set is
%   constant. Each run with M_on >= 1 is cut into time blocks of at most
%   floor(opts.maxBaseBytes/(24*M_on*Npts)) epochs (the engine's dop
%   builds [n x 3 x M_on x Npts] doubles plus several same-size
%   temporaries); each block calls the engine on tau = [tq(blk); tq(blk(end))
%   + dt] with time_samples = numel(blk), so its half-open uniform
%   resampling lands on exactly the block's epochs (every one a spline
%   node; asserted on R.tau_i to 1e-12 ND). One- and two-epoch blocks are
%   fine: interp1 'spline' with 2 or 3 nodes reproduces its nodes. The
%   engine is pointwise in time once positions are given, so blocks are
%   exact. Runs with no spacecraft on station get NaN GDOP and zero
%   visibility without a call. Fewer than four in view gives NaN GDOP
%   from the engine itself.
%
%   SCORES per grid point (verbatim recipe): rak = struct('dt',
%   median(diff(tSec)), 'bin_duration', 86400); R1 =
%   max_min_access_bounds(vis >= 2, rak); R2 = max_min_access_bounds(gdop
%   < 6, rak); R3 = max_min_access_block_bounds(tSec, gdop < 6,
%   struct('min_block_duration', 4*3600, 'total_duration', 9*3600))/3600;
%   R3g7 likewise on gdop < 7. Floors = minimum over grid points; pass at
%   0.90 / 0.40 / 9 h (R3 fail count uses 9 - 1e-9 as the script does).
%   A scored window shorter than one 86400 s bin is refused
%   (mission_gdop:window; the scorers would return -inf).
%
%   VALIDATION ('validate'). On-slots history over one tulip period, the
%   paper grid (or opts.chunkSel) at Nt = 5760, scored, then compared to
%   S2(1).r of proj7/figures/requirements_table_results.mat on the SAME
%   grid points (matched by lat/lon): status PASS iff each of the three
%   floors over the compared points is within tol (R1 0.005, R2 0.005, R3
%   0.10 h) of the reference floor over those points; FAIL otherwise.
%   UNRESOLVED when the file is missing, has no S2, or S2(1).r cannot be
%   matched to this study's inputs: name '6-sat tulip', Nt 5760, mask 0,
%   10,770 points, dt_min equal to tauT*tStar/5760/60 (1e-6 min); and
%   the REFERENCE IDENTITY: the file stores no phasing vector, so its
%   phasing is pinned here as the literal it was computed with (2026-08-13,
%   global-constraints.md) -- the one place the phasing is typed.
%   tulip_phasing() AND M.dtau6 must both equal it to 1e-12, else
%   UNRESOLVED 'reference predates current phasing' (an edited
%   tulip_phasing can then never be compared with the old file).
%   opts.refPhasing replaces the literal (test seam only).
%   ORDER: every identity check above runs on the UNMODIFIED inputs first;
%   only then does opts.dtauOverride move the slots of the VALIDATED
%   history (deliberate-failure seam: an experiment on OUR constellation
%   with the reference identity unchanged), so an override yields FAIL
%   (or PASS), never UNRESOLVED.
%
%% Inputs:
%
%  H                        struct | 'validate'     assemble_histories output
%                                                   (.t [N x 1] ND, .r [N x 3 x
%                                                   6], .phase [N x 6] int8,
%                                                   .tArr [1 x 6], .assign
%                                                   [1 x 6], .tEnd)
%  M                        struct                  mission_inputs() output
%  window                   [1 x 2]                 [t0 t1] clock times [ND]
%                           (validate form: opts)
%  grid                     struct                  mission_grid() output
%  opts                     struct (optional)
%    .score                 logical [true]          Score R1-R3
%    .chunkSel              [1 x c] [all]           Chunk indices
%    .keepSeries            logical [false]         Keep .gdop and .vis
%    .logFile               char ['']               Append per-chunk log
%    .Nt                    scalar []               Epoch count override
%    .maxBaseBytes          scalar [4e8]            Engine base array cap
%                                                   per call [bytes]; the
%                                                   measured process peak is
%                                                   about 7x this (R2026a,
%                                                   paper chunk, Nt 5760:
%                                                   3e8 -> +2.0 GB RSS)
%    .quiet                 logical [false]         No console log
%    .r3g7                  logical [true]          Also score R3 on
%                                                   GDOP < 7 (NaN when
%                                                   false, worst.r3g7
%                                                   [NaN NaN]; the
%                                                   validation skips it)
%    validate form only:
%    .refFile               char                    Reference .mat
%    .dtauOverride          [1 x 6] []              Slot offsets [ND]
%    .refPhasing            [1 x 6] [literal]       Reference identity
%                                                   (test seam)
%
%% Outputs:
%
%  G.t [Nt x 1] ND, G.tSec [Nt x 1] s, G.nOn [Nt x 1]
%  G.lat, G.lon             [1 x Npts]              Grid points [deg]
%  G.poleGdop               [Nt x 1]                GDOP at the lowest
%                                                   latitude, longitude 0
%  G.poleLatLon             [1 x 2]                 That point [deg]
%  G.r1 G.r2 G.r3 G.r3g7    [1 x Npts]              Per-point scores
%                                                   (fraction, fraction, h, h)
%  G.r1Floor G.r2Floor G.r3Floor G.r3g7Floor        Floors
%  G.r1Nfail G.r2Nfail G.r3Nfail                    Points below threshold
%  G.pass                   struct                  .r1 .r2 .r3 logical
%  G.worst                  struct                  .r1 .r2 .r3 .r3g7 =
%                                                   [lat lon] deg of the floor
%  G.runs                   [nRuns x 3]             first, last epoch, M_on
%  G.grid, G.window, G.wallS, G.nCalls
%  G.timing                 [C x 4]                 per chunk: engine, R1+R2,
%                                                   R3, R3g7 wall [s]
%  G.gdop [Nt x Npts], G.vis [Nt x Npts] uint8      when keepSeries
%
%  V.status ('PASS'|'FAIL'|'UNRESOLVED'), V.reason, V.ours / V.ref
%  (.r1Floor .r2Floor .r3Floor; ref over the compared points), V.refFile
%  floors in V.refFull, V.tol, V.perPointMaxDiff (.r1 .r2 .r3 .nPts, .nDiff
%  = points differing [R1 R2 R3]),
%  V.G (the scored result without series), V.wallS
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo: coarse grid, on-slots, two days
if nargin == 0
    M0 = mission_inputs();
    H0 = assemble_histories(M0, [], 1:6, struct('onSlots', true, ...
                            'tEnd', 2*86400/M0.tStar));
    G0 = mission_gdop(H0, M0, [0 2*86400/M0.tStar], mission_grid('coarse'));
    fprintf('coarse, 2 d: R1 %.4f  R2 %.4f  R3 %.3f h  (%.1f s)\n', ...
            G0.r1Floor, G0.r2Floor, G0.r3Floor, G0.wallS);
    return;
end

if ischar(H) || isstring(H)
    if ~strcmpi(char(H), 'validate')
        error('mission_gdop:mode', 'unknown mode "%s"', char(H));
    end
    if nargin < 3 || isempty(window), window = struct(); end
    G = validate_mode(M, window);
    return;
end
if nargin < 5, opts = struct(); end
G = core(H, M, window, grid, opts);
end

%% ======================================================================
function G = core(H, M, window, grid, opts)
wall0     = tic;
doScore   = getopt(opts, 'score', true);
chunkSel  = getopt(opts, 'chunkSel', 1:numel(grid.chunks));
keepSer   = getopt(opts, 'keepSeries', false);
logFile   = getopt(opts, 'logFile', '');
maxBytes  = getopt(opts, 'maxBaseBytes', 4e8);
quiet     = getopt(opts, 'quiet', false);
doG7      = getopt(opts, 'r3g7', true);
nSc       = size(H.r, 3);

t0 = window(1);  t1 = window(2);
if isfield(opts, 'Nt') && ~isempty(opts.Nt)
    Nt = opts.Nt;
    dt = (t1 - t0)/Nt;
else
    dt = grid.dtSec/M.tStar;
    Nt = round((t1 - t0)/dt);
end
if Nt < 1 || t0 < 0 || t1 > H.tEnd + 1e-9*dt
    error('mission_gdop:window', 'window [%.6g %.6g] ND outside [0 %.6g] or empty', ...
          t0, t1, H.tEnd);
end
if doScore && Nt < round(86400/(dt*M.tStar))
    error('mission_gdop:window', ['a scored window needs at least one full ' ...
          '86400 s bin: %d epochs of %.1f s'], Nt, dt*M.tStar);
end
tq = t0 + (0:Nt-1).'*dt;

%% On-station set from tArr; H.phase and H.r must agree with it
tArr = H.tArr(:).';
onAt = @(tt) isnan(tArr) | tt >= tArr;               % [n x nSc]
onH  = onAt(H.t(:));
if ~isequal(onH, H.phase == 3)
    [ii, kk] = find(onH ~= (H.phase == 3), 1);
    error('mission_gdop:phase', ['H.phase disagrees with H.tArr: spacecraft %d ' ...
          'at H.t(%d) = %.6g ND (phase %d)'], kk, ii, H.t(ii), H.phase(ii, kk));
end
slotPos = @(tt, kk) slot_positions(M, H.assign(kk), tt);
for kk = 1:nSc
    rows = find(onH(:, kk));
    if isempty(rows), continue; end
    d = max(abs(H.r(rows, :, kk) - slotPos(H.t(rows), kk)), [], 'all');
    if d > 1e-9
        error('mission_gdop:slot', ['H.r of spacecraft %d differs from its slot ' ...
              '%d state by %.3g ND on station'], kk, H.assign(kk), d);
    end
end
onSet = onAt(tq);                                    % [Nt x nSc]
nOn   = sum(onSet, 2);

%% Runs of constant on-station set
key   = onSet*(2.^(0:nSc-1)).';
edges = [0; find(diff(key) ~= 0); Nt];
runs  = [edges(1:end-1) + 1, edges(2:end), nOn(edges(2:end))];

%% Chunks
nCh    = numel(chunkSel);
latC   = cell(1, nCh);  lonC = latC;  gdC = latC;  visC = latC;
scC    = cell(1, nCh);
timing = zeros(nCh, 4);
nCalls = 0;
% pole series: the lowest latitude over the selected chunks, longitude 0
% (or the smallest longitude evaluated)
chMin  = cellfun(@(c) c(1), grid.chunks(chunkSel));
[~, poleIc] = min(chMin);
poleGd = [];  poleLL = [NaN NaN];
tSec   = tq*M.tStar;
for ic = 1:nCh
    cc = chunkSel(ic);
    p  = struct('lat_range', grid.chunks{cc}, 'lon_range', grid.lon_range, ...
                'lat_res', grid.lat_res, 'lon_res', grid.lon_res, ...
                'elevation_mask', grid.mask);
    latV = p.lat_range(1):p.lat_res:p.lat_range(2);
    lonV = p.lon_range(1):p.lon_res:p.lon_range(2);
    [La, Lo] = meshgrid(latV, lonV);                 % engine point order
    nPts = numel(La);
    gd   = nan(Nt, nPts);
    vis  = zeros(Nt, nPts, 'uint8');
    tc   = tic;
    for ir = 1:size(runs, 1)
        mOn = runs(ir, 3);
        if mOn == 0, continue; end
        sats = find(onSet(runs(ir, 1), :));
        blk  = max(1, floor(maxBytes/(24*mOn*nPts)));
        for b0 = runs(ir, 1):blk:runs(ir, 2)
            idx = (b0:min(b0 + blk - 1, runs(ir, 2))).';
            tau = [tq(idx); tq(idx(end)) + dt];
            rT  = zeros(numel(tau), 3, mOn);
            for mm = 1:mOn
                rT(:, :, mm) = slotPos(tau, sats(mm));
            end
            p.time_samples = numel(idx);
            R = gdop_grid_from_states(tau, rT, M.muStar, M.lStar, M.tStar, p);
            nCalls = nCalls + 1;
            if max(abs(R.tau_i - tq(idx))) > 1e-12 || ...
               max(abs(R.Lat_grid(:) - La(:)*pi/180)) > 1e-12
                error('mission_gdop:engineGrid', ...
                      'engine resampling or grid does not match the requested epochs');
            end
            gd(idx, :)  = R.gdop_time_series;
            vis(idx, :) = uint8(R.visibility_count);
        end
    end
    timing(ic, 1) = toc(tc);
    latC{ic} = La(:).';  lonC{ic} = Lo(:).';
    if doScore
        [scC{ic}, timing(ic, 2:4)] = score_points(gd, vis, tSec, doG7);
    end
    if ic == poleIc
        cand = find(La(:) == min(La(:)));
        [~, kmin] = min(Lo(cand));
        poleGd = gd(:, cand(kmin));
        poleLL = [La(cand(kmin)), Lo(cand(kmin))];
    end
    if keepSer, gdC{ic} = gd;  visC{ic} = vis; end
    msg = sprintf(['mission_gdop: chunk %d (lat %.1f..%.1f, %d pts, Nt %d): ' ...
                   'engine %.1f s, R1+R2 %.1f s, R3 %.1f s, R3g7 %.1f s'], cc, p.lat_range, ...
                   nPts, Nt, timing(ic, :));
    if ~quiet, fprintf('%s\n', msg); end
    if ~isempty(logFile)
        fid = fopen(logFile, 'a');
        if fid > 0, fprintf(fid, '%s %s\n', char(datetime('now')), msg); fclose(fid); end
    end
end

%% Assemble
G.t      = tq;
G.tSec   = tSec;
G.nOn    = nOn;
G.lat    = [latC{:}];
G.lon    = [lonC{:}];
G.runs   = runs;
G.grid   = grid;
G.window = [t0 t1];
G.poleGdop   = poleGd;
G.poleLatLon = poleLL;
if keepSer
    G.gdop = [gdC{:}];
    G.vis  = [visC{:}];
end
if doScore
    S = [scC{:}];
    G.r1   = [S.r1];   G.r2 = [S.r2];   G.r3 = [S.r3];   G.r3g7 = [S.r3g7];
    [G.r1Floor, w1] = min(G.r1);
    [G.r2Floor, w2] = min(G.r2);
    [G.r3Floor, w3] = min(G.r3);
    [G.r3g7Floor, w4] = min(G.r3g7);
    G.r1Nfail = sum(G.r1 < 0.90);
    G.r2Nfail = sum(G.r2 < 0.40);
    G.r3Nfail = sum(G.r3 < 9 - 1e-9);
    G.pass  = struct('r1', G.r1Floor >= 0.90, 'r2', G.r2Floor >= 0.40, ...
                     'r3', G.r3Floor >= 9 - 1e-9);
    G.worst = struct('r1', [G.lat(w1) G.lon(w1)], 'r2', [G.lat(w2) G.lon(w2)], ...
                     'r3', [G.lat(w3) G.lon(w3)], 'r3g7', [G.lat(w4) G.lon(w4)]);
    if ~doG7, G.worst.r3g7 = [NaN NaN]; end
end
G.timing = timing;
G.nCalls = nCalls;
G.wallS  = toc(wall0);
end

%% ======================================================================
function [S, tm] = score_points(gd, vis, tSec, doG7)
% run_requirements_table's verdict recipe, per grid point.
rak  = struct('dt', median(diff(tSec)), 'bin_duration', 86400);
req  = struct('min_block_duration', 4*3600, 'total_duration', 9*3600);
cov6 = gd < 6;
cov7 = gd < 7;
cR1  = vis >= 2;
nPts = size(gd, 2);
r1 = nan(1, nPts);  r2 = r1;  r3 = r1;  r3g7 = r1;
ta = tic;
for pp = 1:nPts
    r1(pp) = max_min_access_bounds(cR1(:, pp), rak);
    r2(pp) = max_min_access_bounds(cov6(:, pp), rak);
end
tm(1) = toc(ta);
tb = tic;
for pp = 1:nPts
    r3(pp) = max_min_access_block_bounds(tSec, cov6(:, pp), req)/3600;
end
tm(2) = toc(tb);
tc = tic;
if doG7
    for pp = 1:nPts
        r3g7(pp) = max_min_access_block_bounds(tSec, cov7(:, pp), req)/3600;
    end
end
tm(3) = toc(tc);
S = struct('r1', r1, 'r2', r2, 'r3', r3, 'r3g7', r3g7);
end

%% ======================================================================
function r = slot_positions(M, slot, tt)
% Exact slot position [n x 3] ND at clock times tt [n x 1].
rv = M.arrState(mod(M.slotPhase0(slot) + tt(:).'/M.tauT, 1));
r  = rv(1:3, :).';
end

%% ======================================================================
function V = validate_mode(M, opts)
wall0   = tic;
refFile = getopt(opts, 'refFile', fullfile(getenv('HOME'), 'Desktop', 'proj7', ...
                 'figures', 'requirements_table_results.mat'));
grid    = mission_grid('paper');
sel     = getopt(opts, 'chunkSel', 1:numel(grid.chunks));
NtRef   = 5760;
V.status  = 'UNRESOLVED';
V.reason  = '';
V.refFile = refFile;
V.tol     = struct('r1', 0.005, 'r2', 0.005, 'r3', 0.10);
nanF      = struct('r1Floor', NaN, 'r2Floor', NaN, 'r3Floor', NaN);
V.ours    = nanF;  V.ref = nanF;  V.refFull = nanF;

%% Reference: present, has S2, matches this study's inputs
if ~isfile(refFile)
    V.reason = sprintf('reference file not found: %s', refFile);
    return;
end
w = whos('-file', refFile);
if ~any(strcmp({w.name}, 'S2'))
    V.reason = sprintf('reference file has no S2: %s', refFile);
    return;
end
L = load(refFile, 'S2');
r = L.S2(1).r;
need = {'name', 'Nt', 'mask', 'dt_min', 'lat', 'lon', 'map_r1', ...
        'map_r2_bestday', 'map_r3s_pair', 'r1_floor', 'r2_floor', 'r3s_floor'};
miss = need(~isfield(r, need));
if ~isempty(miss)
    V.reason = sprintf('S2(1).r lacks %s', strjoin(miss, ', '));
    return;
end
V.refFull = struct('r1Floor', r.r1_floor, 'r2Floor', r.r2_floor, 'r3Floor', r.r3s_floor);
dtMinWant = M.tauT*M.tStar/NtRef/60;
dtau6Now  = tulip_phasing();
refPh     = getopt(opts, 'refPhasing', [0.778317343114462 1.14555587028854 ...
                   2.51741969625761 2.89299534327950 4.26677262887274 ...
                   4.63164623641942]);        % reference identity (see header)
why = {};
if ~strcmp(r.name, '6-sat tulip'), why{end+1} = sprintf('name "%s"', r.name); end
if r.Nt ~= NtRef,                  why{end+1} = sprintf('Nt %d', r.Nt);       end
if r.mask ~= 0,                    why{end+1} = sprintf('mask %g', r.mask);   end
if numel(r.lat) ~= 10770,          why{end+1} = sprintf('%d points', numel(r.lat)); end
if abs(r.dt_min - dtMinWant) > 1e-6
    why{end+1} = sprintf('dt_min %.8f vs %.8f', r.dt_min, dtMinWant);
end
if ~isempty(why)
    V.reason = ['reference cannot be matched: ' strjoin(why, '; ')];
    return;
end
if numel(dtau6Now) ~= 6 || numel(M.dtau6) ~= 6 || numel(refPh) ~= 6 || ...
   max(abs(dtau6Now(:).' - refPh(:).')) > 1e-12 || ...
   max(abs(M.dtau6(:).' - refPh(:).')) > 1e-12
    V.reason = sprintf(['reference predates current phasing: reference %s, ' ...
                        'tulip_phasing() %s, M.dtau6 %s'], mat2str(refPh(:).', 15), ...
                        mat2str(dtau6Now(:).', 15), mat2str(M.dtau6(:).', 15));
    return;
end

%% Our on-slots history and the scored paper grid
Mv = M;
if isfield(opts, 'dtauOverride') && ~isempty(opts.dtauOverride)
    Mv.dtau6      = opts.dtauOverride(:).';
    Mv.slotPhase0 = mod(Mv.dtau6, Mv.tauT)/Mv.tauT;
end
H = assemble_histories(Mv, [], 1:6, struct('onSlots', true, 'tEnd', Mv.tauT));
G = core(H, Mv, [0 Mv.tauT], grid, struct('chunkSel', sel, 'Nt', NtRef, ...
         'quiet', getopt(opts, 'quiet', false), 'r3g7', false, ...
         'logFile', getopt(opts, 'logFile', ''), ...
         'maxBaseBytes', getopt(opts, 'maxBaseBytes', 4e8)));

%% Same points in the reference
keyR = round([r.lat(:), r.lon(:)]*180/pi*1e6);
keyO = round([G.lat(:), G.lon(:)]*1e6);
[found, loc] = ismember(keyO, keyR, 'rows');
if ~all(found)
    V.reason = sprintf('%d of our grid points are not in the reference', nnz(~found));
    return;
end
ref1 = r.map_r1(loc).';  ref2 = r.map_r2_bestday(loc).';  ref3 = r.map_r3s_pair(loc).';
V.ours = struct('r1Floor', G.r1Floor, 'r2Floor', G.r2Floor, 'r3Floor', G.r3Floor);
V.ref  = struct('r1Floor', min(ref1), 'r2Floor', min(ref2), 'r3Floor', min(ref3));
V.perPointMaxDiff = struct('r1', max(abs(G.r1 - ref1)), 'r2', max(abs(G.r2 - ref2)), ...
                           'r3', max(abs(G.r3 - ref3)), 'nPts', numel(loc), ...
                           'nDiff', [nnz(G.r1 ~= ref1), nnz(abs(G.r2 - ref2) > 1e-12), ...
                                     nnz(abs(G.r3 - ref3) > 1e-9)]);
d = [abs(V.ours.r1Floor - V.ref.r1Floor), abs(V.ours.r2Floor - V.ref.r2Floor), ...
     abs(V.ours.r3Floor - V.ref.r3Floor)];
tl = [V.tol.r1, V.tol.r2, V.tol.r3];
V.floorDiff = d;
if all(d <= tl)
    V.status = 'PASS';
else
    V.status = 'FAIL';
end
pd = V.perPointMaxDiff;
V.reason = sprintf(['%d points, chunks %s: floor |diff| R1 %.6f R2 %.6f R3 %.4f h ' ...
                    '(tol %.3f %.3f %.2f h); per point %d R1 / %d R2 / %d R3 ' ...
                    'points differ, max R2 %.4f, max R3 %.2f h'], numel(loc), ...
                    mat2str(sel), d, tl, pd.nDiff, pd.r2, pd.r3);
V.G = G;
V.wallS = toc(wall0);
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
