function [outStem, info] = costate_catalog_extremes_movies(cat_, filt, metric, outStem, opts)
%% Purpose:
%
%   Renders the two EXTREME transfers found by costate_catalog_extremes as a
%   side-by-side MOVIE (.mp4 and .gif) in the pumpkynPie house style: the
%   lit dark CR3BP scene (textured Moon, Earth, colored terminal orbits,
%   colored trail, satellite marker) built by
%   pumpkynPie.plot.SatelliteAnimator -- the same scene as
%   dro_scene_movie.m and the pumpkynPie demos.
%
%   Left panel = the minimum, right = the maximum, separated by a vertical
%   divider, and BOTH RUN ON ONE CLOCK: frames advance on a shared time
%   base, so the shorter transfer arrives early and freezes on its final
%   frame while the longer one plays out. The clock (days) is shown in both
%   panels with identical values.
%
%  ASSUMPTIONS / NOTES:
%
% • SELF-CONTAINED on pumpkyn: the only helpers it needs are
%   costate_catalog_extremes (when no pair is given), status_key and
%   get_family_orbit. The trajectories are FLOWN with pumpkyn's own
%   propagator (pumpkyn.cr3bp.tfMinProp): from the stored multiple-shooting
%   junction states, segment by segment, when the row has them (the first
%   junction IS the departure state, so no endpoint interpolation enters);
%   otherwise in one flight from z8, started at the departure phase of the
%   departure orbit (spline on the orbit table). What you watch is the
%   actual PMP solution, an alternative's included.
% • Junctions sit on the uniform grid linspace(0, t_f, K+1): column k is the
%   augmented state [r; v; m; lam_r; lam_v; lam_m] at the start of segment k.
% • Orbits come from the sheet's own recipe (dep_family / dep_params,
%   arr_family / arr_params); a sheet without those keys (the deliverable-3
%   catalog) uses the DRO / tulip rule costate_catalog_extremes uses.
% • Aero-toolbox guard (ShowEarth trap): planetEphemeris without its data
%   pack is hidden from the path for this render only.
% • Each output frame is EXACTLY 1280 x 720 (two 637-px panels + a 6-px
%   divider), enforced by index resampling -- the known cure for the H.264
%   diagonal-shear artifact. The .gif is written at the same size.
% • Needs pumpkynPie.plot.SatelliteAnimator; if it is not on the path the
%   function warns and returns '' (same graceful degradation as
%   dro_scene_movie).
%
%% Inputs:
%
%  cat_                     struct                  A compact costate catalog
%                                                   (any family: dro_tulip,
%                                                   halo_tulip, ...)
%
%  filt                     struct                  Same filter as
%                                                   costate_catalog_extremes
%                                                   (.tauDRO .Np .thrustN,
%                                                   any subset; [] = all)
%
%  metric                   char                    'deltaV' (default) or
%                                                   'time'
%
%  outStem                  char                    Output stem for
%                                                   '<stem>.mp4' and
%                                                   '<stem>.gif'  [default
%                                                   'catalog_extremes_movie']
%
%  opts                     struct (optional)       .pair [] 2-element cell
%                                                   of transfers_at rows to
%                                                   film INSTEAD of the
%                                                   extremes (sheet 1, the
%                                                   library rung); .labels
%                                                   [] 2 chars (panel
%                                                   labels); .nFrames [240]
%                                                   frame cap; .fly ['auto']
%                                                   | 'junctions' | 'z8'
%
%% Outputs:
%
%  outStem                  char                    The stem written, or ''
%                                                   if the animator was
%                                                   unavailable
%
%  info                     struct                  .F {1x2} the two flights
%                                                   (.t .Y .tfNd .tfDays
%                                                   .dvKms .method .status
%                                                   .statusName) .methods
%                                                   {F1.method, F2.method}
%                                                   .nFrames
%
%% Revision History:
%  M. Casey                                                   (c) 08/06/2026
%  M. Casey  per-panel B (own sheet and thrust), full titles, guards (c) 10/13/2026
%  M. Casey  flights through fly_catalog_entry; opts.pair/.labels/
%            .nFrames/.fly; second output info                (c) 10/13/2026
%  M. Casey  self-demo uses the catalog's rung nearest 5 N (a
%            one-rung 70 mN catalog has no 5 N entries)   (c) 10/13/2026
%  M. Casey  ephemeris guard hides/restores only a folder under
%            matlabroot (ruling R11)                         (c) 10/13/2026
%  M. Casey  self-contained on pumpkyn: flights by tfMinProp from the
%            junctions (else z8), orbits by get_family_orbit at the
%            catalog's muStar; no flight-layer files (the user
%            kit, deliverable 9)                            (c) 10/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: phasing extremes on the mid sheet at 5 N. Runs on whichever
   %compact catalog sits in the current folder:
     F = dir('costate_catalog_*.mat');
     assert(~isempty(F), ...
            'demo: no costate_catalog_*.mat in the current folder');
       L = load(F(1).name);
      fn = fieldnames(L);
    cat_ = L.(fn{1});
   fprintf('demo catalog: %s\n', F(1).name);
      sh = cat_.sheets(ceil(numel(cat_.sheets)/2));
   [~,kDemo] = min(abs(cat_.rungs_N - 5));  tDemo = cat_.rungs_N(kDemo);  % 5 N when stored
   if isnan(sh.Np)
       fl = struct('tauDRO',sh.tauDRO,'tauArr',sh.tau_arr,'thrustN',tDemo);
   else
       fl = struct('tauDRO',sh.tauDRO,'Np',sh.Np,'thrustN',tDemo);
   end
 outStem = costate_catalog_extremes_movies(cat_, fl);
     return;
end
if ~exist('filt','var'),    filt = struct(); end
if ~exist('metric','var') || isempty(metric), metric = 'deltaV'; end
if ~exist('outStem','var'), outStem = 'catalog_extremes_movie'; end
if ~exist('opts','var') || isempty(opts), opts = struct(); end
if ~isfield(opts,'pair'),     opts.pair = [];      end
if ~isfield(opts,'labels'),   opts.labels = [];    end
if ~isfield(opts,'nFrames'),  opts.nFrames = [];   end
if ~isfield(opts,'fly') || isempty(opts.fly), opts.fly = 'auto'; end
assert(any(strcmp(opts.fly, {'auto', 'junctions', 'z8'})), ...
       'costate_catalog_extremes_movies:fly', 'opts.fly must be auto, junctions or z8');
info = struct('F', {{}}, 'methods', {{}}, 'nFrames', 0);

if isempty(which('pumpkynPie.plot.SatelliteAnimator'))
    warning('costate_catalog_extremes_movies:noAnimator', ...
        'pumpkynPie.plot.SatelliteAnimator not on the path -- no movie.');
    outStem = '';  return
end

%% Ephemeris guard: if planetEphemeris is on the path but its data package
%  is not installed, hide the toolbox folder for this render so pumpkyn's
%  analytic Earth-Moon model takes over (restored on exit; the Aerospace
%  Toolbox folder is not a repo folder):
if exist('planetEphemeris', 'file')
    try
        planetEphemeris(juliandate(datetime(2030,1,1)), 'Earth', 'Moon');
    catch
        aeroDir = fileparts(which('planetEphemeris'));
        if startsWith(aeroDir, matlabroot)   % only a MathWorks toolbox folder is ever hidden
            rmpath(aeroDir);
            cleanupPath = onCleanup(@() addpath(aeroDir));
            warning('costate_catalog_extremes_movies:noEphemData', ...
                'planetEphemeris found but its data package is missing; using pumpkyn''s analytic Earth-Moon model for this render.');
        end
    end
end

%% The two rows to film (the extremes, or opts.pair), their sheets and thrust:
   tStar = cat_.constants.tStar_s;
  params = struct('muStar',cat_.constants.muStar, 'lStar',cat_.constants.lStar_km, 'tStar',tStar);
   epoch = juliandate(datetime(2030,1,1,0,0,0));
if ~isempty(opts.pair)
    assert(iscell(opts.pair) && numel(opts.pair) == 2, ...
           'costate_catalog_extremes_movies:pair', 'opts.pair must be a 2-element cell of transfers_at rows');
    R = opts.pair;  ks = [1 1];  TN = cat_.rungs_N([1 1]);  ttl = {'', ''};
    if ~isempty(opts.labels)
        assert((iscell(opts.labels) || isstring(opts.labels) || ischar(opts.labels)) && numel(opts.labels) == 2, ...
               'costate_catalog_extremes_movies:labels', 'opts.labels must hold exactly 2 labels');
        lab = {char(opts.labels(1)), char(opts.labels(2))};
        if iscell(opts.labels), lab = {char(opts.labels{1}), char(opts.labels{2})}; end
    else,                     lab = {'A','B'}; end
else
    [eMin, eMax] = costate_catalog_extremes(cat_, filt, false, metric);
    if strncmpi(metric,'t',1), lab = {'FASTEST','SLOWEST'};
    else,                      lab = {'MIN \DeltaV','MAX \DeltaV'}; end
    R = cell(1,2);  ks = [0 0];  TN = [0 0];  ttl = cell(1,2);
    for k = 1:2
        e = {eMin, eMax};  e = e{k};
        R{k} = struct('z8',e.z8, 'sD',e.dep_frac, 'sA',e.arr_frac, 'junctions',e.junctions, ...
                      'status',e.status, 'statusName',e.statusName, 'reason',e.reason);
        ks(k) = find_sheet(cat_, e);  TN(k) = e.thrustN;
        ttl{k} = sprintf('%s \\tau=%.2f, %s, %.2f N', upper(e.dep_family), e.tauDRO, e.arrLab, e.thrustN);
    end
end
Fl = cell(1,2);  orb = cell(1,2);
for k = 1:2
    orb{k} = sheet_orbits(cat_.sheets(ks(k)), cat_.constants.muStar);
    Fl{k}  = fly_row(cat_, R{k}, orb{k}, TN(k), opts.fly);
end

for k = 1:2
    F = Fl{k};
    tauD = orb{k}.tD;  rvD = orb{k}.rvD;  tauT = orb{k}.tT;  rvT = orb{k}.rvT;
    [tu, iu] = unique(F.t);
    yj = F.Y;

    %% One pumpkyn-style lit scene per panel:
    S(k).fig = figure('Color','k', 'Position',[60 60 640 720], ...
                      'Visible','off', 'InvertHardcopy','off');
    S(k).anim = pumpkynPie.plot.SatelliteAnimator(S(k).fig, ...
        'Epoch',              epoch, ...
        'CR3BPParameters',    params, ...
        'ShowEarth',          true, ...
        'FigureTheme',        'dark', ...
        'BodyLighting',       'lambertian', ...
        'FitMode',            'trajectory', ...
        'TrailMode',          'distance', ...
        'EnableTrailColor',   true, ...
        'EnableSparkles',     false, ...
        'EnableSatelliteMarker', true, ...
        'CameraUpVector',     [0 0 1], ...
        'View',               [-35 22]);
    S(k).anim.addOrbit(epoch + tauD*tStar/86400, rvD, ...
        'Name','Departure orbit', 'Color',[0.25 0.60 0.35], 'LineWidth',0.9);
    S(k).anim.addOrbit(epoch + tauT*tStar/86400, rvT, ...
        'Name','Arrival orbit', 'Color',[0.75 0.28 0.28], 'LineWidth',0.9);
    S(k).anim.addSatellite(epoch + tu*tStar/86400, yj(iu,1:6), ...
        'Name', sprintf('%s transfer', lab{k}), ...
        'MarkerColor',[1.00 0.92 0.55], 'TrailColor',[0.35 0.72 1.00]);
    S(k).anim.initialize();
    camzoom(1.4);
    annotation(S(k).fig, 'textbox', [0.02 0.90 0.96 0.09], 'String', ...
        {sprintf('%s:  %s', lab{k}, strtrim([ttl{k} '  ' sprintf('[%s, %s flight]', F.statusName, F.method)])), ...
         sprintf('\\DeltaV %.3f km/s,  t_f %.3f d', F.dvKms, F.tfDays)}, ...
        'Color','w', 'EdgeColor','none', 'FontSize',11, ...
        'HorizontalAlignment','center');
    S(k).clock = annotation(S(k).fig, 'textbox', [0.02 0.02 0.96 0.05], ...
        'String','', 'Color',[1 0.92 0.55], 'EdgeColor','none', ...
        'FontSize',13, 'FontWeight','bold', 'HorizontalAlignment','center');
    S(k).tf = F.tfNd;
end
tEnd = max(S(1).tf, S(2).tf);            % the shared clock's last tick

%% Frames: shared clock; shorter panel freezes at its final pose:
 nFrames = 240;  if ~isempty(opts.nFrames), nFrames = opts.nFrames; end
     fps = 30;
      pw = 637;                          % panel width; 2*637 + 6 = 1280
 divider = uint8(90*ones(720, 1280-2*pw, 3));
vw = VideoWriter([outStem '.mp4'], 'MPEG-4');
vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
cleanupVW = onCleanup(@() close(vw));
gifFile = [outStem '.gif'];
for kf = 1:nFrames
    tNow = tEnd*(kf-1)/max(nFrames-1, 1);
    for k = 1:2
        S(k).anim.update(epoch + min(tNow, S(k).tf)*tStar/86400);
        set(S(k).clock, 'String', sprintf('t = %.2f days', tNow*tStar/86400));
    end
    drawnow limitrate
    img = [panel_frame(S(1).fig, pw), divider, panel_frame(S(2).fig, pw)];
    writeVideo(vw, img);
    if mod(kf,3) == 1                    % gif at 10 fps, same 1280x720 size
        [A, cmap] = rgb2ind(img, 256);
        if kf == 1
            imwrite(A, cmap, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, gifFile, 'gif', 'WriteMode','append', 'DelayTime', 0.1);
        end
    end
end
close(S(1).fig);  close(S(2).fig);
info.F = Fl;  info.methods = {Fl{1}.method, Fl{2}.method};  info.nFrames = nFrames;
fprintf('  movies -> %s.mp4 / .gif  (%d frames, 1280x720, shared clock to %.2f d)\n', ...
        outStem, nFrames, tEnd*tStar/86400);
end

% ---------------------------------------------------------------------------
function O = sheet_orbits(sh, muStar)
% SHEET_ORBITS  One period of the sheet's departure and arrival orbits, by
% the sheet's own recipe (dep_family / arr_family keys), else the DRO /
% tulip rule of the deliverable-3 catalog, at the catalog's own mass ratio.
% INPUTS: sh (sheet struct); muStar [scalar].  OUTPUTS: O (.tD .rvD .tT .rvT;
% [M x 1], [M x 6]).
if isfield(sh, 'dep_family'), df = sh.dep_family;  dp = sh.dep_params;
else,                         df = 'dro';          dp = struct('tau', sh.tauDRO); end
if isfield(sh, 'arr_family'), af = sh.arr_family;  ap = sh.arr_params;
else,                         af = 'tulip';        ap = struct('Np', sh.Np, 'pm', -1); end
dp.muStar = muStar;  ap.muStar = muStar;
[O.tD, O.rvD] = get_family_orbit(df, dp);
[O.tT, O.rvT] = get_family_orbit(af, ap);
end

% ---------------------------------------------------------------------------
function F = fly_row(cat_, row, O, thrustN, method)
% FLY_ROW  Fly one library row with pumpkyn.cr3bp.tfMinProp: segment by
% segment from its junction states (14 x K on the uniform grid) when it has
% them and the method allows, else in one flight from z8 started at the
% departure phase of the departure orbit.
% INPUTS: cat_ (catalog); row (z8, sD, junctions, status fields); O
% (sheet_orbits); thrustN [scalar]; method 'auto' | 'junctions' | 'z8'.
% OUTPUTS: F (.t .Y .tfNd .tfDays .dvKms .method .status .statusName).
assert(isfield(row, 'z8') && numel(row.z8) == 8, 'costate_catalog_extremes_movies:row', ...
       'a row needs a z8 of 8 elements');
z8 = row.z8(:);
J  = [];  if isfield(row, 'junctions'), J = row.junctions; end
hasJ = ~isempty(J) && size(J, 1) == 14 && all(isfinite(J(:)));
if strcmp(method, 'junctions') && ~hasJ
    error('costate_catalog_extremes_movies:noJunctions', 'the row carries no 14 x K junction states');
end
useJ = hasJ && ~strcmp(method, 'z8');
tStar = cat_.constants.tStar_s;  lStar = cat_.constants.lStar_km;
m0kg  = cat_.thruster.m0_kg;     cnd   = cat_.thruster.c_nd;   mu = cat_.constants.muStar;
Tnd   = (thrustN/m0kg)*tStar^2/(lStar*1000);
tf    = z8(8);
if useJ
    K = size(J, 2);  h = tf/K;  t = [];  Y = [];
    for k = 1:K
        [tk, Yk] = pumpkyn.cr3bp.tfMinProp(h, J(:, k), Tnd, cnd, mu);
        if k > 1, tk = tk(2:end);  Yk = Yk(2:end, :); end
        t = [t; (k - 1)*h + tk];  Y = [Y; Yk];
    end
    F.method = 'junctions';
else
    assert(isfield(row, 'sD'), 'costate_catalog_extremes_movies:row', 'a z8 flight needs the departure phase sD');
    rv0 = interp1(O.tD/O.tD(end), O.rvD, mod(row.sD, 1), 'spline');
    [t, Y] = pumpkyn.cr3bp.tfMinProp(tf, [rv0(1:6).'; 1; z8(1:7)], Tnd, cnd, mu);
    F.method = 'z8';
end
mf = Y(end, 7);
F.t = t;  F.Y = Y;  F.tfNd = tf;  F.tfDays = tf*tStar/86400;
F.dvKms = cnd*log(1/mf)*lStar/tStar;
F.status = NaN;  F.statusName = 'not annotated';
if isfield(row, 'status') && ~isempty(row.status) && isfinite(row.status)
    F.status = double(row.status);
    Kst = status_key();  kk = find(Kst.codes == F.status, 1);
    if ~isempty(kk), F.statusName = Kst.names{kk}; end
end
if isfield(row, 'statusName') && ~isempty(row.statusName), F.statusName = row.statusName; end
end

% ---------------------------------------------------------------------------
function img = panel_frame(fig, w)
% PANEL_FRAME  getframe forced to 720 x w by index resampling (the H.264
% shear cure from dro_scene_movie, generalized to a panel width).
% INPUTS: fig; w target width.  OUTPUTS: img [720 x w x 3 uint8].
F = getframe(fig);
ri = max(1, min(size(F.cdata,1), round(linspace(1, size(F.cdata,1), 720))));
ci = max(1, min(size(F.cdata,2), round(linspace(1, size(F.cdata,2), w))));
img = F.cdata(ri, ci, :);
end

% ---------------------------------------------------------------------------
function ks = find_sheet(cat_, e)
% FIND_SHEET  Index of the sheet an extreme came from (matched on tauDRO,
% Np and the arrival family/period).
% INPUTS: catalog; extreme struct.  OUTPUTS: ks [scalar] sheet index.
ks = 0;
for kk = 1:numel(cat_.sheets)
    sh = cat_.sheets(kk);
    if abs(sh.tauDRO - e.tauDRO) > 1e-9 || ~isequaln(sh.Np, e.Np), continue, end
    if isfield(sh, 'arr_params') && isfield(sh.arr_params, 'tau') && isfield(e.arr_params, 'tau') ...
            && abs(sh.arr_params.tau - e.arr_params.tau) > 1e-9, continue, end
    ks = kk;  return
end
error('costate_catalog_extremes_movies:sheet', 'no sheet matches the extreme (tauDRO %.4f)', e.tauDRO);
end
