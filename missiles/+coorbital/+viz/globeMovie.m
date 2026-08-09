function mv = globeMovie(traj,veh,env,opts)
%% Purpose:
%
%  Write an MP4 of the trajectory DEVELOPING over the Earth: the planet, the
%  path growing one frame at a time, the vehicle marked at its current
%  position with the flown track trailing behind it, and a time and altitude
%  readout that moves with it. The companion to coorbital.viz.globe3D, which
%  draws the same scene as one still arc. A still says where the vehicle went;
%  this says in what order and how fast, which is the part of a boost-glide
%  chain that a static picture cannot carry.
%
%  The planet is drawn by the package-private earthSurface, the SAME routine
%  coorbital.viz.globe3D uses, so the still and the movie cannot end up showing
%  two different Earths. The phase colours are lines() in phase order, again as
%  in globe3D, so a reader can put the two figures side by side.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun. Reads
%                                               t [N x 1] (s), x(:,1) radius
%                                               (m), x(:,2) longitude (rad),
%                                               x(:,3) latitude (rad), and
%                                               phaseIdx [N x 1]
%
%  veh              Struct                      Vehicle parameters. NOT READ.
%                                               Present so that every function
%                                               in coorbital.viz takes the same
%                                               (traj,veh,env,opts) arguments,
%                                               exactly as in
%                                               coorbital.viz.globe3D
%
%  env              Struct                      Environment model handles. NOT
%                                               READ. The arc is the state's
%                                               own geometry; no model enters
%
%  opts             Struct, optional            All fields optional:
%                                               File      char, output path.
%                                                         MUST end .mp4.
%                                                         Default is
%                                                         globeMovie.mp4 in
%                                                         tempdir -- see Notes
%                                               NFrame    [1 x 1] frames to
%                                                         render (-), default
%                                                         240. See Notes for
%                                                         what that costs
%                                               FrameRate [1 x 1] playback rate
%                                                         (frames/s), default
%                                                         20
%                                               Size      [1 x 2] figure size
%                                                         [width height]
%                                                         (pixels), default
%                                                         [1280 720]
%                                               SpinDeg   [1 x 1] azimuth swept
%                                                         by the camera over
%                                                         the WHOLE movie
%                                                         (deg), default 30.
%                                                         0 holds the camera
%                                                         still. The sweep is
%                                                         CENTRED on the arc
%                                                         midpoint, not started
%                                                         there
%                                               ViewOffsetDeg
%                                                         [1 x 1] camera swing
%                                                         ALONG the arc, away
%                                                         from its midpoint
%                                                         (deg), default 0.
%                                                         See the arc-frame
%                                                         note below
%                                               ViewTiltDeg
%                                                         [1 x 1] camera lift
%                                                         OUT of the
%                                                         trajectory plane
%                                                         (deg), default 32.
%                                                         0 puts the lens in
%                                                         the plane of the
%                                                         flight, which
%                                                         flattens the ground
%                                                         track onto a line
%                                               Inset     [1 x 1] logical,
%                                                         default true. Draws
%                                                         the growing
%                                                         altitude-vs-downrange
%                                                         panel, Tag
%                                                         'insetAxes'. false
%                                                         leaves the globe on
%                                                         its own
%                                               AltScale  [1 x 1] altitude
%                                                         exaggeration (-),
%                                                         default 1. Same
%                                                         meaning, and same
%                                                         warning, as in
%                                                         coorbital.viz.globe3D
%                                               Title     char, axes title. A
%                                                         non-unity AltScale is
%                                                         appended to it
%                                               PhaseName Cellstr, one name per
%                                                         phase for the readout.
%                                                         Default 'phase 1' ...
%                                               Texture   'auto' (default) or
%                                                         'plain'. See
%                                                         DEGRADATION below
%                                               Sky       'auto' (default) or
%                                                         'black'. See
%                                                         DEGRADATION below
%                                               FrameFcn  Function handle
%                                                         fn(hAx,kFrame) called
%                                                         once per frame. See
%                                                         Notes
%
%% Outputs:
%
%  mv               Struct                      What was written:
%                                               file       char, full path of
%                                                          the MP4
%                                               nFrame     [1 x 1] frames
%                                                          written (-)
%                                               frameIdx   [nFrame x 1] index
%                                                          of the trajectory
%                                                          SAMPLE shown in each
%                                                          frame (-)
%                                               frameRate  [1 x 1] playback
%                                                          rate (frames/s)
%                                               frameSize  [1 x 2] pixels
%                                                          actually written
%                                                          [width height]
%                                               texture    char, the Earth
%                                                          texture used, or
%                                                          'plain'
%                                               background char, the starfield
%                                                          image used, or
%                                                          'black'
%
%% Notes:
%
%  PLANET-FIXED, NOT INERTIAL. x = r cos(lat) cos(lon), y = r cos(lat) sin(lon),
%  z = r sin(lat), with longitude taken from the state in whatever frame the run
%  was flown in. Every run this library ships is flown with env.omegaE = 0, so
%  those are PLANET-FIXED coordinates and the arc sits over the ground it
%  actually flew over. A reader who assumes they are ECI will be wrong: nothing
%  here precesses the frame, and with a rotating Earth the two differ by the
%  sidereal angle accumulated over the flight.
%
%  DEGRADATION -- WHAT YOU GET WHEN THE PRETTY ASSETS ARE NOT THERE. This
%  function needs no network, ever, and runs headless under matlab -batch. The
%  Earth texture and the starfield are OPTIONAL DECORATION taken from the
%  pumpkyn toolbox when that toolbox is on the MATLAB path, and there are two
%  fallbacks, both of which are a fully working movie:
%
%    Texture 'auto'.  pumpkyn's earth-clouds-4k.jpg, texture-mapped onto a
%                     180-cell sphere, when it can be found.
%    ...falls back to A PLAIN PALE-GREY SPHERE with a 10-degree graticule --
%                     earthSurface's own look, exactly what globe3D draws. If
%                     your movie shows a plain grey ball and you wanted the
%                     blue marble, pumpkyn is not on your path: add
%                     <pumpkyn>/src to it, or pass Texture 'plain' to say you
%                     meant it.
%
%    Sky 'auto'.      pumpkyn.util.stars3D with starmap_4k.jpg on a celestial
%                     sphere behind the scene, when both resolve.
%    ...falls back to A FLAT BLACK BACKGROUND, which is what stars3D's own
%                     'black' mode produces. Nothing is missing from the data;
%                     the sky is a backdrop, not a measurement.
%
%  Both fallbacks are recorded in the returned mv.texture and mv.background, so
%  a script that cares can check which picture it got instead of guessing.
%
%  WHY NOT pumpkyn.util.earth3D, which renders a far prettier planet. Three
%  measured reasons, 2026-08-07, pumpkyn at proj7/external/pumpkyn: it CREATES
%  ITS OWN FIGURE, which would break this package's rule that private/vizParent
%  is the only place allowed to call figure(); only its 'day' texture resolves
%  locally, its 'night' and 'clouds' images not being in the tree; and its
%  stars option calls an undefined star3D and throws. Its texture FILE is
%  excellent and is what 'auto' above loads. Its renderer is not needed.
%
%  FRAMES LAND ON SAMPLES. Frame k shows the last trajectory sample at or
%  before linspace(t(1),t(end),NFrame), and the readout states THAT SAMPLE's
%  time, not the frame's nominal one. Nothing is interpolated, so what the
%  movie draws is always trajectory the propagator actually produced. Frame 1
%  is sample 1 and the last frame is sample N whatever NFrame is.
%
%  COST, AND WHY THE DEFAULT IS 240 FRAMES. MEASURED, not budgeted, on this
%  machine on 2026-08-07 over the full boost-glide-descent chain with the
%  texture and the starfield: 0.08 to 0.10 s per frame at the default
%  1280 x 720, and 0.17 s per frame at 1920 x 1080. The scene is one textured
%  sphere and a handful of lines, which is an order of magnitude cheaper than
%  the multi-body scenes that gave the 2 to 3 s per frame figure this library's
%  plan quoted. SO THE DEFAULT NFrame OF 240 COSTS 19 TO 24 s AT 1280 x 720, OR
%  ABOUT 40 s AT 1920 x 1080, and plays for 12 s at the default 20 frames/s --
%  both measured, not extrapolated. The first cut of this function defaulted to
%  90, which bought a 4.5 s movie for 8 s of work; that is the wrong trade when
%  the whole point of an animation is that a reader can follow the flight. For
%  a smoother-looking flight raise NFrame, not FrameRate: more frames is a
%  finer walk along the trajectory, a higher frame rate just plays the same
%  walk faster.
%
%  READING THREE PHASES OUT OF A PICTURE OF ONE. The three legs of a
%  boost-glide chain are nothing like the same length -- HGV/run_boost_glide
%  flies 114 km of boost, 7528 km of glide and 30 km of descent -- so drawn at
%  one width on a whole globe the descent is a few pixels and the movie reads
%  as two phases and a rounding error. It was measured at FIVE yellow pixels in
%  a 1280 x 720 frame. Three things are done about it, all of them visible in
%  the frame rather than in the data: the LINE WIDTH is graded by phase LENGTH,
%  so the shortest phase is drawn boldest and the longest at the base width; a
%  HANDOFF MARKER, Tag 'handoffMarker', is dropped at each junction in the
%  colour of the phase that starts there, and appears only once the flight
%  reaches it; and a LEGEND names the phases in their own colours. A reader who
%  wants the descent bigger than that wants a second movie zoomed on it, which
%  this function does not do.
%
%  MARKERS RIDE A SHELL 2 PER CENT ABOVE THE SPHERE, AND WHY IT IS THAT BIG. A
%  marker drawn AT the surface comes out as a HALF-DISC: MATLAB draws a marker
%  as a screen-aligned quad at one depth, the planet under it recedes across
%  that quad, and the far half loses the depth test. So every marker is drawn
%  at radius max(plotted radius, 1.02 rE).
%
%  The 2 per cent was MEASURED, not chosen. The lift a marker needs is its own
%  screen radius converted to kilometres times the tangent of its angle from
%  the centre of the disc -- at the default 1280 x 720 this scene runs about
%  23 km per pixel, so the 10 pt vehicle marker 9 deg off centre needs ~19 km
%  and the 14 pt launch ring 55 deg off centre needs ~230 km. A shell was swept
%  at 1.002, 1.01, 1.02 and 1.04 rE with the frames inspected each time:
%  1.002 (12.8 km) left BOTH markers halved, 1.01 cured the vehicle and not the
%  launch ring, 1.02 (128 km) cures both, 1.04 costs visible displacement and
%  buys nothing further. One shell serves the LARGER frames because the markers
%  grow with the frame: marker pixels rise and kilometres-per-pixel falls in
%  the same proportion, so the kilometres a marker spans stop changing with
%  Size once past the floors.
%
%  What it costs, stated plainly: the shell is radial, so a marker's apparent
%  GROUND position shifts by up to 128*sin(angle from disc centre) km, about
%  105 km and four pixels at the worst point of the shipped run. Nothing moves
%  along the ground in the DATA -- coorbital.viz.globe3D draws the same markers
%  exactly on the surface and is the figure to read a position off. And the
%  rule cannot be perfect: the required lift grows without bound towards the
%  limb, and the marker floors mean a very small frame carries relatively
%  larger markers than the sweep assumed, so a marker near the limb of a very
%  small frame can still come out clipped. Measured residual: at 480 x 320 the
%  launch ring 55 deg off centre is still cut, and it is cut at the DEFAULT
%  1280 x 720 as well once the spin has carried it far enough towards the limb
%  -- it is visible in the final frame of the default render. The residual is
%  not confined to the largest or the smallest frame size. That is the known
%  edge of this fix, not an unnoticed one.
%
%  The launch marker is additionally drawn hollow and larger than the vehicle,
%  so at t = 0, when the two are the same point, the picture is a white dot
%  inside a green ring rather than one white dot with the ring hidden under it.
%
%  THE DEFAULT OUTPUT GOES TO tempdir, not to the current folder. A movie is a
%  multi-megabyte binary and a self-demo that dropped one into whatever
%  directory happened to be current is a self-demo that dirties a repository.
%  Pass File to put it somewhere you chose.
%
%  FrameFcn is the seam that makes this function testable. Everything a movie
%  draws is gone by the time it returns, so a caller -- tests/test_viz among
%  them -- gets one call per frame, fn(hAx,kFrame), made AFTER that frame's
%  graphics are updated and BEFORE it is captured, and can read the drawn
%  objects back out by Tag or add an overlay of its own.
%
%  Tags on the drawn objects, for anyone reading them back. In the AXES:
%  'earthSurface' the planet, 'globeTrack' one line per phase with the phase
%  index in its UserData, 'handoffMarker' one per junction with the index of
%  the phase it starts in its UserData, 'vehicleMarker' the current position,
%  'launchMarker' the start point. In the FIGURE: 'titleText' the caption,
%  'hudText' the readout -- both annotations rather than axes children, see the
%  note beside them for the measured reason -- and 'insetAxes' the altitude
%  panel, which is a SECOND AXES rather than a line and is present only when
%  Inset is true. INSIDE that panel: 'insetBackdrop' its semi-transparent
%  card, 'insetTrack' one line per phase with the phase index in its UserData,
%  and 'insetVehicle' the current point on the profile.
%
%  The axes are invisible, as pumpkyn.util.earth3D leaves them, because a tick
%  grid around a cinematic globe reads as clutter. The units are therefore
%  stated where the numbers are: the readout labels its own seconds and
%  kilometres, and the caption carries any altitude exaggeration.
%
%  THE TEXT IS SIZED FROM THE FRAME, not fixed in points. Absolute point sizes
%  against fractional boxes work at one frame size and break at every other:
%  the first cut of this function used 14 pt for the caption and 12 pt for the
%  readout, which are right at 1280 x 720 and wrong at the 480 x 320 the
%  shipped test renders, where the caption wrapped onto the globe and the
%  readout put its 's' and its 'km' on lines of their own. Both are now
%  round(Size(2)/50) and round(Size(2)/60), floored at 8 pt, so the layout
%  holds from 480 x 320 to 1920 x 1080. The test cannot see this -- it reads
%  the String property, which was correct throughout -- so it was found, and
%  must be re-checked, by rendering a frame and LOOKING at it.
%
%  THE SAME DEFECT RECURRED IN THE INSET, twice, and both cures are geometric
%  rather than cosmetic. First, the inset's x label was CLIPPED off the bottom
%  of every frame under about 500 px high -- which is both the 480 x 320 the
%  shipped test renders and the 640 x 360 the self-demo does -- because the
%  panel sat at a FRACTIONAL height of 0.055 while the label it has to clear is
%  an absolute stack of tick text and axis label about 26 px tall. The origin
%  is therefore max(0.055,26/Size(2)): fractional where there is room, pixels
%  where there is not.
%
%  Second, the inset was drawn OPAQUE, and an opaque panel over a globe hides
%  whatever the flight put behind it. On the shipped run_target case the glide
%  descended behind its left edge and re-emerged above its top -- a movie of a
%  trajectory concealing the trajectory. The panel's axes Color is now 'none'
%  and its backdrop is a semi-transparent patch inside it, so the track shows
%  through: the arc stays readable, the curve stays legible against the
%  planet, and nothing that flew is hidden. Both were found by rendering and
%  looking, and both have to be re-checked the same way.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A synthetic two-phase arc, so the demo costs no integration, and
%% a small cheap movie, so it costs no patience either:
if nargin == 0
                 c = coorbital.util.missileConst();
              nDem = 300;
              lam  = linspace(0,deg2rad(70),nDem)';
           demo.x  = [c.rE + 80e3.*sin(pi.*lam./max(lam)), ...
                      lam, ...
                      deg2rad(35).*sin(lam), ...
                      6000.*ones(nDem,1), ...
                      zeros(nDem,1), ...
                      zeros(nDem,1)];
           demo.t  = linspace(0,2400,nDem)';
           demo.u  = zeros(nDem,2);
    demo.phaseIdx  = 1 + (lam > deg2rad(35));
    demo.junction  = repmat(struct('t',[],'x',[]),1,1);
          demoOpt  = struct('NFrame',24,'Size',[640 360],'AltScale',30, ...
                            'PhaseName',{{'glide','descent'}}, ...
                            'Title','coorbital.viz.globeMovie self-demo');
                mv = coorbital.viz.globeMovie(demo, ...
                         coorbital.util.vehicleDefaults(),struct(),demoOpt);
    fprintf('globeMovie self-demo wrote %d frames to %s\n',mv.nFrame,mv.file);
    return;
end

%% Options, every one of them defaulted before anything is drawn:
    if nargin < 4 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
           outFile = vizOption(opts,'File',fullfile(tempdir,'globeMovie.mp4'));
            nFrame = vizOption(opts,'NFrame',240);
            fpsOut = vizOption(opts,'FrameRate',20);
            sizPix = vizOption(opts,'Size',[1280 720]);
            spinDg = vizOption(opts,'SpinDeg',30);
           viewOff = vizOption(opts,'ViewOffsetDeg',0);
           viewTlt = vizOption(opts,'ViewTiltDeg',32);
           insetOn = vizOption(opts,'Inset',true);
          altScale = vizOption(opts,'AltScale',1);
            titTxt = vizOption(opts,'Title','Trajectory over the Earth');
            phName = vizOption(opts,'PhaseName',{});
            texOpt = vizOption(opts,'Texture','auto');
            skyOpt = vizOption(opts,'Sky','auto');
            frmFcn = vizOption(opts,'FrameFcn',[]);

    assert(isscalar(nFrame) && nFrame >= 2 && nFrame == fix(nFrame), ...
        'NFrame must be a whole number of at least 2; %s given.',mat2str(nFrame));
    assert(isscalar(fpsOut) && fpsOut > 0, ...
        'FrameRate must be a positive scalar; %s given.',mat2str(fpsOut));
    assert(numel(sizPix) == 2 && all(sizPix >= 64) && all(sizPix == fix(sizPix)), ...
        'Size must be [width height] in whole pixels, both at least 64.');
    assert(isscalar(spinDg) && isfinite(spinDg), ...
        'SpinDeg must be a finite scalar; %s given.',mat2str(spinDg));

%% The three camera-and-inset options are checked exactly as hard as the rest.
%% They were shipped unchecked, and an unchecked angle is the worst of the
%% three: ViewTiltDeg = NaN produces a NaN camera direction, view() silently
%% ignores it, and the movie comes out on MATLAB's default viewpoint with
%% nothing anywhere saying the option was thrown away:
    assert(isscalar(viewOff) && isfinite(viewOff), ...
        'ViewOffsetDeg must be a finite scalar; %s given.',mat2str(viewOff));
    assert(isscalar(viewTlt) && isfinite(viewTlt), ...
        'ViewTiltDeg must be a finite scalar; %s given.',mat2str(viewTlt));
    assert(isscalar(insetOn) && (islogical(insetOn) || ...
                                 (isnumeric(insetOn) && isfinite(insetOn))), ...
        'Inset must be a logical scalar; %s given.',mat2str(insetOn));
              insetOn = logical(insetOn);

    assert(isscalar(altScale) && altScale > 0, ...
        'AltScale must be a positive scalar; %s given.',mat2str(altScale));
    assert(isempty(frmFcn) || isa(frmFcn,'function_handle'), ...
        'FrameFcn must be a function handle taking (hAx,kFrame), or [].');

%% Text sized FROM THE FRAME. A point size is absolute and a text box is
%% fractional, so any fixed size is right at one Size and wrong at every other:
%% at the 480 x 320 the shipped test renders, the 14 pt caption and 12 pt
%% readout this function first carried wrapped onto the globe and onto extra
%% lines. Height over 50 and over 60, floored at 8 pt so that a very small
%% frame stays legible rather than becoming correct and unreadable:
            fntCap = max(8,round(sizPix(2)./50));
            fntHud = max(8,round(sizPix(2)./60));

%% The markers GROW with the frame and never shrink below the sizes that were
%% checked by eye at 480 x 320 and 1280 x 720. Both halves of that were
%% measured. Letting them scale down as well was tried first, and at 480 x 320
%% it turned the descent handoff into a speck -- the count of descent-coloured
%% pixels outside the legend fell from 50 to 8, which is most of the way back
%% to the defect this is here to fix. Letting them stay fixed instead leaves a
%% 10 pt dot on a 1080p frame, which is a speck of a different kind. So:
%% floors at the values that read at the smallest frame, growth above them:
            mkVeh  = max(10,round(sizPix(2)./72));
            mkLau  = max(14,round(1.4.*mkVeh));
            mkHand = max( 9,round(0.9.*mkVeh));

%% ...and the vehicle dot inside the altitude inset by the same rule. It is
%% HALF the globe marker, because the inset is a fifth of the frame and a
%% marker sized for the planet covers a visible slice of the profile; and it is
%% derived from mkVeh rather than fixed at 5 pt, because a 5 pt dot that reads
%% at 720 is invisible at 1080 -- which is the whole reason nothing else in
%% this file carries an absolute point size:
            mkIns  = max( 5,round(0.5.*mkVeh));

%% Texture and Sky are refused rather than quietly downgraded when they are not
%% one of the two things they can be. A misspelt 'atuo' that silently produced
%% the grey ball would be indistinguishable from an absent pumpkyn, which is the
%% one confusion the DEGRADATION note exists to prevent:
    assert((ischar(texOpt) || isstring(texOpt)) && ...
           any(strcmpi(texOpt,{'auto','plain'})), ...
        'Texture must be ''auto'' or ''plain''.');
    assert((ischar(skyOpt) || isstring(skyOpt)) && ...
           any(strcmpi(skyOpt,{'auto','black'})), ...
        'Sky must be ''auto'' or ''black''.');

%% The output path, checked BEFORE a single frame is rendered. A movie that
%% renders for ten minutes and then cannot be written is ten minutes lost:
    [outDir,~,outExt] = fileparts(outFile);
    assert(strcmpi(outExt,'.mp4'), ...
        'File must name an .mp4; "%s" ends "%s".',outFile,outExt);
    if ~isempty(outDir)
        assert(isfolder(outDir), ...
            'the output folder "%s" does not exist.',outDir);
    end

%% The trajectory, and the geometry drawn from it. ONE conversion to
%% kilometres, and the exaggeration multiplies the ALTITUDE and never the
%% sphere, so the planet stays the right size -- see coorbital.viz.globe3D:
                 c = coorbital.util.missileConst();
             rEkm  = c.rE./1000;
                tS = traj.t(:);
                nS = numel(tS);
    assert(nS >= 2,'a movie needs at least two samples; %d given.',nS);
    assert(all(diff(tS) >= 0),'traj.t must be non-decreasing.');
             phIdx = traj.phaseIdx(:);
    assert(numel(phIdx) == nS, ...
        'phaseIdx has %d entries against %d samples.',numel(phIdx),nS);
            phList = unique(phIdx).';
              nPh  = numel(phList);

              rPlt = c.rE + altScale.*(traj.x(:,1) - c.rE);
               lon = traj.x(:,2);
               lat = traj.x(:,3);
             xTrk  = rPlt.*cos(lat).*cos(lon)./1000;
             yTrk  = rPlt.*cos(lat).*sin(lon)./1000;
             zTrk  = rPlt.*sin(lat)./1000;
             altKm = (traj.x(:,1) - c.rE)./1000;

%% ...and the same points again on a shell 2 per cent above the sphere, for the
%% MARKERS only. A marker sitting on the surface is cut in half by the depth
%% buffer and renders as a half-disc; 1.02 rE, 128 km here, is the smallest
%% lift measured to clear BOTH the vehicle marker and the larger launch ring in
%% the shipped scene -- see the shell note in the header for the sweep and for
%% what the lift costs. A marker already above the shell is left where the
%% trajectory put it. The TRACK is never lifted -- it is the data:
            shellF = 1.02;
              rMk  = max(rPlt,shellF.*c.rE);
             xMk   = rMk.*cos(lat).*cos(lon)./1000;
             yMk   = rMk.*cos(lat).*sin(lon)./1000;
             zMk   = rMk.*sin(lat)./1000;

%% Which SAMPLE each frame shows. Uniform in time rather than in index, so the
%% movie plays at a constant rate whatever the integrator's step history did:
            tFrame = linspace(tS(1),tS(end),nFrame)';
          frameIdx = zeros(nFrame,1);
    for kf = 1:nFrame
      frameIdx(kf) = find(tS <= tFrame(kf),1,'last');
    end

%% The optional pumpkyn assets, resolved once. Absent ones cost nothing but the
%% look; see DEGRADATION in the header:
            texFil = '';
    if strcmpi(texOpt,'auto')
            texFil = pumpkynAsset('earth-clouds-4k.jpg');
    end
            skyFil = '';
    if strcmpi(skyOpt,'auto')
            skyFil = pumpkynAsset('starmap_4k.jpg');
    end

%% The figure. Created through the one routine in this package allowed to
%% create one, invisible because a movie is captured rather than watched, and
%% deleted before this function returns whatever happens below:
        [hFig,hAx] = vizParent([],1,'Globe movie','off');
    set(hFig,'Color','k','Position',[80 80 sizPix(1) sizPix(2)]);
    hold(hAx,'on');

%% The planet. A 180-cell mesh under a texture so the map is not visibly
%% faceted, and earthSurface's own 36-cell graticule sphere without one:
    if ~isempty(texFil)
             hEart = earthSurface(hAx,180);
        set(hEart,'FaceColor','texturemap', ...
                  'CData',flipud(imread(texFil)), ...
                  'EdgeColor','none');
    else
        earthSurface(hAx,36);
    end

%% Limits, aspect and camera, all frozen before the first frame. A movie whose
%% axes rescale as the track grows is a movie of a wobbling planet:
               lim = 1.02.*max(rEkm,max(rPlt)./1000);
    axis(hAx,'equal');
    set(hAx,'XLim',[-lim lim],'YLim',[-lim lim],'ZLim',[-lim lim], ...
            'Color','k','Clipping','off','Visible','off');
    axis(hAx,'vis3d');

%% ...and then zoom in by root three. MATLAB frames a three-dimensional axes so
%% that the whole DATA CUBE fits whatever way it is turned, which means it fits
%% the cube's DIAGONAL and the planet is left occupying 1/sqrt(3) of the frame
%% with the corners empty. Everything drawn here lies inside the SPHERE of
%% radius lim, not just inside the cube, so that factor is pure waste and
%% recovering it clips nothing, at any spin angle.
%%
%% WHAT 1.62 ACTUALLY BOUGHT, stated honestly because the first version of this
%% comment overclaimed. It is held under the full sqrt(3) = 1.7321 to leave the
%% top of the frame less crowded, and it is strictly safer than the full
%% factor: everything sits 6.5 per cent further from the caption than it would.
%% It does NOT keep the limb out of the caption band, and no zoom factor above
%% one can. The caption occupies the top tenth of the figure and the data
%% sphere of radius lim is drawn to fill the frame, so on a polar route -- the
%% case this was changed for -- the arc still reaches into the text. Clearing
%% it properly means shrinking the AXES RECTANGLE to leave a caption strip, not
%% easing the zoom, and that would re-frame every movie this library has
%% rendered. Left as a known, measured edge.
    camzoom(hAx,1.62);

%% A viewpoint looking down on the middle of the arc: MATLAB's view(az,el)
%% points at longitude L and latitude B for az = L + 90, el = B.
%%
%% The middle is the GREAT-CIRCLE MIDPOINT OF THE ENDPOINTS, not the median
%% sample. ode45 clusters samples where the dynamics are fast -- the skips and
%% the terminal dive -- so on a long flight the median sample sits well down
%% the track, and aiming there swings the launch end round the limb and out of
%% frame. On a 97 deg intercontinental arc that hid the launch point entirely:
%% Build an ARC FRAME from the two endpoints: m points at the arc midpoint, n
%% is the trajectory-plane normal, and t runs along-track at the midpoint.
%% BOTH DEGENERACIES ARE GUARDED, and they are different degeneracies. The
%% plane normal u1 x u2 vanishes when the endpoints are coincident AND when
%% they are antipodal; the midpoint u1 + u2 vanishes only in the antipodal
%% case, and it vanishes into ROUNDOFF rather than into zero -- measured at
%% 2.0e-16 for 10N 20E against 10S 160W, which normalises to the arbitrary
%% direction [0.5547 -0.8321 0]. A camera aimed by roundoff points somewhere
%% reproducible and meaningless, with nothing in the picture saying so. So the
%% midpoint is checked as well as the normal, and antipodal endpoints aim at
%% the LAUNCH end, which is a defensible answer rather than an accidental one.
%%
%% The fallback frame is built to be ORTHONORMAL, which the first cut was not:
%% it took nP = [0;0;1] whatever m was, so m . n came out at 0.5 for a 30 deg
%% case and ViewTiltDeg was then not the angle out of any plane at all. The
%% replacement takes whichever axis is least aligned with m and removes m from
%% it, so m, n and t are a genuine right-handed triad in every branch and the
%% tilt means what the header says it means:
               uv = @(a,b) [cos(a).*cos(b); cos(a).*sin(b); sin(a)];
               u1 = uv(lat(1),lon(1));
               u2 = uv(lat(end),lon(end));
               nP = cross(u1,u2);
               mM = u1 + u2;
    if sqrt(sum(mM.^2)) < 1e-8
               mM = u1;
    end
               mM = mM./sqrt(sum(mM.^2));
    if sqrt(sum(nP.^2)) < 1e-12
               ax = [0;0;1];
        if abs(mM(3)) > 0.9
               ax = [1;0;0];
        end
               nP = cross(mM,ax);
    end
               nP = nP./sqrt(sum(nP.^2));

%% n is orthogonal to m by construction in both branches -- the cross product
%% of the endpoints is perpendicular to their sum, and the fallback cross
%% product is perpendicular to m outright -- so t closes a right-handed triad
%% without any further orthogonalisation:
               tT = cross(nP,mM);
               tT = tT./sqrt(sum(tT.^2));

%% ALTITUDE IS ONLY VISIBLE AWAY FROM THE SUB-CAMERA POINT. A trajectory point
%% sitting at angular distance theta from directly-under-the-camera shows its
%% height as h*sin(theta): straight underneath, the altitude vector points at
%% the lens and the whole vertical profile collapses to nothing. Aiming at the
%% arc midpoint therefore produces a beautiful ground track and a flat-looking
%% flight. So swing the camera ViewOffsetDeg along the arc, which throws the
%% far half of the track out towards the limb where its height projects, and
%% lift it ViewTiltDeg out of the trajectory plane so the ground track does not
%% collapse onto a line at the same time.
               cD = cosd(viewTlt).*(cosd(viewOff).*mM + sind(viewOff).*tT) ...
                    + sind(viewTlt).*nP;
               cD = cD./sqrt(sum(cD.^2));
               az0 = rad2deg(atan2(cD(2),cD(1))) + 90;
               el0 = rad2deg(asin(max(-1,min(1,cD(3)))));
    view(hAx,az0,el0);

%% The sky, drawn AFTER the camera is set because stars3D sizes its celestial
%% sphere from the camera it finds. The epoch below only orients the backdrop;
%% it is a display choice and not a physical constant, so it does not belong in
%% missileConst any more than the exaggeration factor does:
             jdSky = 2451545.0;
    if ~isempty(skyFil)
        pumpkyn.util.stars3D(hFig,hAx,'stars',skyFil,jdSky,'Quality','interactive');
            skyTxt = skyFil;
    else
        set(hFig,'Color','k');
        set(hAx,'Color','k');
            skyTxt = 'black';
    end

%% Which samples each phase is drawn over, and what it is called. Each phase is
%% extended backward to the junction so the growing track has no holes in it --
%% the same contract coorbital.viz.globe3D documents:
               col = lines(max(nPh,7));
            selAll = cell(1,nPh);
           nmPhase = cell(1,nPh);
    for kp = 1:nPh
        selAll{kp} = find(phIdx == phList(kp));
        if kp > 1
        selAll{kp} = [find(phIdx == phList(kp-1),1,'last'); selAll{kp}];
        end
       nmPhase{kp} = sprintf('phase %d',phList(kp));
        if numel(phName) >= kp
       nmPhase{kp} = phName{kp};
        end
    end

%% HOW WIDE each phase is drawn, graded by how LONG it is. The three legs of a
%% boost-glide chain differ by two orders of magnitude in length -- 114 km,
%% 7528 km and 30 km on this library's shipped run -- and at one width on a
%% whole globe the short ones vanish: the terminal descent measured FIVE
%% coloured pixels in a 1280 x 720 frame. So the longest phase is drawn at the
%% base width and a phase of vanishing length at the base plus the whole of the
%% extra, which puts the boldest stroke on exactly the leg a globe-wide view
%% would otherwise swallow. The lengths are the DRAWN lengths, chord by chord
%% through the plotted points, so the grading follows the exaggerated picture
%% the reader is actually looking at:
            lwBase = 2.0;
            lwXtra = 3.5;
             lenPh = zeros(1,nPh);
    for kp = 1:nPh
               sel = selAll{kp};
        if numel(sel) >= 2
         lenPh(kp) = sum(sqrt(diff(xTrk(sel)).^2 + ...
                              diff(yTrk(sel)).^2 + ...
                              diff(zTrk(sel)).^2));
        end
    end
            lenMax = max([lenPh,eps]);
              lwPh = lwBase + lwXtra.*(1 - lenPh./lenMax);

%% One line per phase, each empty until the frame loop fills it:
            hTrack = gobjects(1,nPh);
    for kp = 1:nPh
        hTrack(kp) = line(hAx,NaN,NaN,NaN, ...
                          'Color',col(kp,:),'LineWidth',lwPh(kp), ...
                          'Tag','globeTrack','UserData',phList(kp));
    end

%% A marker at each junction, in the colour of the phase that STARTS there, and
%% carrying that phase's index. Colour and width tell a reader that the stroke
%% changed; a marker tells them WHERE, which is the thing a 30 km terminal leg
%% cannot say for itself. Each is held off-screen until the flight reaches it,
%% in the frame loop, for the same reason there is no impact marker:
              nJct = max(nPh - 1,0);
            hHand  = gobjects(1,nJct);
            idxJct = zeros(1,nJct);
    for kp = 2:nPh
       idxJct(kp-1) = selAll{kp}(1);
       hHand(kp-1)  = line(hAx,NaN,NaN,NaN, ...
                           'LineStyle','none','Marker','d', ...
                           'MarkerSize',mkHand, ...
                           'MarkerFaceColor',col(kp,:), ...
                           'MarkerEdgeColor',[0 0 0],'LineWidth',1.2, ...
                           'Tag','handoffMarker','UserData',phList(kp));
    end

%% The launch point, on the marker shell over the ground the flight began on
%% and never exaggerated -- it marks ground, and lifting it by the exaggeration
%% would put it somewhere the vehicle was not. HOLLOW and larger than the
%% vehicle marker, so that at t = 0, when the two are the same point, the
%% picture is a white dot inside a green ring rather than one white dot. There
%% is no impact marker: the vehicle has not got there yet, and drawing where it
%% will end up would give the ending away in frame one:
    line(hAx,shellF.*rEkm.*cos(lat(1)).*cos(lon(1)), ...
             shellF.*rEkm.*cos(lat(1)).*sin(lon(1)), ...
             shellF.*rEkm.*sin(lat(1)), ...
         'LineStyle','none','Marker','o','MarkerSize',mkLau,'LineWidth',2.0, ...
         'Color',[0.2 1 0.4],'Tag','launchMarker');

%% The vehicle. WHITE, and deliberately not one of the phase colours: lines()
%% reaches yellow at phase 3, and a yellow vehicle riding the end of a yellow
%% descent is a marker a reader cannot find:
              hVeh = line(hAx,NaN,NaN,NaN, ...
                          'LineStyle','none','Marker','o','MarkerSize',mkVeh, ...
                          'MarkerFaceColor',[1 1 1], ...
                          'MarkerEdgeColor',[0 0 0],'LineWidth',1.2, ...
                          'Tag','vehicleMarker');

%% The legend, naming the phases in their own colours and at their own widths.
%% AutoUpdate off and an explicit handle list, so the markers created above and
%% anything a FrameFcn adds later stay out of it. It is the one part of the
%% frame that says what the colours MEAN; the readout names only the phase the
%% vehicle is in at that instant:
              hLeg = legend(hAx,hTrack,nmPhase,'AutoUpdate','off', ...
                            'Location','southwest','Interpreter','none');
    set(hLeg,'TextColor',[1 1 1],'Color',[0.06 0.06 0.06], ...
             'EdgeColor',[0.45 0.45 0.45],'FontSize',fntHud,'Box','on');

%% The caption and the readout, both FIGURE-level annotations rather than an
%% axes title and an axes text. Measured on R2025b, 2026-08-07: in an axes
%% carrying axis vis3d, one normalized unit of a text child spans the projected
%% DATA CUBE and not the axes rectangle -- 637 px inside a 540 px figure in the
%% case checked -- so a caption placed at 0.97 lands off the top of the frame,
%% and the axes title, positioned the same way, goes with it. Both were
%% silently clipped out of the first cut of this function. An annotation is
%% positioned in FIGURE coordinates and cannot be moved by the camera:
            titFul = titTxt;
    if altScale ~= 1
            titFul = sprintf('%s (altitude exaggerated %gx)',titTxt,altScale);
    end
%% The caption carries the SAME semi-transparent backdrop the inset panel and
%% the legend already use, and for the same measured reason. White bold text
%% at the top centre of the frame crosses whatever the camera has put there,
%% and on a high-latitude arc that is the polar ice cap: rendered bare, the
%% middle of the caption is white on white and simply cannot be read. Measured
%% on the Plesetsk-to-New-York arc, 2026-08-09 -- the launch coordinates were
%% illegible for the whole flight while the two ends of the same line, over
%% sky, were crisp. The alpha is low enough that the limb still reads through
%% the strip, which is why this is a backdrop and not an opaque bar:
    annotation(hFig,'textbox',[0.03 0.90 0.94 0.09], ...
               'String',titFul,'Tag','titleText', ...
               'Color',[1 1 1],'FontSize',fntCap,'FontWeight','bold', ...
               'HorizontalAlignment','center','VerticalAlignment','top', ...
               'BackgroundColor',[0.06 0.06 0.06],'FaceAlpha',0.55, ...
               'EdgeColor','none','FitBoxToText','off','Interpreter','none');
              hHud = annotation(hFig,'textbox',[0.015 0.68 0.42 0.22], ...
                          'String',{''},'Tag','hudText', ...
                          'Color',[1 1 1],'FontName','FixedWidth', ...
                          'FontSize',fntHud, ...
                          'HorizontalAlignment','left','VerticalAlignment','top', ...
                          'EdgeColor','none','FitBoxToText','off','Interpreter','none');

%% THE ALTITUDE INSET. The globe shows where the vehicle went; it cannot show
%% how high, because at the sub-camera point the altitude vector points at the
%% lens and the vertical profile collapses. A small two-dimensional panel is a
%% far better answer than tilting the camera until both are equally bad: this
%% one plots TRUE altitude, unexaggerated, against downrange, and grows with
%% the flight exactly as the track above it does.
              hIns = [];
              hInL = [];
              hInV = [];
    if insetOn
%% Downrange along the great circle from the launch point, per sample. In
%% KILOMETRES, the same unit the label states and the same unit the altitude
%% is drawn in, because this panel's whole claim is that it is at true scale:
              dRng = coorbital.util.greatCircle(lat(1),lon(1),lat,lon) ...
                     .*(rEkm);

%% Fixed limits, set from the WHOLE flight, so the curve grows into a frame
%% that never rescales -- a rescaling axis makes a growing line look static.
%% Both are FLOORED at one kilometre: a trajectory that never leaves the
%% surface, or never leaves the pad, gives a zero span, and MATLAB refuses a
%% zero-width axis limit with an error. That error used to be thrown here with
%% the video writer already open, so a flat run left a figure on screen and a
%% part-written MP4 on disk; with Inset false the same run rendered perfectly:
             xLimI = [0 max(max(dRng).*1.02,1)];
             yLimI = [0 max(max(altKm).*1.12,1)];

%% The panel itself. Its origin is max(0.055,26/height) rather than a flat
%% 0.055 because the x label under it is an absolute 26 px of tick text and
%% axis label: below about 500 px of figure height a fractional origin puts the
%% label off the bottom of the frame. And its Color is 'none' rather than
%% black, so the globe and the flown track stay visible THROUGH it -- see the
%% inset note in the header for what an opaque panel hid:
              yIns = max(0.055,26./sizPix(2));
              fntIns = max(7,round(sizPix(2)./85));
              hIns = axes(hFig,'Position',[0.615 yIns 0.345 0.235], ...
                          'Color','none','XColor',[1 1 1],'YColor',[1 1 1], ...
                          'GridColor',[0.6 0.6 0.6],'GridAlpha',0.35, ...
                          'XLim',xLimI,'YLim',yLimI,'FontSize',fntIns, ...
                          'Box','on','Tag','insetAxes');
        hold(hIns,'on'); grid(hIns,'on');

%% A semi-transparent backdrop, drawn FIRST so everything else sits on top of
%% it. Transparent enough that the track behind the panel is still followable,
%% dark enough that the profile curve reads against a bright Earth texture.
%%
%% It is drawn with CLIPPING OFF and overhangs the axes box, because the panel
%% is wider than its plot area: the tick labels, the two axis labels and the
%% title all sit OUTSIDE the box, and white text on the pale limb of the planet
%% is unreadable. The overhangs are set in PIXELS and converted to data units,
%% for the same reason the origin above is -- text is absolute and the axes is
%% fractional, so a margin quoted as a fraction of the data span is right at
%% one frame size and wrong at every other:
             axWpx = 0.345.*sizPix(1);
             axHpx = 0.235.*sizPix(2);
             padLx = (3.5.*fntIns + 6)./axWpx.*xLimI(2);
             padRx = 6./axWpx.*xLimI(2);
             padBy = max(26,2.2.*fntIns + 10)./axHpx.*yLimI(2);
             padTy = (1.9.*max(7,round(sizPix(2)./80)) + 6)./axHpx.*yLimI(2);
        patch(hIns,'XData',[-padLx, xLimI(2) + padRx, ...
                            xLimI(2) + padRx, -padLx], ...
                   'YData',[-padBy, -padBy, ...
                            yLimI(2) + padTy, yLimI(2) + padTy], ...
              'FaceColor',[0 0 0],'FaceAlpha',0.55,'EdgeColor','none', ...
              'Clipping','off','Tag','insetBackdrop');

%% One line per phase so the inset colours match the globe:
              hInL = gobjects(nPh,1);
        for kp = 1:nPh
            hInL(kp) = plot(hIns,NaN,NaN,'-', ...
                            'Color',col(kp,:),'LineWidth',lwPh(kp), ...
                            'Tag','insetTrack','UserData',phList(kp));
        end
              hInV = plot(hIns,NaN,NaN,'o','MarkerSize',mkIns, ...
                          'MarkerFaceColor',[1 1 1], ...
                          'MarkerEdgeColor',[0 0 0],'LineWidth',0.5, ...
                          'Tag','insetVehicle');
        xlabel(hIns,'downrange (km)','FontSize',fntIns);
        ylabel(hIns,'altitude (km)','FontSize',fntIns);
        title(hIns,'altitude profile (true scale)', ...
              'Color',[1 1 1],'FontSize',max(7,round(sizPix(2)./80)));
    end

%% The writer, opened LAST so that nothing can fail after a part-written file
%% exists on disk. "Last" is the invariant the comment always claimed and did
%% not hold: the inset block above was inserted between the open() and the try,
%% where a throw stranded both an open writer and a figure. It is above the
%% open() again, and the limits it sets are floored so the case that threw
%% cannot:
                vw = VideoWriter(outFile,'MPEG-4');
      vw.FrameRate = fpsOut;
        vw.Quality = 95;
    open(vw);

%% The frame loop. Everything that can throw is inside the try, because a
%% half-rendered movie must still leave the writer closed and no figure open:
           frmSize = [];
    try
        for kf = 1:nFrame
              kEnd = frameIdx(kf);

%% The track, grown to this frame and no further. Each phase keeps its own
%% colour, so boost, glide and descent stay distinguishable as they appear:
            for kp = 1:nPh
               sel = selAll{kp};
               sel = sel(sel <= kEnd);
                set(hTrack(kp),'XData',xTrk(sel), ...
                               'YData',yTrk(sel), ...
                               'ZData',zTrk(sel));
            end

%% The junction markers, each appearing only in the frame that reaches it. A
%% handoff drawn before the vehicle gets there would announce the ending, which
%% is the same reason there is no impact marker:
            for kj = 1:nJct
                if idxJct(kj) <= kEnd
                    set(hHand(kj),'XData',xMk(idxJct(kj)), ...
                                  'YData',yMk(idxJct(kj)), ...
                                  'ZData',zMk(idxJct(kj)));
                else
                    set(hHand(kj),'XData',NaN,'YData',NaN,'ZData',NaN);
                end
            end

%% The vehicle and the readout, at the sample this frame shows. The marker sits
%% on the shell, so it stays a whole disc when the vehicle is on the deck:
            set(hVeh,'XData',xMk(kEnd),'YData',yMk(kEnd),'ZData',zMk(kEnd));
               kNow = find(phList == phIdx(kEnd),1);
            set(hHud,'String',{sprintf('t   = %7.1f s',tS(kEnd)), ...
                               sprintf('alt = %7.2f km',altKm(kEnd)), ...
                               nmPhase{kNow}});

%% The inset grows on the same sample index as the globe track, so the two can
%% never drift apart:
            if insetOn
                for kp = 1:nPh
                   selI = selAll{kp};
                   selI = selI(selI <= kEnd);
                    set(hInL(kp),'XData',dRng(selI),'YData',altKm(selI));
                end
                set(hInV,'XData',dRng(kEnd),'YData',altKm(kEnd));
            end

%% The camera, sweeping SpinDeg of azimuth across the whole movie, CENTRED on
%% the arc midpoint rather than starting there -- a one-sided sweep walks a
%% long arc off the limb by the end of the movie:
            view(hAx,az0 - spinDg./2 + spinDg.*(kf - 1)./(nFrame - 1),el0);

%% The caller's per-frame hook, before the capture so that whatever it draws
%% is in the frame:
            if ~isempty(frmFcn)
                frmFcn(hAx,kf);
            end

%% drawnow BEFORE getframe, every single frame. Without it getframe captures
%% whatever the renderer had finished, which is a partially drawn figure:
            drawnow;
               frm = getframe(hFig);

%% Every frame must be the size of the first, and MPEG-4 needs both dimensions
%% even, so the first frame sets an even size and the rest are cropped to it:
            if isempty(frmSize)
           frmSize = 2.*floor([size(frm.cdata,1),size(frm.cdata,2)]./2);
            end
            assert(size(frm.cdata,1) >= frmSize(1) && ...
                   size(frm.cdata,2) >= frmSize(2), ...
                ['frame %d came back %d x %d against the first frame''s ' ...
                 '%d x %d; the figure was resized mid-render.'], ...
                kf,size(frm.cdata,1),size(frm.cdata,2),frmSize(1),frmSize(2));
            writeVideo(vw,frm.cdata(1:frmSize(1),1:frmSize(2),:));
        end
    catch err
        close(vw);
        delete(hFig);
        rethrow(err);
    end
    close(vw);
    delete(hFig);

%% What was written, including WHICH picture came out, so a caller that cares
%% whether it got the blue marble or the grey ball can ask instead of guess:
           mv.file = outFile;
         mv.nFrame = nFrame;
       mv.frameIdx = frameIdx;
      mv.frameRate = fpsOut;
      mv.frameSize = [frmSize(2),frmSize(1)];
       mv.texture  = 'plain';
    if ~isempty(texFil)
       mv.texture  = texFil;
    end
    mv.background  = skyTxt;
end

function fPath = pumpkynAsset(fileName)
%% Purpose:
%
%  Locate one image shipped inside pumpkyn's +util folder, or return empty when
%  pumpkyn is not on this MATLAB path. This is the whole of the optional
%  dependency: the missile library never calls a pumpkyn function to draw
%  anything, it only borrows two pictures when they happen to be reachable.
%
%  which() is used rather than exist(). Measured on R2025b, 2026-08-07:
%  exist('pumpkyn.util.stars3D','file') returns 0 even with the toolbox on the
%  path and the function callable, because exist does not resolve
%  package-qualified names. which() does, and returns '' when it cannot.
%
%% Inputs:
%
%  fileName         Char [1 x n]                Bare image file name, e.g.
%                                               'earth-clouds-4k.jpg'
%
%% Outputs:
%
%  fPath            Char [1 x m]                Full path to the file, or ''
%                                               when pumpkyn is absent or the
%                                               file is not in its tree
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             fPath = '';
            anchor = which('pumpkyn.util.stars3D');
    if isempty(anchor)
        return;
    end
              cand = fullfile(fileparts(anchor),fileName);
    if isfile(cand)
             fPath = cand;
    end
end
