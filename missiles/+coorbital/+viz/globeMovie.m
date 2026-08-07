function mv = globeMovie(traj,opts)
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
%  opts             Struct, optional            All fields optional:
%                                               File      char, output path.
%                                                         MUST end .mp4.
%                                                         Default is
%                                                         globeMovie.mp4 in
%                                                         tempdir -- see Notes
%                                               NFrame    [1 x 1] frames to
%                                                         render (-), default
%                                                         90. See Notes for the
%                                                         cost
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
%                                                         still
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
%  COST, AND HOW TO MAKE IT LONGER. MEASURED, not budgeted: 0.13 s per frame at
%  1920 x 1080 with the texture and the starfield, 90 frames of the full
%  boost-glide-descent chain in 11.6 s wall clock, on this machine on
%  2026-08-07. The scene is one textured sphere and a handful of lines, which
%  is an order of magnitude cheaper than the multi-body scenes that gave the
%  2 to 3 s per frame figure this library's plan quoted. So make it as long as
%  you like: NFrame 300 at 1080p is under a minute. For a smoother-looking
%  flight raise NFrame, not FrameRate -- more frames is a finer walk along the
%  trajectory, a higher frame rate just plays the same walk faster.
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
%  index in its UserData, 'vehicleMarker' the current position, 'launchMarker'
%  the start point on the surface. In the FIGURE: 'titleText' the caption and
%  'hudText' the readout, both annotations rather than axes children -- see the
%  note beside them for the measured reason.
%
%  The axes are invisible, as pumpkyn.util.earth3D leaves them, because a tick
%  grid around a cinematic globe reads as clutter. The units are therefore
%  stated where the numbers are: the readout labels its own seconds and
%  kilometres, and the caption carries any altitude exaggeration.
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
                mv = coorbital.viz.globeMovie(demo,demoOpt);
    fprintf('globeMovie self-demo wrote %d frames to %s\n',mv.nFrame,mv.file);
    return;
end

%% Options, every one of them defaulted before anything is drawn:
    if nargin < 2 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
           outFile = vizOption(opts,'File',fullfile(tempdir,'globeMovie.mp4'));
            nFrame = vizOption(opts,'NFrame',90);
            fpsOut = vizOption(opts,'FrameRate',20);
            sizPix = vizOption(opts,'Size',[1280 720]);
            spinDg = vizOption(opts,'SpinDeg',30);
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
    assert(isscalar(altScale) && altScale > 0, ...
        'AltScale must be a positive scalar; %s given.',mat2str(altScale));
    assert(isempty(frmFcn) || isa(frmFcn,'function_handle'), ...
        'FrameFcn must be a function handle taking (hAx,kFrame), or [].');

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
%% recovering it clips nothing, at any spin angle:
    camzoom(hAx,sqrt(3));

%% A viewpoint looking down on the middle of the arc, as in globe3D: MATLAB's
%% view(az,el) points at longitude L and latitude B for az = L + 90, el = B:
              kMid = ceil(nS./2);
               az0 = rad2deg(lon(kMid)) + 90;
               el0 = rad2deg(lat(kMid));
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

%% One line per phase, each empty until the frame loop fills it, and each
%% extended backward to the junction so the growing track has no holes in it --
%% the same contract coorbital.viz.globe3D documents:
               col = lines(max(nPh,7));
            hTrack = gobjects(1,nPh);
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
        hTrack(kp) = line(hAx,NaN,NaN,NaN, ...
                          'Color',col(kp,:),'LineWidth',2.5, ...
                          'Tag','globeTrack','UserData',phList(kp));
    end

%% The launch point, on the surface and never exaggerated -- it marks ground,
%% and lifting it would put it somewhere the vehicle was not. There is no
%% impact marker: the vehicle has not got there yet, and drawing where it will
%% end up would give the ending away in frame one:
    line(hAx,rEkm.*cos(lat(1)).*cos(lon(1)), ...
             rEkm.*cos(lat(1)).*sin(lon(1)), ...
             rEkm.*sin(lat(1)), ...
         'LineStyle','none','Marker','o','MarkerSize',9,'LineWidth',1.5, ...
         'Color',[0.2 1 0.4],'Tag','launchMarker');

%% The vehicle. WHITE, and deliberately not one of the phase colours: lines()
%% reaches yellow at phase 3, and a yellow vehicle riding the end of a yellow
%% descent is a marker a reader cannot find:
              hVeh = line(hAx,NaN,NaN,NaN, ...
                          'LineStyle','none','Marker','o','MarkerSize',10, ...
                          'MarkerFaceColor',[1 1 1], ...
                          'MarkerEdgeColor',[0 0 0],'LineWidth',1.2, ...
                          'Tag','vehicleMarker');

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
    annotation(hFig,'textbox',[0.05 0.90 0.90 0.08], ...
               'String',titFul,'Tag','titleText', ...
               'Color',[1 1 1],'FontSize',14,'FontWeight','bold', ...
               'HorizontalAlignment','center','VerticalAlignment','top', ...
               'EdgeColor','none','FitBoxToText','off','Interpreter','none');
              hHud = annotation(hFig,'textbox',[0.015 0.76 0.30 0.14], ...
                          'String',{''},'Tag','hudText', ...
                          'Color',[1 1 1],'FontName','FixedWidth','FontSize',12, ...
                          'HorizontalAlignment','left','VerticalAlignment','top', ...
                          'EdgeColor','none','FitBoxToText','off','Interpreter','none');

%% The writer, opened last so that nothing can fail after a part-written file
%% exists on disk:
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

%% The vehicle and the readout, at the sample this frame shows:
            set(hVeh,'XData',xTrk(kEnd),'YData',yTrk(kEnd),'ZData',zTrk(kEnd));
               kNow = find(phList == phIdx(kEnd),1);
            set(hHud,'String',{sprintf('t   = %7.1f s',tS(kEnd)), ...
                               sprintf('alt = %7.2f km',altKm(kEnd)), ...
                               nmPhase{kNow}});

%% The camera, sweeping SpinDeg of azimuth across the whole movie:
            view(hAx,az0 + spinDg.*(kf - 1)./(nFrame - 1),el0);

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
