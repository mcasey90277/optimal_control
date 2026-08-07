function hFig = globe3D(traj,veh,env,opts)
%% Purpose:
%
%  Draw the reference sphere in three dimensions with the trajectory arc over
%  it, as a single still figure. The companion to coorbital.viz.groundTrack:
%  the ground track answers where the vehicle went, this answers what the
%  flight looks like over a round planet, which is the picture a flat
%  lat/lon grid cannot give at intercontinental range.
%
%  Everything is drawn in kilometres and all three axis labels say so. The
%  planet itself is rendered by the package-private earthSurface, which
%  coorbital.viz.globeMovie renders with as well, so the still and the movie
%  cannot end up showing two different Earths.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun. Reads
%                                               x(:,1) radius (m), x(:,2)
%                                               longitude (rad), x(:,3)
%                                               latitude (rad) and phaseIdx
%                                               [N x 1]
%
%  veh              Struct                      Vehicle parameters. NOT READ.
%                                               Present so that every function
%                                               in coorbital.viz takes the same
%                                               (traj,veh,env,opts) arguments
%
%  env              Struct                      Environment model handles. NOT
%                                               READ. The arc is the state's
%                                               own geometry; no model enters
%
%  opts             Struct, optional            All fields optional:
%                                               Parent    axes to draw into, or
%                                                         a figure/panel to
%                                                         create one in; [] or
%                                                         absent makes a new
%                                                         figure
%                                               Target    [2 x 1] aim point
%                                                         [lat; lon] (RAD),
%                                                         marked on the surface
%                                                         when given
%                                               AltScale  [1 x 1] altitude
%                                                         exaggeration (-),
%                                                         default 1. See Notes
%                                               Title     char, axes title. A
%                                                         non-unity AltScale is
%                                                         appended to whatever
%                                                         is given
%                                               Visible   'on' (default) or
%                                                         'off'; applies only
%                                                         to a figure this
%                                                         function creates
%
%% Outputs:
%
%  hFig             [1 x 1] figure              The figure the globe was drawn
%                                               in: the one created here, or
%                                               the one the Parent lives in
%
%% Notes:
%
%  ALTITUDE EXAGGERATION, AND WHY IT DEFAULTS TO OFF. A 60 km glide over a
%  6378 km sphere stands one part in 106 off the surface. Drawn true to
%  scale the arc lies on the planet and tells a reader nothing about the
%  vertical profile, which is why almost every published globe plot
%  exaggerates it. This one will too, but only when asked: AltScale multiplies
%  the ALTITUDE ABOVE THE SPHERE and leaves the sphere alone, so the plotted
%  radius is rE + AltScale*(r - rE). The default of 1 draws the truth. Any
%  other value is stated in the axis title, because a picture whose vertical
%  scale is a hundred times its horizontal one and does not say so is a lie
%  told with a plot.
%
%  The exaggeration moves the arc only. It does not move the launch or impact
%  markers, which sit on the surface where the flight began and ended, nor the
%  target marker.
%
%  PLANET-FIXED, NOT INERTIAL. The state's longitude is measured in whatever
%  frame the run was flown in -- planet-fixed when env.omegaE is nonzero, and
%  a non-rotating frame otherwise, the two coinciding for the omegaE = 0 runs
%  this library ships. The sphere is drawn in that same frame, so the arc sits
%  over the ground it actually flew over either way; what changes is what the
%  fixed x-axis means.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A synthetic two-phase arc, so the demo costs no integration:
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
          demoOpt  = struct('AltScale',20, ...
                            'Title','coorbital.viz.globe3D self-demo');
              hFig = coorbital.viz.globe3D(demo, ...
                         coorbital.util.vehicleDefaults(),struct(),demoOpt);
    return;
end

%% Options, every one of them defaulted before anything is drawn:
    if nargin < 4 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of options.');
            parent = vizOption(opts,'Parent',[]);
            target = vizOption(opts,'Target',[]);
          altScale = vizOption(opts,'AltScale',1);
            titTxt = vizOption(opts,'Title','Trajectory over the Earth');
            visTxt = vizOption(opts,'Visible','on');
    assert(isscalar(altScale) && altScale > 0, ...
        'AltScale must be a positive scalar; %s given.',mat2str(altScale));

                 c = coorbital.util.missileConst();
             phIdx = traj.phaseIdx;
            phList = unique(phIdx(:)).';
              nPh  = numel(phList);

%% The arc, in kilometres. The exaggeration multiplies the altitude above the
%% sphere and never the sphere, so the planet stays the right size:
              rPlt = c.rE + altScale.*(traj.x(:,1) - c.rE);
               lon = traj.x(:,2);
               lat = traj.x(:,3);
             xTrk  = rPlt.*cos(lat).*cos(lon)./1000;
             yTrk  = rPlt.*cos(lat).*sin(lon)./1000;
             zTrk  = rPlt.*sin(lat)./1000;

%% Axes to draw into, and the only place a figure can be created:
        [hFig,hAx] = vizParent(parent,1,'Globe',visTxt);
           wasHeld = ishold(hAx);
    hold(hAx,'on');

%% The planet, from the one routine in this package allowed to draw it:
    earthSurface(hAx,36);

%% One line per phase, extended backward to the junction so the arc has no
%% holes in it -- the same contract coorbital.viz.groundTrack documents:
               col = lines(max(nPh,7));
    for kp = 1:nPh
               sel = find(phIdx == phList(kp));
        if kp > 1
               sel = [find(phIdx == phList(kp-1),1,'last'); sel];
        end
        line(hAx,xTrk(sel),yTrk(sel),zTrk(sel), ...
             'Color',col(kp,:),'LineWidth',2, ...
             'Tag','globeTrack','UserData',phList(kp));
    end

%% Launch and impact, on the surface. Not exaggerated: they mark ground
%% positions, and lifting them would put them somewhere the vehicle was not:
    surfaceMarker(hAx,c.rE,lat(1)  ,lon(1)  ,'o',[0 0.45 0],'launchMarker');
    surfaceMarker(hAx,c.rE,lat(end),lon(end),'s',[0.75 0 0],'impactMarker');
    if ~isempty(target)
        assert(numel(target) == 2, ...
            'Target must be [lat; lon] in radians, %d element(s) given.',numel(target));
        surfaceMarker(hAx,c.rE,target(1),target(2),'p',[0 0 0.8],'targetMarker');
    end

%% A viewpoint looking straight down on the middle of the arc, so the flight
%% faces the reader whatever hemisphere it is in. MATLAB's view(az,el) places
%% the camera along [sin(az)cos(el), -cos(az)cos(el), sin(el)], which points
%% at longitude L and latitude B for az = L + 90 and el = B:
              kMid = ceil(numel(lon)./2);
    view(hAx,rad2deg(lon(kMid)) + 90,rad2deg(lat(kMid)));
    axis(hAx,'equal');
    grid(hAx,'on');

%% Labels. All three carry the unit the coordinates were converted into, and
%% the title states any exaggeration rather than leaving it to be discovered:
    xlabel(hAx,'x (km)');
    ylabel(hAx,'y (km)');
    zlabel(hAx,'z (km)');
    if altScale == 1
        title(hAx,titTxt);
    else
        title(hAx,sprintf('%s (altitude exaggerated %gx)',titTxt,altScale));
    end

    if ~wasHeld
        hold(hAx,'off');
    end
end

function hMk = surfaceMarker(hAx,rE,latRad,lonRad,mkr,col,tagTxt)
%% Purpose:
%
%  Place one marker on the surface of the sphere at a given latitude and
%  longitude. Factored out because the launch, impact and target markers
%  differ only in symbol, colour and tag, and writing the spherical-to-
%  Cartesian conversion three times is three chances to transpose it.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to draw into
%
%  rE               [1 x 1]                     Sphere radius (m)
%
%  latRad           [1 x 1]                     Latitude (rad)
%
%  lonRad           [1 x 1]                     Longitude (rad)
%
%  mkr              Char [1 x 1]                MATLAB marker symbol
%
%  col              [1 x 3]                     Marker colour, RGB in [0,1]
%
%  tagTxt           Char [1 x n]                Tag for the created line, so a
%                                               test can find it by name
%
%% Outputs:
%
%  hMk              [1 x 1] line                The created marker
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               xMk = rE.*cos(latRad).*cos(lonRad)./1000;
               yMk = rE.*cos(latRad).*sin(lonRad)./1000;
               zMk = rE.*sin(latRad)./1000;
               hMk = line(hAx,xMk,yMk,zMk, ...
                          'LineStyle','none','Marker',mkr,'MarkerSize',10, ...
                          'LineWidth',1.5,'Color',col,'Tag',tagTxt);
end
