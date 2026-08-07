function test_viz()
%% Purpose:
%
%  Pin what the coorbital.viz figures SAY. A plot test that draws a figure and
%  checks that nothing threw is worth almost nothing: every mistake worth
%  catching -- plotting the wrong column, converting metres to kilometres
%  twice, labelling an axis without its unit, marking a launch point that is
%  not where the flight began -- produces a figure that draws perfectly and is
%  wrong. So this file reads the drawn objects back out and compares them
%  against the trajectory the figure claims to be showing.
%
%  Four things are asserted for every figure in the package:
%
%    1. The axes count is what the function promises.
%    2. Every axes carries a non-empty XLabel, YLabel and Title.
%    3. Every axis label names the UNIT that axis is displaying, which is the
%       only thing standing between a reader and a silent unit error.
%    4. The XData and YData of the primary line equal the trajectory arrays,
%       recomputed here from traj rather than read back from the function.
%
%  The fourth is the one that makes the file a test. The recomputations below
%  duplicate the arithmetic in the plotting routines ON PURPOSE, in the same
%  units, so that a wrong column, a stale array or a doubled conversion moves
%  one side and not the other.
%
%  THE TRAJECTORIES ARE SYNTHETIC, and deliberately so. The figures do not
%  integrate anything, so a real propagation would add seconds to the suite
%  and buy nothing; what a synthetic state buys instead is control over
%  phaseIdx, which is how the phase-segment contract below gets exercised at
%  all. The composition of viz with the real propagator is covered where it
%  belongs, in the entry-script tests.
%
%  Every figure is created with 'Visible','off' and closed again, and the
%  figure count is asserted back to its starting value at the end, so the
%  suite neither displays nor leaks.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
             nFig0 = numel(findall(groot,'Type','figure'));

%% The environment the derived panels are computed from, and three vehicles
%% that are genuinely different from one another. The differences are the
%% point: a load-factor panel that ignored the per-phase vehicle list and used
%% the single veh argument throughout would agree with the expectation below
%% only if the vehicles happened to match, so they are made not to:
         env.atmos = @coorbital.atmos.expAtmos;
          env.aero = @coorbital.aero.constLD;
          env.grav = @coorbital.grav.sphereGrav;
        env.omegaE = 0;
%%
%% The boost vehicle is the real coorbital.util.boosterDefaults struct rather
%% than a doctored copy of the glide vehicle, and that is deliberate: a
%% booster has NO mass field -- its mass is split into massDry and massProp,
%% and only the chain that assembles them knows the total. A load-factor panel
%% that reached for veh.mass unconditionally would throw on it, which is
%% exactly the defect this argument caught during construction:
            vehOne = coorbital.util.vehicleDefaults();
            vehBst = coorbital.util.boosterDefaults();
    assert(~isfield(vehBst,'mass'), ...
        ['boosterDefaults has grown a mass field; this test relied on its ' ...
         'absence to prove profilePlot does not need one.']);
          vehPhase = {vehBst,vehOne,vehOne};

%% A three-phase, seven-state trajectory, and a one-phase six-state one. Both
%% shapes ship in this library and they exercise different branches: the mass
%% column, and the per-phase vehicle list:
            trjCha = chainTraj(c.rE);
            trjOne = singleTraj(c.rE);

%% ---------------------------------------------------------------------
%% 1. groundTrack -- structure, labels and units
%% ---------------------------------------------------------------------
            target = [deg2rad(35.5); deg2rad(-77.0)];
            hFigGT = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Visible','off','Target',target, ...
                                'PhaseName',{{'boost','glide','descent'}}, ...
                                'Title','Ground track, 7663 km great-circle range'));
             axGT  = axesInOrder(hFigGT);
    assert(numel(axGT) == 1, ...
        'groundTrack drew %d axes; it promises exactly 1.',numel(axGT));
    assertLabelled(axGT,'groundTrack');
    assertUnit(axGT(1).XLabel.String,'deg','groundTrack x axis');
    assertUnit(axGT(1).YLabel.String,'deg','groundTrack y axis');
    assertUnit(axGT(1).XLabel.String,'longitude','groundTrack x axis');
    assertUnit(axGT(1).YLabel.String,'latitude','groundTrack y axis');

%% ---------------------------------------------------------------------
%% 2. groundTrack -- the plotted data IS the trajectory
%% ---------------------------------------------------------------------
%% One line per phase, each covering its own samples with the preceding
%% phase's last sample prepended, which is the contract the function's header
%% states and the thing that keeps the track from breaking at every junction.
%% Reconstructed here from phaseIdx alone, so a segment drawn over the wrong
%% index range fails even though the union of the segments still looks right:
            hTrack = lineByTag(axGT(1),'groundTrack');
            phList = unique(trjCha.phaseIdx(:)).';
    assert(numel(hTrack) == numel(phList), ...
        'groundTrack drew %d track segments for %d phases.', ...
        numel(hTrack),numel(phList));
    for kp = 1:numel(phList)
              hSeg = hTrack([hTrack.UserData] == phList(kp));
        assert(isscalar(hSeg), ...
            'expected exactly one track segment tagged for phase %d.',phList(kp));
              sel  = find(trjCha.phaseIdx == phList(kp));
        if kp > 1
              sel  = [find(trjCha.phaseIdx == phList(kp-1),1,'last'); sel];
        end
        assertSeries(hSeg.XData(:),rad2deg(trjCha.x(sel,2)), ...
            sprintf('groundTrack phase %d longitude (deg)',phList(kp)));
        assertSeries(hSeg.YData(:),rad2deg(trjCha.x(sel,3)), ...
            sprintf('groundTrack phase %d latitude (deg)',phList(kp)));
    end

%% ...and the markers, each against the sample it claims to mark:
    assertMarker(axGT(1),'launchMarker', ...
        rad2deg(trjCha.x(1,2)),rad2deg(trjCha.x(1,3)),'launch');
    assertMarker(axGT(1),'impactMarker', ...
        rad2deg(trjCha.x(end,2)),rad2deg(trjCha.x(end,3)),'impact');
    assertMarker(axGT(1),'targetMarker', ...
        rad2deg(target(2)),rad2deg(target(1)),'target');

%% The launch and impact markers must not be the same point, or the two
%% assertions above would both pass on a figure that marked one thing twice:
             dMark = abs(rad2deg(trjCha.x(1,2)) - rad2deg(trjCha.x(end,2)));
    assert(dMark > 1, ...
        'the synthetic track starts and ends within %.3f deg of longitude; the marker checks cannot discriminate.',dMark);

    close(hFigGT);

%% ---------------------------------------------------------------------
%% 3. groundTrack -- no target means NO target marker
%% ---------------------------------------------------------------------
%% The aim point is optional, and an optional marker that is drawn anyway --
%% at the origin, or at the impact point -- is worse than no marker at all,
%% because a reader would read the miss distance off it:
            hFigNT = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Visible','off'));
             axNT  = axesInOrder(hFigNT);
    assert(isempty(findobj(axNT(1),'Tag','targetMarker')), ...
        'groundTrack drew a target marker when no Target was supplied.');
    assert(~isempty(findobj(axNT(1),'Tag','launchMarker')) && ...
           ~isempty(findobj(axNT(1),'Tag','impactMarker')), ...
        'launch and impact must be marked whether or not a target is given.');
    close(hFigNT);

%% ---------------------------------------------------------------------
%% 4. profilePlot -- structure, labels and units
%% ---------------------------------------------------------------------
            hFigPP = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Visible','off','VehPhase',{vehPhase}));
             axPP  = axesInOrder(hFigPP);
    assert(numel(axPP) == 5, ...
        'profilePlot drew %d axes; it promises exactly 5.',numel(axPP));
    assertLabelled(axPP,'profilePlot');

%% Every panel shares one time axis, and every panel's own unit must appear in
%% its own label. The list is written out rather than looped over a pattern so
%% that a panel silently relabelled into another panel's unit still fails:
            unitsY = {'km','m/s','(-)','kPa','(g)'};
    for ka = 1:5
        assertUnit(axPP(ka).XLabel.String,'(s)', ...
            sprintf('profilePlot panel %d x axis',ka));
        assertUnit(axPP(ka).YLabel.String,unitsY{ka}, ...
            sprintf('profilePlot panel %d y axis',ka));
    end

%% ---------------------------------------------------------------------
%% 5. profilePlot -- every panel's data, recomputed from traj
%% ---------------------------------------------------------------------
    [hExp,vExp,machExp,qExp,nExp] = profileChannels(trjCha,vehPhase,c);
             wantY = {hExp,vExp,machExp,qExp,nExp};
             wantN = {'altitude (km)','speed (m/s)','Mach number', ...
                      'dynamic pressure (kPa)','aerodynamic load factor (g)'};
    for ka = 1:5
               hLn = lineByTag(axPP(ka),'profile');
        assert(isscalar(hLn), ...
            'profilePlot panel %d holds %d primary lines; expected 1.',ka,numel(hLn));
        assertSeries(hLn.XData(:),trjCha.t,sprintf('profilePlot panel %d time (s)',ka));
        assertSeries(hLn.YData(:),wantY{ka},sprintf('profilePlot %s',wantN{ka}));
    end

%% The five channels must be five DIFFERENT things. Without this, a function
%% that plotted altitude into all five panels could still pass every check
%% above if the expectations were ever built from the drawn data:
    for ka = 1:5
        for kb = ka+1:5
            assert(max(abs(wantY{ka} - wantY{kb})) > 1e-6, ...
                'panels %d and %d carry the same series; the checks cannot discriminate.',ka,kb);
        end
    end

%% Phase boundaries are drawn on every panel, one fewer than the phases:
    for ka = 1:5
              hBnd = findobj(axPP(ka),'Tag','phaseBoundary');
        assert(numel(hBnd) == numel(phList) - 1, ...
            'profilePlot panel %d drew %d phase boundaries for %d phases.', ...
            ka,numel(hBnd),numel(phList));
    end
    close(hFigPP);

%% ---------------------------------------------------------------------
%% 6. profilePlot -- the per-phase vehicle list is actually read
%% ---------------------------------------------------------------------
%% Run the same trajectory again with the single-vehicle default. The load
%% factor must MOVE, because the boost phase is then flown on the wrong
%% airframe. If it does not, VehPhase is being ignored and section 5 proved
%% nothing about it:
            hFigSV = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Visible','off'));
             axSV  = axesInOrder(hFigSV);
             nSing = lineByTag(axSV(5),'profile');
    assert(max(abs(nSing.YData(:) - nExp)) > 1e-3, ...
        ['the load factor is unchanged when the per-phase vehicle list is ' ...
         'withheld; profilePlot is ignoring VehPhase.']);
    close(hFigSV);

%% ---------------------------------------------------------------------
%% 7. profilePlot -- a six-state, single-phase trajectory
%% ---------------------------------------------------------------------
%% No mass column and no staging: the mass must then come from the vehicle,
%% and there must be no phase boundary to draw:
            hFig6  = coorbital.viz.profilePlot(trjOne,vehOne,env, ...
                         struct('Visible','off'));
             ax6   = axesInOrder(hFig6);
    assert(numel(ax6) == 5,'profilePlot drew %d axes on a six-state run.',numel(ax6));
    [h6,v6,m6,q6,n6] = profileChannels(trjOne,{vehOne},c);
            want6  = {h6,v6,m6,q6,n6};
    for ka = 1:5
               hLn = lineByTag(ax6(ka),'profile');
        assertSeries(hLn.YData(:),want6{ka}, ...
            sprintf('six-state profilePlot panel %d',ka));
    end
    assert(isempty(findobj(ax6(1),'Tag','phaseBoundary')), ...
        'a single-phase run has no phase boundary to draw.');
    close(hFig6);

%% ---------------------------------------------------------------------
%% 8. globe3D -- structure, labels, units and the arc
%% ---------------------------------------------------------------------
             altSc = 30;
            hFigG3 = coorbital.viz.globe3D(trjCha,vehOne,env, ...
                         struct('Visible','off','AltScale',altSc, ...
                                'Target',target));
             axG3  = axesInOrder(hFigG3);
    assert(numel(axG3) == 1,'globe3D drew %d axes; it promises exactly 1.',numel(axG3));
    assertLabelled(axG3,'globe3D');
    assertUnit(axG3(1).XLabel.String,'km','globe3D x axis');
    assertUnit(axG3(1).YLabel.String,'km','globe3D y axis');
    assertUnit(axG3(1).ZLabel.String,'km','globe3D z axis');

%% A non-unity exaggeration MUST be stated in the title. A vertical scale
%% thirty times the horizontal one, unannounced, is a lie told with a plot:
            titG3  = labelText(axG3(1).Title.String);
    assert(contains(titG3,'exaggerat') && contains(titG3,sprintf('%g',altSc)), ...
        ['globe3D was given AltScale = %g and its title does not say so. ' ...
         'Title was "%s".'],altSc,titG3);

%% The planet is drawn, once, by the shared renderer:
            hEarth = findobj(axG3(1),'Tag','earthSurface');
    assert(isscalar(hEarth),'globe3D drew %d Earth surfaces; expected 1.',numel(hEarth));
    assert(abs(max(hEarth.XData(:)) - c.rE./1000) < 1e-6, ...
        'the Earth surface has radius %.3f km against rE = %.3f km.', ...
        max(hEarth.XData(:)),c.rE./1000);

%% The arc, phase by phase, against the spherical-to-Cartesian conversion done
%% here independently. This is where an exaggeration applied to the sphere
%% instead of to the altitude, or a latitude/longitude transposition in the
%% conversion, is caught:
            hArc   = lineByTag(axG3(1),'globeTrack');
    assert(numel(hArc) == numel(phList), ...
        'globe3D drew %d arc segments for %d phases.',numel(hArc),numel(phList));
              rPlt = c.rE + altSc.*(trjCha.x(:,1) - c.rE);
              lonC = trjCha.x(:,2);
              latC = trjCha.x(:,3);
             xWant = rPlt.*cos(latC).*cos(lonC)./1000;
             yWant = rPlt.*cos(latC).*sin(lonC)./1000;
             zWant = rPlt.*sin(latC)./1000;
    for kp = 1:numel(phList)
              hSeg = hArc([hArc.UserData] == phList(kp));
        assert(isscalar(hSeg),'expected one globe arc segment for phase %d.',phList(kp));
              sel  = find(trjCha.phaseIdx == phList(kp));
        if kp > 1
              sel  = [find(trjCha.phaseIdx == phList(kp-1),1,'last'); sel];
        end
        assertSeries(hSeg.XData(:),xWant(sel),sprintf('globe3D phase %d x (km)',phList(kp)));
        assertSeries(hSeg.YData(:),yWant(sel),sprintf('globe3D phase %d y (km)',phList(kp)));
        assertSeries(hSeg.ZData(:),zWant(sel),sprintf('globe3D phase %d z (km)',phList(kp)));
    end

%% The arc really is lifted off the surface by the exaggeration, so the check
%% above is discriminating rather than comparing two things that coincide:
            liftKm = max(rPlt - c.rE)./1000;
    assert(liftKm > 100, ...
        'the exaggerated arc stands only %.3f km off the sphere.',liftKm);

%% ...and the surface markers are NOT lifted. They mark ground positions, and
%% an exaggerated launch marker sits somewhere the vehicle never was:
    assertSurfaceMarker(axG3(1),'launchMarker',c.rE,latC(1)  ,lonC(1)  ,'launch');
    assertSurfaceMarker(axG3(1),'impactMarker',c.rE,latC(end),lonC(end),'impact');
    assertSurfaceMarker(axG3(1),'targetMarker',c.rE,target(1),target(2),'target');
    close(hFigG3);

%% A unity exaggeration is the default and must NOT be announced, or every
%% honest figure would carry a caveat that does not apply to it:
            hFigG1 = coorbital.viz.globe3D(trjCha,vehOne,env,struct('Visible','off'));
             axG1  = axesInOrder(hFigG1);
    assert(~contains(labelText(axG1(1).Title.String),'exaggerat'), ...
        'globe3D announced an exaggeration it did not apply.');
    close(hFigG1);

%% ---------------------------------------------------------------------
%% 9. 'Parent' composes instead of creating
%% ---------------------------------------------------------------------
%% The whole point of the option: handed somewhere to draw, none of these may
%% open a figure of their own. Counted, not assumed:
             hHost = figure('Visible','off');
              axG  = axes('Parent',hHost);
           nBefore = numel(findall(groot,'Type','figure'));
              fGT  = coorbital.viz.groundTrack(trjCha,vehOne,env, ...
                         struct('Parent',axG));
    assert(numel(findall(groot,'Type','figure')) == nBefore, ...
        'groundTrack opened a figure although it was handed an axes.');
    assert(isequal(fGT,hHost), ...
        'groundTrack returned a figure other than the one its Parent lives in.');
    assert(~isempty(lineByTag(axG,'groundTrack')), ...
        'groundTrack drew nothing into the axes it was given.');

%% The globe goes into its own axes in the SAME figure -- a lat/lon grid and a
%% sphere sharing one set of limits would be unreadable, and the point being
%% made is only that no new figure appears:
              axG2 = axes('Parent',hHost);
              fG3  = coorbital.viz.globe3D(trjCha,vehOne,env, ...
                         struct('Parent',axG2));
    assert(~isempty(lineByTag(axG2,'globeTrack')), ...
        'globe3D drew nothing into the axes it was given.');
    assert(numel(findall(groot,'Type','figure')) == nBefore, ...
        'globe3D opened a figure although it was handed an axes.');
    assert(isequal(fG3,hHost), ...
        'globe3D returned a figure other than the one its Parent lives in.');

%% profilePlot needs five axes, so it is handed five:
            hHost5 = figure('Visible','off');
             ax5   = gobjects(1,5);
    for ka = 1:5
           ax5(ka) = axes('Parent',hHost5, ...
                          'OuterPosition',[0 (5-ka)./5 1 1./5]);
    end
          nBefore5 = numel(findall(groot,'Type','figure'));
            fPP    = coorbital.viz.profilePlot(trjCha,vehOne,env, ...
                         struct('Parent',ax5,'VehPhase',{vehPhase}));
    assert(numel(findall(groot,'Type','figure')) == nBefore5, ...
        'profilePlot opened a figure although it was handed five axes.');
    assert(isequal(fPP,hHost5), ...
        'profilePlot returned a figure other than the one its Parent lives in.');
    for ka = 1:5
        assert(isscalar(lineByTag(ax5(ka),'profile')), ...
            'profilePlot did not draw into supplied axes %d.',ka);
    end

%% ...and the wrong number of axes must be refused rather than half-filled:
             threw = false;
    try
        coorbital.viz.profilePlot(trjCha,vehOne,env,struct('Parent',ax5(1:3)));
    catch
             threw = true;
    end
    assert(threw,'profilePlot accepted three axes for a five-panel figure.');

    close(hHost);
    close(hHost5);

%% ---------------------------------------------------------------------
%% 10. No leaked figures
%% ---------------------------------------------------------------------
             nFig1 = numel(findall(groot,'Type','figure'));
    assert(nFig1 == nFig0, ...
        'test_viz leaked %d figure(s); the suite must end with the %d it started with.', ...
        nFig1 - nFig0,nFig0);
end

function traj = chainTraj(rE)
%% Purpose:
%
%  A synthetic three-phase, seven-state trajectory: a climb, a long shallow
%  glide and a dive, with the mass falling through phase 1 and constant
%  afterwards. It is not a solution of anything, and it does not need to be --
%  what a plotting test needs is a state history whose every column is
%  distinct and whose phase structure is known exactly.
%
%% Inputs:
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%% Outputs:
%
%  traj             Struct                      Fields t [N x 1] (s),
%                                               x [N x 7], u [N x 2] (rad),
%                                               phaseIdx [N x 1], junction
%                                               [2 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                n1 = 20;
                n2 = 40;
                n3 = 15;
                nS = n1 + n2 + n3;
            traj.t = linspace(0,2100,nS)';

%% Altitude: up to 47 km through boost, a decaying skip through the glide,
%% then down to the deck:
               hKm = [linspace(0.2,47,n1), ...
                      40 + 8.*exp(-linspace(0,2.2,n2)).*cos(linspace(0,7,n2)), ...
                      linspace(15,0.1,n3)]';
               spd = [linspace(10,4600,n1), ...
                      linspace(4600,900,n2), ...
                      linspace(900,400,n3)]';
               lon = deg2rad(linspace(-155,-77.6,nS))';
               lat = deg2rad(linspace(20,34,nS))';
               gam = deg2rad([linspace(89,2,n1), ...
                              2.*cos(linspace(0,7,n2)), ...
                              linspace(-5,-32,n3)])';
               psi = deg2rad(linspace(60,80,nS))';
               mas = [linspace(9400,4300,n1), 900.*ones(1,n2 + n3)]';

            traj.x = [rE + 1000.*hKm, lon, lat, spd, gam, psi, mas];
            traj.u = zeros(nS,2);
     traj.phaseIdx = [ones(n1,1); 2.*ones(n2,1); 3.*ones(n3,1)];
     traj.junction = repmat(struct('t',[],'x',[]),2,1);
    for kj = 1:2
        traj.junction(kj).t = traj.t(find(traj.phaseIdx == kj,1,'last'));
        traj.junction(kj).x = traj.x(find(traj.phaseIdx == kj,1,'last'),:)';
    end
end

function traj = singleTraj(rE)
%% Purpose:
%
%  A synthetic one-phase, six-state trajectory -- the shape HGV/run_glide
%  returns. Exists so the mass-column and per-phase-vehicle branches of
%  coorbital.viz.profilePlot are exercised on the side where neither is
%  present.
%
%% Inputs:
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%% Outputs:
%
%  traj             Struct                      Fields t [N x 1] (s),
%                                               x [N x 6], u [N x 2] (rad),
%                                               phaseIdx [N x 1] all ones,
%                                               junction [0 x 1] struct
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = 50;
            traj.t = linspace(0,2073,nS)';
               hKm = linspace(60,5,nS)';
               spd = linspace(6000,321,nS)';
               lon = deg2rad(linspace(0,62.76,nS))';
               lat = zeros(nS,1);
               gam = deg2rad(linspace(-1,-12.6,nS))';
               psi = deg2rad(90).*ones(nS,1);
            traj.x = [rE + 1000.*hKm, lon, lat, spd, gam, psi];
            traj.u = zeros(nS,2);
     traj.phaseIdx = ones(nS,1);
     traj.junction = repmat(struct('t',[],'x',[]),0,1);
end

function [hKm,V,mach,qKPa,nAero] = profileChannels(traj,vehPhase,c)
%% Purpose:
%
%  Recompute the five channels coorbital.viz.profilePlot draws, in the units
%  it draws them in, from the trajectory and the per-phase vehicles. The
%  duplication of the plotting routine's arithmetic is the whole mechanism:
%  the two are compared, so a wrong column or a doubled conversion moves one
%  and not the other.
%
%% Inputs:
%
%  traj             Struct                      Trajectory to recompute from
%
%  vehPhase         Cell [1 x P]                One vehicle per phase
%
%  c                Struct                      Constants from missileConst;
%                                               uses rE (m) and g0 (m/s^2)
%
%% Outputs:
%
%  hKm              [N x 1]                     Altitude (km)
%
%  V                [N x 1]                     Speed (m/s)
%
%  mach             [N x 1]                     Mach number (-)
%
%  qKPa             [N x 1]                     Dynamic pressure (kPa)
%
%  nAero            [N x 1]                     Aerodynamic load factor (g)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = numel(traj.t);
                hM = traj.x(:,1) - c.rE;
                 V = traj.x(:,4);
    [rho,~,~,aSnd] = coorbital.atmos.expAtmos(hM);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;
               hKm = hM./1000;
              qKPa = qbar./1000;
            phList = unique(traj.phaseIdx(:)).';

    if size(traj.x,2) >= 7
             massV = traj.x(:,7);
    else
             massV = zeros(nS,1);
        for kp = 1:numel(phList)
               sel = traj.phaseIdx == phList(kp);
        massV(sel) = vehPhase{kp}.mass;
        end
    end

             aLift = zeros(nS,1);
             aDrag = zeros(nS,1);
    for ks = 1:nS
              kPh  = find(phList == traj.phaseIdx(ks),1);
              vehK = vehPhase{kPh};
         [CLk,CDk] = coorbital.aero.constLD(traj.u(ks,1),mach(ks),vehK);
         aLift(ks) = qbar(ks).*vehK.Sref.*CLk./massV(ks);
         aDrag(ks) = qbar(ks).*vehK.Sref.*CDk./massV(ks);
    end
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;
end

function hAx = axesInOrder(hFig)
%% Purpose:
%
%  Return a figure's axes ordered top to bottom on the page. findobj hands
%  them back in reverse creation order, which is an implementation detail of
%  the graphics system rather than anything a test should depend on; the
%  panels of a stacked profile figure have to be identified by WHERE they are.
%  Legends and colorbars are not axes and do not appear here.
%
%% Inputs:
%
%  hFig             [1 x 1] figure              Figure to search
%
%% Outputs:
%
%  hAx              [1 x M] axes                Axes, topmost first
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hAx = findobj(hFig,'Type','axes');
    assert(~isempty(hAx),'the figure carries no axes at all.');
              yTop = zeros(1,numel(hAx));
    for ka = 1:numel(hAx)
              oPos = hAx(ka).OuterPosition;
          yTop(ka) = oPos(2);
    end
            [~,ord] = sort(yTop,'descend');
               hAx = reshape(hAx(ord),1,[]);
end

function hLn = lineByTag(hAx,tagTxt)
%% Purpose:
%
%  All line objects in one axes carrying a given tag, as a row.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag to match
%
%% Outputs:
%
%  hLn              [1 x M] line                Matching lines, possibly empty
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hLn = findobj(hAx,'Type','line','Tag',tagTxt);
               hLn = reshape(hLn,1,[]);
end

function txt = labelText(str)
%% Purpose:
%
%  Flatten an axis label or title into one char row. MATLAB allows a label to
%  be char, string or a cellstr of stacked lines, and a test that assumed one
%  of the three would fail for the wrong reason.
%
%% Inputs:
%
%  str              Char / string / cellstr     Label as read from the object
%
%% Outputs:
%
%  txt              Char [1 x n]                One-line rendering
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if iscell(str) || (isstring(str) && ~isscalar(str))
               txt = char(strjoin(string(str),' '));
    else
               txt = char(str);
    end
end

function assertLabelled(hAx,what)
%% Purpose:
%
%  Every axes in a figure must carry a non-empty XLabel, YLabel and Title. An
%  unlabelled axis is not a stylistic lapse; it is a plot whose reader cannot
%  tell what is being shown or in what unit.
%
%% Inputs:
%
%  hAx              [1 x M] axes                Axes to check
%
%  what             Char [1 x n]                Name of the function under
%                                               test, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    for ka = 1:numel(hAx)
              want = {'XLabel','YLabel','Title'};
        for kw = 1:numel(want)
               txt = labelText(hAx(ka).(want{kw}).String);
            assert(~isempty(strtrim(txt)), ...
                '%s axes %d has an empty %s.',what,ka,want{kw});
        end
    end
end

function assertUnit(str,tok,what)
%% Purpose:
%
%  Assert that an axis label contains a required token -- normally the unit
%  the axis is displaying. This is the check that catches a panel showing
%  kilometres under a label that says only "altitude", which is the defect a
%  reader has no way of detecting from the figure.
%
%% Inputs:
%
%  str              Char / cellstr              Label as read from the object
%
%  tok              Char [1 x n]                Token that must appear
%
%  what             Char [1 x m]                Name of the axis, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               txt = labelText(str);
    assert(contains(txt,tok), ...
        '%s is labelled "%s", which does not state "%s".',what,txt,tok);
end

function assertMarker(hAx,tagTxt,xWant,yWant,what)
%% Purpose:
%
%  Assert that a tagged single-point marker exists in an axes and sits at the
%  expected coordinates. Existence alone is not enough: a launch marker drawn
%  at the impact point, or at the origin, is present and wrong.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag of the marker
%
%  xWant            [1 x 1]                     Expected x coordinate, in the
%                                               axis's displayed unit
%
%  yWant            [1 x 1]                     Expected y coordinate
%
%  what             Char [1 x m]                Name of the marker, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hMk = findobj(hAx,'Tag',tagTxt);
    assert(isscalar(hMk),'expected exactly one %s marker, found %d.',what,numel(hMk));
    assert(isscalar(hMk.XData) && isscalar(hMk.YData), ...
        'the %s marker holds %d points; a marker is one point.',what,numel(hMk.XData));
    assert(abs(hMk.XData - xWant) < 1e-9, ...
        'the %s marker is at x = %.9f against the trajectory''s %.9f', ...
        what,hMk.XData,xWant);
    assert(abs(hMk.YData - yWant) < 1e-9, ...
        'the %s marker is at y = %.9f against the trajectory''s %.9f', ...
        what,hMk.YData,yWant);
end

function assertSurfaceMarker(hAx,tagTxt,rE,latRad,lonRad,what)
%% Purpose:
%
%  Assert that a tagged marker in a three-dimensional globe figure sits on the
%  sphere at the given latitude and longitude, in kilometres. Written against
%  the spherical-to-Cartesian conversion done here rather than read back from
%  the figure, so a transposition of latitude and longitude fails.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to search
%
%  tagTxt           Char [1 x n]                Tag of the marker
%
%  rE               [1 x 1]                     Sphere radius (m)
%
%  latRad           [1 x 1]                     Expected latitude (rad)
%
%  lonRad           [1 x 1]                     Expected longitude (rad)
%
%  what             Char [1 x m]                Name of the marker, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hMk = findobj(hAx,'Tag',tagTxt);
    assert(isscalar(hMk),'expected exactly one %s marker, found %d.',what,numel(hMk));
             xWant = rE.*cos(latRad).*cos(lonRad)./1000;
             yWant = rE.*cos(latRad).*sin(lonRad)./1000;
             zWant = rE.*sin(latRad)./1000;
               got = [hMk.XData hMk.YData hMk.ZData];
              want = [xWant yWant zWant];
               err = max(abs(got - want));
    assert(err < 1e-6, ...
        ['the %s marker is at [%.4f %.4f %.4f] km against [%.4f %.4f %.4f] km ' ...
         'from its latitude and longitude; %.3e km apart.'], ...
        what,got(1),got(2),got(3),want(1),want(2),want(3),err);
end

function assertSeries(got,want,what)
%% Purpose:
%
%  Assert that a series read back out of a figure equals the series recomputed
%  from the trajectory. The budget is relative to the largest magnitude in the
%  reference, so it is meaningful for a channel that passes through zero, and
%  it is tight enough that no real defect can hide under it: the loosest thing
%  it has to tolerate is a reassociation of the same arithmetic, and the
%  smallest thing it has to catch -- a doubled kilometre conversion -- is off
%  by three orders of magnitude.
%
%% Inputs:
%
%  got              [N x 1]                     Series read from the figure
%
%  want             [N x 1]                     Series recomputed from traj
%
%  what             Char [1 x n]                Name of the series, for the
%                                               message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(numel(got) == numel(want), ...
        '%s: the figure holds %d points against %d in the trajectory.', ...
        what,numel(got),numel(want));
             scale = max(abs(want));
    if scale == 0
             scale = 1;
    end
               err = max(abs(got(:) - want(:)))./scale;
    assert(err < 1e-12, ...
        '%s: the plotted series differs from the trajectory by %.3e relative.', ...
        what,err);
end
